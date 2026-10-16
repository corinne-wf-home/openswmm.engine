/**
 * @file RomAdditiveChannel.hpp
 * @brief PR H19: full-order linearized additive-uncertainty channel.
 *
 * One head-deviation field Psi per additive source, shared by all members:
 *
 *     S dPsi/dt = -J Psi + dq
 *
 * J is the flow-response Jacobian (per conduit dQ = Ju*Psi_up + Jd*Psi_dn, per
 * structure dQ = d*Psi_up, oriented by flow), so extra water is carried
 * downstream step by step, one-way where the flow response is one-way. Backward
 * Euler: (S/dt + J) Psi_new = (S/dt) Psi_old + dq. The matrix is a strictly
 * column-diagonally-dominant M-matrix (column sums = S/dt plus outfall losses:
 * volume conservation), so the step is unconditionally stable and Gauss-Seidel
 * converges. See docs/uncertainty/ADDITIVE_CHANNEL.md.
 */
#ifndef OPENSWMM_ENGINE_UNCERTAINTY_ROM_ADDITIVE_CHANNEL_HPP
#define OPENSWMM_ENGINE_UNCERTAINTY_ROM_ADDITIVE_CHANNEL_HPP

#include "RomFlowResponse.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace openswmm::uncertainty {

/// The backward-Euler system for one routing step, shared by every plane.
struct AdditiveSystem {
    int n = 0;
    double dt = 0.0;
    std::vector<double> diag;    ///< S/dt + outgoing conductance + downstream-end coupling
    std::vector<double> s_dt;    ///< S/dt
    std::vector<int> start, col; ///< neighbours per row (CSR)
    std::vector<double> w;       ///< x_i = (rhs_i + sum w * x_col) / diag_i, w >= 0
    std::vector<int> order;      ///< upstream-first sweep order
    std::vector<char> pinned;    ///< nodes held at zero (spilling)
    int last_sweeps = 0;

    /// Assemble from the oriented flow-response edges, storage S (> 0) and dt.
    void build(int n_nodes, const std::vector<FlowEdge>& edges, const double* S,
               double dt_s, const std::vector<char>* pinned_nodes = nullptr) {
        n = n_nodes;
        dt = dt_s;
        const auto un = static_cast<std::size_t>(std::max(n, 0));
        diag.assign(un, 0.0);
        s_dt.assign(un, 0.0);
        pinned.assign(un, 0);
        if (pinned_nodes)
            for (std::size_t i = 0; i < un && i < pinned_nodes->size(); ++i) pinned[i] = (*pinned_nodes)[i];
        for (std::size_t i = 0; i < un; ++i) {
            s_dt[i] = S[i] / dt;
            diag[i] = s_dt[i];
        }
        std::vector<std::vector<std::pair<int, double>>> nb(un);
        for (const auto& e : edges) {
            const bool up_ok = e.up >= 0 && e.up < n;
            const bool dn_ok = e.dn >= 0 && e.dn < n;
            if (!up_ok) continue;
            const auto uu = static_cast<std::size_t>(e.up);
            if (e.conduit) {
                diag[uu] += e.ju;                       // outflow grows with the upstream head
                if (dn_ok) {
                    const auto ud = static_cast<std::size_t>(e.dn);
                    nb[uu].emplace_back(e.dn, -e.jd);   // ...and falls as the downstream head rises
                    diag[ud] += -e.jd;                  // inflow to dn falls as its own head rises
                    nb[ud].emplace_back(e.up, e.ju);    // ...and grows with the upstream head
                }
            } else {
                diag[uu] += e.d;
                if (dn_ok) nb[static_cast<std::size_t>(e.dn)].emplace_back(e.up, e.d);
            }
        }
        start.assign(un + 1, 0);
        for (std::size_t i = 0; i < un; ++i) start[i + 1] = start[i] + static_cast<int>(nb[i].size());
        col.resize(static_cast<std::size_t>(start[un]));
        w.resize(static_cast<std::size_t>(start[un]));
        for (std::size_t i = 0; i < un; ++i)
            for (std::size_t k = 0; k < nb[i].size(); ++k) {
                col[static_cast<std::size_t>(start[i]) + k] = nb[i][k].first;
                w[static_cast<std::size_t>(start[i]) + k] = nb[i][k].second;
            }
        // Upstream-first order (Kahn on the flow graph); cycle members appended.
        std::vector<int> indeg(un, 0);
        std::vector<std::vector<int>> next(un);
        for (const auto& e : edges)
            if (e.up >= 0 && e.up < n && e.dn >= 0 && e.dn < n) {
                ++indeg[static_cast<std::size_t>(e.dn)];
                next[static_cast<std::size_t>(e.up)].push_back(e.dn);
            }
        order.clear();
        order.reserve(un);
        std::vector<char> seen(un, 0);
        for (int i = 0; i < n; ++i) if (indeg[static_cast<std::size_t>(i)] == 0) { order.push_back(i); seen[static_cast<std::size_t>(i)] = 1; }
        for (std::size_t h = 0; h < order.size(); ++h)
            for (int d : next[static_cast<std::size_t>(order[h])])
                if (--indeg[static_cast<std::size_t>(d)] == 0 && !seen[static_cast<std::size_t>(d)]) {
                    order.push_back(d);
                    seen[static_cast<std::size_t>(d)] = 1;
                }
        for (int i = 0; i < n; ++i) if (!seen[static_cast<std::size_t>(i)]) order.push_back(i);
    }

    /// One backward-Euler step of one plane: psi (in: old, out: new), dq in cfs.
    /// Returns the Gauss-Seidel sweeps used.
    int step(const double* dq, double* psi, int max_sweeps = 200, double tol = 1.0e-10) {
        const auto un = static_cast<std::size_t>(std::max(n, 0));
        rhs_.resize(un);
        for (std::size_t i = 0; i < un; ++i) rhs_[i] = s_dt[i] * psi[i] + dq[i];
        auto relax = [&](int i) -> double {
            const auto ui = static_cast<std::size_t>(i);
            if (pinned[ui]) { const double d = std::fabs(psi[ui]); psi[ui] = 0.0; return d; }
            double v = rhs_[ui];
            for (int k = start[ui]; k < start[ui + 1]; ++k)
                v += w[static_cast<std::size_t>(k)] * psi[static_cast<std::size_t>(col[static_cast<std::size_t>(k)])];
            v /= diag[ui];
            const double d = std::fabs(v - psi[ui]);
            psi[ui] = v;
            return d;
        };
        int sweeps = 0;
        for (; sweeps < max_sweeps; ) {
            double change = 0.0;
            if (sweeps % 2 == 0) { for (int i : order) change = std::max(change, relax(i)); }
            else { for (auto it = order.rbegin(); it != order.rend(); ++it) change = std::max(change, relax(*it)); }
            ++sweeps;
            double scale = 0.0;
            for (std::size_t i = 0; i < un; ++i) scale = std::max(scale, std::fabs(psi[i]));
            if (change <= tol * scale || change == 0.0) break;
        }
        for (std::size_t i = 0; i < un; ++i) if (!std::isfinite(psi[i])) psi[i] = 0.0;
        last_sweeps = sweeps;
        return sweeps;
    }

private:
    std::vector<double> rhs_;
};

}  // namespace openswmm::uncertainty

#endif  // OPENSWMM_ENGINE_UNCERTAINTY_ROM_ADDITIVE_CHANNEL_HPP
