/**
 * @file RomDirectionalOperator.hpp
 * @brief PR H14 — Froude-gated directional reduced operator for the 1D ROM.
 *
 * @details The 1D ROM's operator is the weighted graph Laplacian built from
 *          the Picard `dqdh` (one coefficient per conduit, applied to both
 *          endpoints), so it is symmetric by construction and carries no
 *          flow direction. In Saint-Venant a supercritical reach transmits
 *          nothing upstream; the symmetric operator transmits everything.
 *          H5b measured the consequence on its surcharged chain: a pool's
 *          10 ft deviation diffused upstream across a 5 % supercritical
 *          reach and left the free-surface nodes above it 7–15× over-wide
 *          (VALIDATION.md, "H5b" §4).
 *
 *          This header assembles the reduced k×k operator
 *
 *              M = Pᵀ · L_dir · P,
 *
 *          where L_dir is the same weighted, outfall-grounded Laplacian with
 *          one change: on each conduit whose Froude number exceeds the gate,
 *          the UPSTREAM node's coupling to the downstream node is removed
 *          (the entry L[up, down]), while the downstream node still sees the
 *          upstream one. The upstream row's diagonal is unchanged, so that
 *          node still relaxes at its own rate — it just no longer follows
 *          what happens below it. The gate ramps smoothly over
 *          [fr_lo, fr_hi] so a reach hovering near critical does not
 *          switch the operator on and off step to step:
 *
 *              L_dir[up, down] = −w · (1 − g(Fr)),   g = clamp((Fr − fr_lo)/(fr_hi − fr_lo), 0, 1).
 *
 *          A second, node-level signal catches what a mid-depth Froude number
 *          cannot: a steep reach whose DOWNSTREAM end is drowned by a pool
 *          reads as subcritical (the mid-depth is large), yet its upstream
 *          node sits far above the pool's water surface and is beyond any
 *          backwater reach. The "drop gate" opens when the water-surface drop
 *          across the conduit exceeds the upstream node's own depth:
 *
 *              r = (h_up − h_dn) / d_up,   g_drop = clamp((r − drop_lo)/(drop_hi − drop_lo), 0, 1),
 *
 *          fully coupled below one depth of drop, fully one-way above two.
 *          The applied gate is max(g_Fr, g_drop). Measured on H5b's chain
 *          the drop gate barely moved J3 (above the pool, reach C3 drowned at
 *          its lower end: 1.20 → 1.07 ft against MC 0.108) and narrowed the
 *          subcritical soft-rain fixtures (CL-1e 0.827 → 0.762, red), so it
 *          ships OFF; J3's residual is not relayed through the operator at
 *          all (see VALIDATION.md "H14" for the truncation measurement).
 *
 *          With every gate at 0, M equals Pᵀ L P, which is diag(λ) up to the
 *          eigensolver's Ritz residual — the diagonal path the ROM has always
 *          integrated. M is not symmetric once a gate opens, so the ROM
 *          integrates it with the dense matrix exponential
 *          (RomDensePropagator.hpp), exactly as the 2D ROM integrates its
 *          anisotropic/advective operator (W3) — the 1D operator is small
 *          (k ≈ 20) so the per-member exponential is a few µs.
 *
 *          Pure functions, no engine dependencies (RomDiagTrust.hpp pattern):
 *          everything is testable with hand-built arrays.
 *
 * @ingroup engine_uncertainty
 */

#ifndef OPENSWMM_ENGINE_UNCERTAINTY_ROM_DIRECTIONAL_OPERATOR_HPP
#define OPENSWMM_ENGINE_UNCERTAINTY_ROM_DIRECTIONAL_OPERATOR_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace openswmm::uncertainty {

/// Dials for the Froude gate. Not parser-exposed (same policy as H5's
/// SurchargeAttenuationConfig and H11's PhaseConfig).
struct DirectionalOperatorConfig {
    /// Off switch: false ⇒ the engine never installs a directional operator
    /// and the ROM stays on its pre-H14 diagonal path, bit-identical.
    bool   enabled = true;
    /// Gate ramp: fully coupled at Fr ≤ fr_lo, fully one-way at Fr ≥ fr_hi.
    double fr_lo   = 0.8;
    double fr_hi   = 1.2;
    /// Drop gate (node level): ramp on (h_up − h_dn)/d_up. OFF BY DEFAULT —
    /// measured and rejected 2026-10-07: it moved H5b's drowned-outlet J3
    /// only from 1.20 to 1.07 ft (MC 0.108) while gating steep-but-subcritical
    /// reaches on the soft-rain fixtures, which took CL-1e from 0.827 to 0.762
    /// member coverage (red) and SR-5 from 0.881 to 0.869. Kept as a dial so
    /// the measurement is reproducible (H14_DROP=1 in the MC harness), not as
    /// a recommendation. Disabled when drop_hi ≤ drop_lo or use_drop is false.
    bool   use_drop = false;
    double drop_lo  = 1.0;
    double drop_hi  = 2.0;
};

/// Fraction of the upstream→downstream coupling REMOVED on a conduit whose
/// water-surface drop is `drop_ratio` upstream depths: 0 below drop_lo, 1
/// above drop_hi, linear between; 0 when the drop gate is disabled.
inline double dropGate(double drop_ratio, const DirectionalOperatorConfig& cfg) {
    if (!cfg.use_drop || !(drop_ratio > 0.0)) return 0.0;
    const double band = cfg.drop_hi - cfg.drop_lo;
    if (band <= 0.0) return 0.0;
    return std::clamp((drop_ratio - cfg.drop_lo) / band, 0.0, 1.0);
}

/// Fraction of the upstream→downstream coupling REMOVED on a conduit at
/// Froude number `fr`: 0 below fr_lo, 1 above fr_hi, linear between.
inline double upstreamGate(double fr, const DirectionalOperatorConfig& cfg) {
    if (!(fr > 0.0)) return 0.0;          // NaN / negative ⇒ no gate
    const double band = cfg.fr_hi - cfg.fr_lo;
    if (band <= 0.0) return fr >= cfg.fr_hi ? 1.0 : 0.0;
    return std::clamp((fr - cfg.fr_lo) / band, 0.0, 1.0);
}

/**
 * @brief Assemble the reduced directional operator M = Pᵀ·L_dir·P.
 *
 * @param P            Eigenbasis, k × n_active row-major (`P[j*n + i]`), the
 *                     same layout GraphEigenBasis stores.
 * @param k            Retained modes.
 * @param n_active     Active (non-outfall) node count.
 * @param n_conduits   Conduit count.
 * @param a1, a2       Per-conduit ACTIVE node indices of the two endpoints
 *                     (−1 for an outfall endpoint, which is grounded exactly
 *                     as NetworkLaplacian1D::buildWeighted grounds it).
 * @param weights      Per-conduit Laplacian weights, already normalized and
 *                     floored the way the basis build normalized them.
 * @param flow_sign    Per-conduit flow direction: +1 means a1 → a2 (the
 *                     conduit's nominal orientation), −1 the reverse, 0
 *                     unknown/still (no gate applied).
 * @param froude       Per-conduit Froude number (|Fr|; sign is ignored).
 * @param cfg          Gate dials.
 * @param drop_ratio   Optional per-conduit (h_up − h_dn)/d_up, oriented by
 *                     flow_sign (upstream = the sending end). Null disables
 *                     the drop gate; the applied gate is max(g_Fr, g_drop).
 * @param M_out        Resized to k×k row-major.
 * @param n_gated_out  Optional: number of conduits whose gate is > 0.
 * @return false on a dimension problem (k, n_active ≤ 0, null inputs).
 *
 * Assembly is per edge in O(E·k²): each conduit contributes
 *   M[p][q] += w·(P[p,u]·P[q,u] + P[p,d]·P[q,d] − (1−g)·P[p,u]·P[q,d] − P[p,d]·P[q,u])
 * for an interior conduit with upstream u and downstream d (flow_sign
 * oriented), and w·P[p,i]·P[q,i] for the interior endpoint i of a grounded
 * outfall conduit (the outfall head is fixed; no gating is needed there).
 */
inline bool assembleDirectionalReduced(const double* P, int k, int n_active,
                                       int n_conduits,
                                       const int* a1, const int* a2,
                                       const double* weights,
                                       const int* flow_sign,
                                       const double* froude,
                                       const DirectionalOperatorConfig& cfg,
                                       std::vector<double>& M_out,
                                       int* n_gated_out = nullptr,
                                       const double* drop_ratio = nullptr) {
    if (!P || k <= 0 || n_active <= 0 || n_conduits <= 0 ||
        !a1 || !a2 || !weights) return false;
    const auto uk = static_cast<std::size_t>(k);
    const auto un = static_cast<std::size_t>(n_active);
    M_out.assign(uk * uk, 0.0);
    int n_gated = 0;

    for (int c = 0; c < n_conduits; ++c) {
        const auto uc = static_cast<std::size_t>(c);
        const int i1 = a1[uc], i2 = a2[uc];
        const double w = weights[uc];
        if (w <= 0.0) continue;

        if (i1 >= 0 && i2 >= 0) {
            // Orient by flow: u = upstream (sends), d = downstream (receives).
            int u = i1, d = i2;
            double g = 0.0;
            if (flow_sign && (froude || drop_ratio)) {
                const int sgn = flow_sign[uc];
                if (sgn != 0) {
                    if (sgn < 0) std::swap(u, d);
                    if (froude)     g = upstreamGate(std::fabs(froude[uc]), cfg);
                    if (drop_ratio) g = std::max(g, dropGate(drop_ratio[uc], cfg));
                }
            }
            if (g > 0.0) ++n_gated;
            const double keep = 1.0 - g;   // surviving fraction of L[u,d]
            const auto uu = static_cast<std::size_t>(u), ud = static_cast<std::size_t>(d);
            for (std::size_t p = 0; p < uk; ++p) {
                const double Ppu = P[p * un + uu], Ppd = P[p * un + ud];
                double* Mp = &M_out[p * uk];
                for (std::size_t q = 0; q < uk; ++q) {
                    const double Pqu = P[q * un + uu], Pqd = P[q * un + ud];
                    Mp[q] += w * (Ppu * Pqu + Ppd * Pqd
                                  - keep * Ppu * Pqd - Ppd * Pqu);
                }
            }
        } else if (i1 >= 0 || i2 >= 0) {
            // Interior–outfall: grounded diagonal at the interior endpoint.
            const auto ui = static_cast<std::size_t>(i1 >= 0 ? i1 : i2);
            for (std::size_t p = 0; p < uk; ++p) {
                const double Ppi = P[p * un + ui];
                double* Mp = &M_out[p * uk];
                for (std::size_t q = 0; q < uk; ++q)
                    Mp[q] += w * Ppi * P[q * un + ui];
            }
        }
        // else: outfall–outfall, no active endpoint.
    }
    if (n_gated_out) *n_gated_out = n_gated;
    return true;
}

} // namespace openswmm::uncertainty

#endif // OPENSWMM_ENGINE_UNCERTAINTY_ROM_DIRECTIONAL_OPERATOR_HPP
