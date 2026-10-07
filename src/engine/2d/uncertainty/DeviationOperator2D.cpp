/**
 * @file DeviationOperator2D.cpp
 * @brief Reduced deviation operator — assembly and dense propagator.
 *
 * @see DeviationOperator2D.hpp
 * @ingroup engine_2d
 */

#include "DeviationOperator2D.hpp"

#include "uncertainty/GraphEigenBasis.hpp"
#include "uncertainty/RomDensePropagator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace openswmm::twoD {

using openswmm::uncertainty::CsrGraph;
using openswmm::uncertainty::coo_to_csr;
using openswmm::uncertainty::csr_matvec;

// ============================================================================
// assemble
// ============================================================================

bool DeviationOperator2D::assemble(const MeshData& mesh,
                                   const MeshEigenBasis& basis,
                                   double D_scale,
                                   const double* h_cell,
                                   const double* cell_u, const double* cell_v,
                                   const double* ground_w,
                                   const double* cond_mult) {
    k = 0;
    M.clear();
    if (!basis.is_ready()) return false;
    const int n = mesh.n_triangles();
    if (n < 4 || basis.n_triangles != n) return false;
    if (D_scale < 0.0) return false;

    const bool has_flow = (cell_u != nullptr) && (cell_v != nullptr);
    const bool advect   = has_flow && (c_factor != 0.0);

    // Optional depth weighting: per-cell (h/h̄)^{5/3}, floored so a dry cell
    // cannot zero an edge outright (same floor role as the eigenbasis's
    // depth-weighted build).
    std::vector<double> w_cell;
    if (h_cell != nullptr) {
        w_cell.assign(static_cast<std::size_t>(n), 1.0);
        double h_sum = 0.0;
        for (int t = 0; t < n; ++t)
            h_sum += std::max(h_cell[static_cast<std::size_t>(t)], 0.0);
        const double h_mean = h_sum / static_cast<double>(n);
        if (h_mean > 1.0e-9) {
            const double inv_h53 = 1.0 / std::pow(h_mean, 5.0 / 3.0);
            for (int t = 0; t < n; ++t) {
                const double h_t =
                    std::max(h_cell[static_cast<std::size_t>(t)], 0.0);
                w_cell[static_cast<std::size_t>(t)] =
                    std::max(std::pow(h_t, 5.0 / 3.0) * inv_h53, 1.0e-9);
            }
        }
    }

    // ---- Full-mesh operator L_op, physical FV convention -------------------
    // Diffusion row i: +Σ_j (D_edge·w_geo/A_i), off-diagonal −D_edge·w_geo/A_i.
    // Advection (first-order upwind): outflow faces take the local cell value,
    // inflow faces the neighbour's — dissipative by construction, and the part
    // of the physics a symmetric operator cannot represent.
    std::vector<std::vector<std::pair<int, double>>> coo(
        static_cast<std::size_t>(n));

    const double iso_mean = 0.5 * (alpha_par + alpha_perp);

    for (int i = 0; i < n; ++i) {
        const int tri_nbrs[3] = {
            mesh.tri_nbr0[static_cast<std::size_t>(i)],
            mesh.tri_nbr1[static_cast<std::size_t>(i)],
            mesh.tri_nbr2[static_cast<std::size_t>(i)]
        };
        for (int e = 0; e < 3; ++e) {
            const int j = tri_nbrs[e];
            if (j < 0 || j <= i) continue;  // each interior face once (i < j)

            const auto ui = static_cast<std::size_t>(i);
            const auto uj = static_cast<std::size_t>(j);

            const double len = mesh.edge_length[static_cast<std::size_t>(i * 3 + e)];
            const double dx  = mesh.tri_cx[uj] - mesh.tri_cx[ui];
            const double dy  = mesh.tri_cy[uj] - mesh.tri_cy[ui];
            const double d   = std::sqrt(dx * dx + dy * dy);
            if (d < 1e-14) continue;

            const double A_i = mesh.tri_area[ui];
            const double A_j = mesh.tri_area[uj];
            if (A_i < 1e-30 || A_j < 1e-30) continue;

            // Per-cell conductance multipliers (H15): 1 when not supplied.
            const double m_i = cond_mult ? cond_mult[ui] : 1.0;
            const double m_j = cond_mult ? cond_mult[uj] : 1.0;

            // Face flow: mean of the two cell velocities, each scaled by its
            // own cell's Manning factor (u ∝ 1/n).
            double fu = 0.0, fv = 0.0, speed = 0.0;
            if (has_flow) {
                fu = 0.5 * (cell_u[ui] * m_i + cell_u[uj] * m_j);
                fv = 0.5 * (cell_v[ui] * m_i + cell_v[uj] * m_j);
                speed = std::sqrt(fu * fu + fv * fv);
            }

            // Flow-aligned anisotropic conductance factor. cos θ projects the
            // centroid direction onto the flow direction; where the flow is too
            // slow to define a direction, the factor blends to the isotropic
            // mean (vector friction linearizes without a preferred axis at
            // rest — the anisotropy is a property of sustained flow).
            double aniso = iso_mean;
            if (has_flow && speed > u_eps) {
                const double ce = (dx * fu + dy * fv) / (d * speed);  // cos θ
                const double c2 = ce * ce;
                aniso = alpha_par * c2 + alpha_perp * (1.0 - c2);
            } else if (!has_flow) {
                // No velocity supplied at all: pure isotropic operator uses
                // the parallel dial as THE dial (α∥ = α⊥ expected).
                aniso = iso_mean;
            }

            double w_geo = len / d;  // dimensionless geometric conductance
            if (h_cell != nullptr) {
                const double wi = w_cell[ui], wj = w_cell[uj];
                w_geo *= 2.0 * wi * wj / (wi + wj);  // harmonic mean
            }
            if (cond_mult != nullptr) {
                // Series conductance across the face: harmonic mean of the
                // two cells' 1/n factors (reduces to c for uniform c).
                const double den = m_i + m_j;
                if (den <= 0.0) continue;
                w_geo *= 2.0 * m_i * m_j / den;
            }

            const double cond = D_scale * aniso * w_geo;  // m²/s · (len/d)

            // Diffusion, FV: divide each row's flux by that row's cell area.
            coo[ui].emplace_back(i,  cond / A_i);
            coo[ui].emplace_back(j, -cond / A_i);
            coo[uj].emplace_back(j,  cond / A_j);
            coo[uj].emplace_back(i, -cond / A_j);

            // Upwind advection at celerity c⃗ = c_factor·u⃗.
            if (advect && speed > u_eps) {
                const double nx = dx / d, ny = dy / d;   // outward normal of i toward j
                const double un = c_factor * (fu * nx + fv * ny);  // m/s, + = i→j
                const double q  = un * len;                        // m²/s (per unit depth)
                if (un > 0.0) {
                    // Outflow from i carries δh_i; inflow to j receives it.
                    coo[ui].emplace_back(i,  q / A_i);
                    coo[uj].emplace_back(i, -q / A_j);
                } else if (un < 0.0) {
                    // Flow j→i: outflow from j carries δh_j; i receives it.
                    coo[uj].emplace_back(j, -q / A_j);
                    coo[ui].emplace_back(j,  q / A_i);
                }
            }
        }
    }

    // Open-boundary grounding (see the header note): a diagonal-only edge to
    // a zero-deviation ghost, at the mean conductance dial — the boundary face
    // orientation is not carried per cell, and the constant is a calibration
    // dial like the rest.
    if (ground_w != nullptr) {
        for (int i = 0; i < n; ++i) {
            const auto ui = static_cast<std::size_t>(i);
            const double gw = ground_w[ui];
            if (gw <= 0.0) continue;
            const double A_i = mesh.tri_area[ui];
            if (A_i < 1e-30) continue;
            double w = gw;
            if (h_cell != nullptr) w *= w_cell[ui];
            if (cond_mult != nullptr) w *= cond_mult[ui];
            coo[ui].emplace_back(i, D_scale * iso_mean * w / A_i);
        }
    }

    CsrGraph L = coo_to_csr(coo, n);

    // ---- Galerkin projection: M[p][q] = P[:,p]ᵀ · L · P[:,q] ----------------
    k = basis.num_kept;
    const auto nk = static_cast<std::size_t>(k);
    const auto nn = static_cast<std::size_t>(n);
    M.assign(nk * nk, 0.0);

    std::vector<double> y(nn, 0.0);
    for (std::size_t q = 0; q < nk; ++q) {
        const double* Pq = &basis.P[q * nn];
        csr_matvec(L, Pq, y.data());
        for (std::size_t p = 0; p < nk; ++p) {
            const double* Pp = &basis.P[p * nn];
            double dot = 0.0;
            for (std::size_t t = 0; t < nn; ++t)
                dot += Pp[t] * y[t];
            M[p * nk + q] = dot;
        }
    }
    return true;
}

// ============================================================================
// expm / propagate — delegate to the shared dense propagator (PR H14 moved
// the implementation to uncertainty/RomDensePropagator.hpp so the 1D ROM can
// use it; the code is verbatim, so 2D results are bit-identical).
// ============================================================================

void DeviationOperator2D::expm(std::vector<double>& A, int n) {
    openswmm::uncertainty::denseExpm(A, n);
}

void DeviationOperator2D::propagate(const std::vector<double>& M, int k,
                                    double s, double dt,
                                    double* delta_a, const double* g) {
    openswmm::uncertainty::propagateDense(M, k, s, dt, delta_a, g);
}

} // namespace openswmm::twoD
