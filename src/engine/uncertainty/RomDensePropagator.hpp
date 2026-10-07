/**
 * @file RomDensePropagator.hpp
 * @brief Dense k×k matrix exponential and exact deviation step, shared by the
 *        2D reduced operator (DeviationOperator2D) and the 1D directional
 *        operator (PR H14).
 *
 * @details Header-only and dependency-free so the 1D ROM (compiled
 *          unconditionally) can use it without pulling in the 2D sources.
 *          The code is the former file-local implementation from
 *          DeviationOperator2D.cpp, moved verbatim; DeviationOperator2D::expm
 *          and ::propagate now delegate here, so 2D results are bit-identical
 *          to before the move.
 *
 * @ingroup engine_uncertainty
 */

#ifndef OPENSWMM_ENGINE_UNCERTAINTY_ROM_DENSE_PROPAGATOR_HPP
#define OPENSWMM_ENGINE_UNCERTAINTY_ROM_DENSE_PROPAGATOR_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace openswmm::uncertainty {

namespace dense_detail {

// C = A·B, all n×n row-major.
inline void matmul(const std::vector<double>& A, const std::vector<double>& B,
                   std::vector<double>& C, int n) {
    const auto un = static_cast<std::size_t>(n);
    C.assign(un * un, 0.0);
    for (std::size_t i = 0; i < un; ++i) {
        const double* Ai = &A[i * un];
        double* Ci = &C[i * un];
        for (std::size_t p = 0; p < un; ++p) {
            const double a = Ai[p];
            if (a == 0.0) continue;
            const double* Bp = &B[p * un];
            for (std::size_t j = 0; j < un; ++j)
                Ci[j] += a * Bp[j];
        }
    }
}

// Solve D·X = N in place (X returned in N). Dense LU, partial pivoting.
// n is small (ROM modes + 1); this is not a performance path.
inline bool solveInPlace(std::vector<double>& D, std::vector<double>& N, int n) {
    const auto un = static_cast<std::size_t>(n);
    for (int col = 0; col < n; ++col) {
        int piv = col;
        double best = std::fabs(D[static_cast<std::size_t>(col) * un + col]);
        for (int r = col + 1; r < n; ++r) {
            const double v = std::fabs(D[static_cast<std::size_t>(r) * un + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-300) return false;  // singular — cannot happen for Padé D
        if (piv != col) {
            for (int c = 0; c < n; ++c) {
                std::swap(D[static_cast<std::size_t>(piv) * un + c],
                          D[static_cast<std::size_t>(col) * un + c]);
                std::swap(N[static_cast<std::size_t>(piv) * un + c],
                          N[static_cast<std::size_t>(col) * un + c]);
            }
        }
        const double inv = 1.0 / D[static_cast<std::size_t>(col) * un + col];
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            const double f = D[static_cast<std::size_t>(r) * un + col] * inv;
            if (f == 0.0) continue;
            for (int c = 0; c < n; ++c) {
                D[static_cast<std::size_t>(r) * un + c] -=
                    f * D[static_cast<std::size_t>(col) * un + c];
                N[static_cast<std::size_t>(r) * un + c] -=
                    f * N[static_cast<std::size_t>(col) * un + c];
            }
        }
        for (int c = 0; c < n; ++c)
            N[static_cast<std::size_t>(col) * un + c] *= inv;
        for (int c = 0; c < n; ++c)
            D[static_cast<std::size_t>(col) * un + c] *= inv;
    }
    return true;
}

} // namespace dense_detail

/**
 * @brief In-place dense matrix exponential A ← exp(A), n×n row-major.
 *
 * Scaling-and-squaring with a [6/6] Padé approximant — the standard
 * dependency-free algorithm, accurate to near round-off for well-scaled
 * matrices (‖A‖∞ is reduced below 1/2 before the Padé step).
 */
inline void denseExpm(std::vector<double>& A, int n) {
    const auto un = static_cast<std::size_t>(n);
    assert(A.size() == un * un);

    double norm = 0.0;
    for (std::size_t i = 0; i < un; ++i) {
        double rs = 0.0;
        for (std::size_t j = 0; j < un; ++j) rs += std::fabs(A[i * un + j]);
        norm = std::max(norm, rs);
    }
    int s = 0;
    if (norm > 0.5) {
        s = static_cast<int>(std::ceil(std::log2(norm / 0.5)));
        const double f = std::ldexp(1.0, -s);  // 2^{-s}
        for (auto& a : A) a *= f;
    }

    static const double c[7] = {
        1.0, 1.0 / 2.0, 5.0 / 44.0, 1.0 / 66.0,
        1.0 / 792.0, 1.0 / 15840.0, 1.0 / 665280.0
    };

    std::vector<double> Apow(un * un, 0.0);
    for (std::size_t i = 0; i < un; ++i) Apow[i * un + i] = 1.0;
    std::vector<double> Nmat(un * un, 0.0), Dmat(un * un, 0.0), tmp;
    for (std::size_t i = 0; i < un; ++i) {
        Nmat[i * un + i] = c[0];
        Dmat[i * un + i] = c[0];
    }
    double sign = 1.0;
    for (int m = 1; m <= 6; ++m) {
        dense_detail::matmul(Apow, A, tmp, n);
        Apow.swap(tmp);
        sign = -sign;
        for (std::size_t idx = 0; idx < un * un; ++idx) {
            Nmat[idx] += c[m] * Apow[idx];
            Dmat[idx] += sign * c[m] * Apow[idx];
        }
    }

    const bool ok = dense_detail::solveInPlace(Dmat, Nmat, n);
    assert(ok && "Padé denominator singular — matrix not properly scaled");
    (void)ok;

    for (int r = 0; r < s; ++r) {
        dense_detail::matmul(Nmat, Nmat, tmp, n);
        Nmat.swap(tmp);
    }
    A.swap(Nmat);
}

/**
 * @brief One exact deviation step:
 *        δa ← exp(−s·Δt·M)·δa + φ₁(−s·Δt·M)·Δt·g.
 *
 * φ₁ is evaluated via the augmented-matrix identity
 * exp([[A, Δt·g],[0,0]]) = [[e^A, φ₁(A)·Δt·g],[0,1]], well-defined for a
 * singular A (φ₁(0) = I ⇒ the Euler limit δa += g·Δt). For k = 1 this is the
 * scalar exact exponential integrator identically.
 *
 * @param M        k×k row-major operator (1/s).
 * @param k        Dimension.
 * @param s        Per-member operator scaling (1/mm_i for a Manning multiplier).
 * @param dt       Step (s).
 * @param delta_a  In/out: member deviation coefficients (length k).
 * @param g        Constant-over-step modal forcing (length k).
 */
inline void propagateDense(const std::vector<double>& M, int k,
                           double s, double dt,
                           double* delta_a, const double* g) {
    assert(k > 0);
    const auto uk = static_cast<std::size_t>(k);
    assert(M.size() == uk * uk);

    const int n = k + 1;
    const auto un = static_cast<std::size_t>(n);
    std::vector<double> B(un * un, 0.0);
    const double f = -s * dt;
    for (std::size_t i = 0; i < uk; ++i) {
        for (std::size_t j = 0; j < uk; ++j)
            B[i * un + j] = f * M[i * uk + j];
        B[i * un + uk] = dt * g[i];
    }

    denseExpm(B, n);

    std::vector<double> out(uk, 0.0);
    for (std::size_t i = 0; i < uk; ++i) {
        double v = B[i * un + uk];
        for (std::size_t j = 0; j < uk; ++j)
            v += B[i * un + j] * delta_a[j];
        out[i] = v;
    }
    for (std::size_t i = 0; i < uk; ++i) delta_a[i] = out[i];
}

} // namespace openswmm::uncertainty

#endif // OPENSWMM_ENGINE_UNCERTAINTY_ROM_DENSE_PROPAGATOR_HPP
