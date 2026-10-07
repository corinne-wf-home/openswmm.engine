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


// ============================================================================
// PR H14b — batch propagator: many scalars s_i, one operator A
// ============================================================================

/**
 * @brief Exact per-member deviation steps that share the operator's matrix
 *        powers across members.
 *
 * Every member of a ROM ensemble integrates the SAME operator scaled by its
 * own scalar (s_i = 1/mm_i): δa_i ← exp(s_i·A)·δa_i + φ₁(s_i·A)·c_i with
 * A = −Δt·K·M and c_i = Δt·g_i. propagateDense() pays a full Padé
 * exponential per member. This class evaluates Higham's degree-13 Padé
 * approximant (the standard expm algorithm, θ₁₃ = 5.37) on the augmented
 * matrix B_i = [[s_i·A, c_i],[0, 0]], whose powers are
 * B_i^m = [[s_i^m·A^m, s_i^{m−1}·A^{m−1}·c_i],[0, 0]]: the (1,1) blocks are
 * s-weighted combinations of A^m, computed ONCE per prepare(), and the last
 * column needs only matrix–vector products. Per member that leaves the
 * O(k²) combination plus one (k+1)×(k+1) linear solve, versus ~7 dense
 * products per member before. The scaling exponent is chosen once from the
 * LARGEST |s_i| (scaling more than a member needs is harmless); squarings,
 * when any are needed, are per member.
 *
 * Numerically this is the same approximant family as denseExpm (a higher
 * Padé order, no more rounding); results agree with propagateDense to
 * ~1e-13 relative (tested), not bit-for-bit.
 *
 * Usage per step:
 *     batch.prepare(A, k, s_max);            // once
 *     for each member: batch.apply(s_i, c_i, delta_a_i);
 */
class BatchPropagator {
public:
    /// Highest Padé order kept: powers A^1..A^13 are formed.
    static constexpr int kOrder = 13;

    /**
     * @param A      k×k row-major operator (already scaled by Δt and any
     *               rate constant, e.g. −Δt·K1d·M).
     * @param k      Dimension (≥ 1).
     * @param s_max  max_i |s_i| over the members that will call apply().
     */
    void prepare(const std::vector<double>& A, int k, double s_max) {
        k_ = k;
        const auto uk = static_cast<std::size_t>(k);
        assert(A.size() == uk * uk);
        // ‖B‖₁ bound for the scaled augmented matrix: |s|·‖A‖₁ (+ the
        // forcing column, which only ever lowers the needed scaling by a
        // bounded amount; it is re-checked per member in apply()).
        double norm1 = 0.0;
        for (std::size_t j = 0; j < uk; ++j) {
            double cs = 0.0;
            for (std::size_t i = 0; i < uk; ++i) cs += std::fabs(A[i * uk + j]);
            norm1 = std::max(norm1, cs);
        }
        scale_j_ = 0;
        const double bound = s_max * norm1;
        if (bound > kTheta13)
            scale_j_ = static_cast<int>(std::ceil(std::log2(bound / kTheta13)));
        scale_f_ = std::ldexp(1.0, -scale_j_);

        // Powers of the (scaled) operator: pow_[m] = (f·A)^m, m = 1..13.
        pow_.assign(static_cast<std::size_t>(kOrder) + 1, std::vector<double>());
        pow_[1].resize(uk * uk);
        for (std::size_t idx = 0; idx < uk * uk; ++idx) pow_[1][idx] = scale_f_ * A[idx];
        for (int m = 2; m <= kOrder; ++m)
            dense_detail::matmul(pow_[static_cast<std::size_t>(m - 1)], pow_[1],
                                 pow_[static_cast<std::size_t>(m)], k);
        prepared_ = true;
    }

    bool prepared() const noexcept { return prepared_; }
    int  scalingExponent() const noexcept { return scale_j_; }

    /**
     * @brief One member: δa ← exp(s·A)·δa + φ₁(s·A)·c.
     * @param s        Member scalar (|s| ≤ s_max given to prepare()).
     * @param c        Forcing integral column (length k): Δt·g.
     * @param delta_a  In/out (length k).
     */
    void apply(double s, const double* c, double* delta_a) {
        assert(prepared_);
        const int k = k_;
        const auto uk = static_cast<std::size_t>(k);
        const int n = k + 1;
        const auto un = static_cast<std::size_t>(n);

        // Scaled forcing column and its images under the scaled powers:
        // v[m] = (f·A)^m · (f·c), m = 0..12 (we need A^{m-1}c for m ≤ 13).
        // Scratch buffers live on the object: apply() runs M times per step
        // for 10⁵ steps, and per-call allocation measurably dominated the
        // arithmetic in an optimized build.
        std::vector<double>& v = scr_v_;
        v.resize(static_cast<std::size_t>(kOrder) * uk);
        for (std::size_t i = 0; i < uk; ++i) v[i] = scale_f_ * c[i];
        for (int m = 1; m < kOrder; ++m) {
            const double* src = &v[static_cast<std::size_t>(m - 1) * uk];
            double* dst = &v[static_cast<std::size_t>(m) * uk];
            matvec(pow_[1], src, dst, k);
        }

        // Padé[13/13]: U = Σ_{m odd} b_m B^m, V = Σ_{m even} b_m B^m, with
        // B^m = [[s^m A^m, s^{m-1} A^{m-1} c],[0,0]] and B^0 = I.
        std::vector<double>& U = scr_U_;
        std::vector<double>& V = scr_V_;
        U.assign(un * un, 0.0);
        V.assign(un * un, 0.0);
        double sp = 1.0;   // s^m
        for (int m = 0; m <= kOrder; ++m) {
            const double bm = kB13[m];
            std::vector<double>& T = (m % 2 == 1) ? U : V;
            if (m == 0) {
                for (std::size_t i = 0; i < un; ++i) T[i * un + i] += bm;
            } else {
                const std::vector<double>& Am = pow_[static_cast<std::size_t>(m)];
                const double w = bm * sp;
                for (std::size_t i = 0; i < uk; ++i)
                    for (std::size_t j = 0; j < uk; ++j)
                        T[i * un + j] += w * Am[i * uk + j];
                const double wc = bm * sp / s;   // s^{m-1}
                const double* vm = &v[static_cast<std::size_t>(m - 1) * uk];
                for (std::size_t i = 0; i < uk; ++i)
                    T[i * un + uk] += wc * vm[i];
            }
            sp *= s;
        }
        // X = (V − U)^{-1} (V + U)
        std::vector<double>& D = scr_D_;
        std::vector<double>& N = scr_N_;
        D.resize(un * un);
        N.resize(un * un);
        for (std::size_t idx = 0; idx < un * un; ++idx) {
            D[idx] = V[idx] - U[idx];
            N[idx] = V[idx] + U[idx];
        }
        const bool ok = dense_detail::solveInPlace(D, N, n);
        assert(ok && "Padé denominator singular");
        (void)ok;
        // Undo the scaling: square scale_j_ times (per member — the scalar
        // is inside the exponent, so this cannot be shared).
        std::vector<double>& tmp = scr_tmp_;
        for (int r = 0; r < scale_j_; ++r) {
            dense_detail::matmul(N, N, tmp, n);
            N.swap(tmp);
        }
        std::vector<double>& out = scr_out_;
        out.assign(uk, 0.0);
        for (std::size_t i = 0; i < uk; ++i) {
            double y = N[i * un + uk];
            for (std::size_t j = 0; j < uk; ++j) y += N[i * un + j] * delta_a[j];
            out[i] = y;
        }
        for (std::size_t i = 0; i < uk; ++i) delta_a[i] = out[i];
    }

private:
    static constexpr double kTheta13 = 5.371920351148152;
    // Higham (2005) Padé[13/13] coefficients b_0..b_13.
    static constexpr double kB13[14] = {
        64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
        1187353796428800.0, 129060195264000.0, 10559470521600.0,
        670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
        960960.0, 16380.0, 182.0, 1.0 };

    static void matvec(const std::vector<double>& A, const double* x, double* y, int k) {
        const auto uk = static_cast<std::size_t>(k);
        for (std::size_t i = 0; i < uk; ++i) {
            double acc = 0.0;
            const double* Ai = &A[i * uk];
            for (std::size_t j = 0; j < uk; ++j) acc += Ai[j] * x[j];
            y[i] = acc;
        }
    }

    int k_ = 0;
    int scale_j_ = 0;
    double scale_f_ = 1.0;
    bool prepared_ = false;
    std::vector<std::vector<double>> pow_;
    // apply() scratch (reused across members and steps).
    std::vector<double> scr_v_, scr_U_, scr_V_, scr_D_, scr_N_, scr_tmp_, scr_out_;
};

} // namespace openswmm::uncertainty

#endif // OPENSWMM_ENGINE_UNCERTAINTY_ROM_DENSE_PROPAGATOR_HPP
