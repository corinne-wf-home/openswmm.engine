/**
 * @file test_rom_directional_operator.cpp
 * @brief PR H14 — Froude-gated directional 1D operator: pure-function
 *        contracts (gate ramp, assembly vs a hand-built dense L_dir, the
 *        "nothing travels upstream" property) and the SpectralROM1D
 *        reduced-operator path (diag(λ) reproduces the diagonal integrator;
 *        clear/throw contracts).
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "uncertainty/GraphEigenBasis.hpp"
#include "uncertainty/NetworkLaplacian1D.hpp"
#include "uncertainty/RomDensePropagator.hpp"
#include "uncertainty/RomDirectionalOperator.hpp"
#include "uncertainty/SpectralROM1D.hpp"

using namespace openswmm::uncertainty;

namespace {

// Chain J0 → J1 → … → J(n−2) → outfall J(n−1); conduit c joins c and c+1.
struct Chain {
    int n_full, n_conduits, n_active = 0;
    std::vector<int> n1, n2, is_outfall, active_map, full_to_active, a1, a2;
    std::vector<double> w;
    CsrGraph L;
    GraphEigenBasis basis;

    explicit Chain(int n, int k_req, double w_base = 1.0, double w_slope = 0.0) : n_full(n), n_conduits(n - 1) {
        for (int c = 0; c < n_conduits; ++c) {
            n1.push_back(c); n2.push_back(c + 1);
            w.push_back(w_base + w_slope * c);
        }
        is_outfall.assign(static_cast<std::size_t>(n), 0);
        is_outfall.back() = 1;
        L = NetworkLaplacian1D::buildWeighted(n_full, n_conduits, n1.data(), n2.data(),
                                              is_outfall.data(), w.data(),
                                              active_map, full_to_active);
        n_active = static_cast<int>(active_map.size());
        for (int c = 0; c < n_conduits; ++c) {
            a1.push_back(full_to_active[static_cast<std::size_t>(n1[static_cast<std::size_t>(c)])]);
            a2.push_back(full_to_active[static_cast<std::size_t>(n2[static_cast<std::size_t>(c)])]);
        }
        if (!basis.build(L, k_req)) throw std::runtime_error("basis build failed");
    }
};

// Hand-built dense L_dir in ACTIVE space (row-major n×n), the reference the
// reduced assembly is checked against.
std::vector<double> denseLdir(const Chain& ch, const std::vector<int>& sign,
                              const std::vector<double>& fr, const DirectionalOperatorConfig& cfg) {
    const auto n = static_cast<std::size_t>(ch.n_active);
    std::vector<double> L(n * n, 0.0);
    for (int c = 0; c < ch.n_conduits; ++c) {
        const auto uc = static_cast<std::size_t>(c);
        const int i1 = ch.a1[uc], i2 = ch.a2[uc];
        const double w = ch.w[uc];
        if (i1 >= 0 && i2 >= 0) {
            int u = i1, d = i2;
            double g = 0.0;
            if (sign[uc] != 0) { if (sign[uc] < 0) std::swap(u, d); g = upstreamGate(std::fabs(fr[uc]), cfg); }
            const auto uu = static_cast<std::size_t>(u), ud = static_cast<std::size_t>(d);
            L[uu * n + uu] += w; L[ud * n + ud] += w;
            L[uu * n + ud] -= w * (1.0 - g);   // upstream row: gated coupling
            L[ud * n + uu] -= w;               // downstream row: full coupling
        } else if (i1 >= 0 || i2 >= 0) {
            const auto ui = static_cast<std::size_t>(i1 >= 0 ? i1 : i2);
            L[ui * n + ui] += w;
        }
    }
    return L;
}

double maxAbsDiff(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

} // namespace

// ─── Gate ramp ──────────────────────────────────────────────────────────────

TEST(DirectionalGate, RampClosedForms) {
    DirectionalOperatorConfig cfg;   // [0.8, 1.2]
    EXPECT_EQ(upstreamGate(0.0, cfg), 0.0);
    EXPECT_EQ(upstreamGate(0.5, cfg), 0.0);
    EXPECT_EQ(upstreamGate(0.8, cfg), 0.0);
    EXPECT_NEAR(upstreamGate(1.0, cfg), 0.5, 1e-15);
    EXPECT_EQ(upstreamGate(1.2, cfg), 1.0);
    EXPECT_EQ(upstreamGate(3.0, cfg), 1.0);
    EXPECT_EQ(upstreamGate(std::numeric_limits<double>::quiet_NaN(), cfg), 0.0);
    EXPECT_EQ(upstreamGate(-1.5, cfg), 0.0);
    DirectionalOperatorConfig step; step.fr_lo = step.fr_hi = 1.0;   // degenerate band
    EXPECT_EQ(upstreamGate(0.99, step), 0.0);
    EXPECT_EQ(upstreamGate(1.0, step), 1.0);
}

TEST(DirectionalGate, DropGateClosedForms) {
    DirectionalOperatorConfig cfg;   // [1.0, 2.0]
    cfg.use_drop = true;             // ships off; tested on
    EXPECT_EQ(dropGate(0.0, cfg), 0.0);
    EXPECT_EQ(dropGate(0.9, cfg), 0.0);
    EXPECT_EQ(dropGate(1.0, cfg), 0.0);
    EXPECT_NEAR(dropGate(1.5, cfg), 0.5, 1e-15);
    EXPECT_EQ(dropGate(2.0, cfg), 1.0);
    EXPECT_EQ(dropGate(20.0, cfg), 1.0);
    EXPECT_EQ(dropGate(std::numeric_limits<double>::quiet_NaN(), cfg), 0.0);
    DirectionalOperatorConfig off = cfg; off.use_drop = false;
    EXPECT_EQ(dropGate(20.0, off), 0.0);
    EXPECT_FALSE(DirectionalOperatorConfig{}.use_drop) << "drop gate must ship off (measured regression on CL-1e)";
}

TEST(DirectionalAssembly, DropGateOpensADrownedSteepReachTheFroudeGateMisses) {
    // A reach that reads subcritical at mid-depth (Fr 0.3) but whose upstream
    // node sits 10 depths above the downstream water surface: the Froude gate
    // leaves it coupled, the drop gate makes it one-way, and the combined
    // gate equals the drop-only result.
    Chain ch(7, 3);
    const int k = ch.basis.num_kept;
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr(static_cast<std::size_t>(ch.n_conduits), 0.3);
    std::vector<double> drop(static_cast<std::size_t>(ch.n_conduits), 0.2);
    drop[2] = 10.0;   // J2 -> J3 drowned-outlet steep reach
    DirectionalOperatorConfig cfg;
    cfg.use_drop = true;
    std::vector<double> M_fr, M_both, M_drop_only;
    int ng_fr = -1, ng_both = -1;
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr.data(), cfg, M_fr, &ng_fr));
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr.data(), cfg, M_both, &ng_both, drop.data()));
    EXPECT_EQ(ng_fr, 0);
    EXPECT_EQ(ng_both, 1);
    // Drop-only reference: Froude all zero.
    std::vector<double> fr0(static_cast<std::size_t>(ch.n_conduits), 0.0);
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr0.data(), cfg, M_drop_only, nullptr, drop.data()));
    EXPECT_EQ(maxAbsDiff(M_both, M_drop_only), 0.0);
    EXPECT_GT(maxAbsDiff(M_both, M_fr), 1e-3);
}

// ─── Assembly contracts ─────────────────────────────────────────────────────

TEST(DirectionalAssembly, UngatedOperatorIsDiagLambdaUpToRitzResidual) {
    Chain ch(8, 4, 1.0, 0.15);   // 7 active, non-uniform weights
    const int k = ch.basis.num_kept;
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr(static_cast<std::size_t>(ch.n_conduits), 0.3);   // all subcritical
    DirectionalOperatorConfig cfg;
    std::vector<double> M; int ng = -1;
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(),
                                           sign.data(), fr.data(), cfg, M, &ng));
    EXPECT_EQ(ng, 0);
    double lam_max = 0.0;
    for (int j = 0; j < k; ++j) lam_max = std::max(lam_max, ch.basis.eigenvalues[static_cast<std::size_t>(j)]);
    for (int p = 0; p < k; ++p)
        for (int q = 0; q < k; ++q) {
            const double expect = (p == q) ? ch.basis.eigenvalues[static_cast<std::size_t>(p)] : 0.0;
            EXPECT_NEAR(M[static_cast<std::size_t>(p * k + q)], expect, 1e-6 * lam_max) << p << "," << q;
        }
}

TEST(DirectionalAssembly, NullDirectionInputsAreBitIdenticalToAllSubcritical) {
    Chain ch(7, 3);
    const int k = ch.basis.num_kept;
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr(static_cast<std::size_t>(ch.n_conduits), 0.0);
    DirectionalOperatorConfig cfg;
    std::vector<double> Ma, Mb;
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), nullptr, nullptr, cfg, Ma));
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr.data(), cfg, Mb));
    ASSERT_EQ(Ma.size(), Mb.size());
    for (std::size_t i = 0; i < Ma.size(); ++i) EXPECT_EQ(Ma[i], Mb[i]);
}

TEST(DirectionalAssembly, MatchesHandBuiltDenseProjectionWithMixedGates) {
    Chain ch(9, 4, 1.0, 0.2);
    const int k = ch.basis.num_kept;
    const auto n = static_cast<std::size_t>(ch.n_active);
    // Mixed: some reaches supercritical, one reversed, one on the ramp.
    std::vector<int> sign = {1, 1, -1, 1, 0, 1, 1, 1};
    std::vector<double> fr = {1.6, 0.4, 1.9, 1.0, 2.0, 0.0, 1.3, 0.9};
    DirectionalOperatorConfig cfg;
    std::vector<double> M; int ng = -1;
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(),
                                           sign.data(), fr.data(), cfg, M, &ng));
    // Recount by the gate itself so the assertion documents the rule, not a guess.
    int expect_gated = 0;
    for (std::size_t c = 0; c < fr.size(); ++c)
        if (ch.a1[c] >= 0 && ch.a2[c] >= 0 &&     // outfall conduits are grounded, never gated
            sign[c] != 0 && upstreamGate(fr[c], cfg) > 0.0) ++expect_gated;
    EXPECT_EQ(ng, expect_gated);

    const std::vector<double> L = denseLdir(ch, sign, fr, cfg);
    // Reference: M_ref[p][q] = Σ_i Σ_j P[p,i] L[i,j] P[q,j]
    std::vector<double> Mref(static_cast<std::size_t>(k) * static_cast<std::size_t>(k), 0.0);
    for (int p = 0; p < k; ++p)
        for (int q = 0; q < k; ++q) {
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    s += ch.basis.P[static_cast<std::size_t>(p) * n + i] * L[i * n + j] *
                         ch.basis.P[static_cast<std::size_t>(q) * n + j];
            Mref[static_cast<std::size_t>(p * k + q)] = s;
        }
    EXPECT_LT(maxAbsDiff(M, Mref), 1e-12);
    // And it is genuinely non-symmetric once a gate is open.
    double asym = 0.0;
    for (int p = 0; p < k; ++p)
        for (int q = 0; q < k; ++q)
            asym = std::max(asym, std::fabs(M[static_cast<std::size_t>(p * k + q)] - M[static_cast<std::size_t>(q * k + p)]));
    EXPECT_GT(asym, 1e-3);
}

TEST(DirectionalOperator, FullyGatedChainCarriesNothingUpstream) {
    // Dense L_dir on a 6-active chain, every reach supercritical: a deviation
    // injected at the most downstream active node must stay exactly there
    // (its row has no upstream neighbour, and no upstream row sees it),
    // whereas the symmetric Laplacian spreads it upstream within one step.
    Chain ch(7, 3);
    const auto n = static_cast<std::size_t>(ch.n_active);
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr_super(static_cast<std::size_t>(ch.n_conduits), 2.0);
    std::vector<double> fr_sub(static_cast<std::size_t>(ch.n_conduits), 0.0);
    DirectionalOperatorConfig cfg;
    std::vector<double> Ldir = denseLdir(ch, sign, fr_super, cfg);
    std::vector<double> Lsym = denseLdir(ch, sign, fr_sub, cfg);

    std::vector<double> g(n, 0.0);
    std::vector<double> da_dir(n, 0.0), da_sym(n, 0.0);
    da_dir[n - 1] = da_sym[n - 1] = 1.0;    // deviation at the downstream-most active node
    for (int step = 0; step < 5; ++step) {
        propagateDense(Ldir, static_cast<int>(n), 1.0, 0.7, da_dir.data(), g.data());
        propagateDense(Lsym, static_cast<int>(n), 1.0, 0.7, da_sym.data(), g.data());
    }
    for (std::size_t i = 0; i + 1 < n; ++i)
        EXPECT_NEAR(da_dir[i], 0.0, 1e-14) << "upstream node " << i << " received a downstream deviation";
    EXPECT_GT(std::fabs(da_sym[n - 2]), 1e-3) << "symmetric reference must leak upstream (else the test is vacuous)";
    // The downstream node itself still decays (its diagonal is intact).
    EXPECT_LT(da_dir[n - 1], 1.0);
    EXPECT_GT(da_dir[n - 1], 0.0);
}

TEST(DirectionalOperator, UpstreamDeviationStillReachesDownstream) {
    // The gate is one-way: an upstream deviation must still propagate down.
    Chain ch(7, 3);
    const auto n = static_cast<std::size_t>(ch.n_active);
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr(static_cast<std::size_t>(ch.n_conduits), 2.0);
    DirectionalOperatorConfig cfg;
    std::vector<double> Ldir = denseLdir(ch, sign, fr, cfg);
    std::vector<double> g(n, 0.0), da(n, 0.0);
    da[0] = 1.0;
    propagateDense(Ldir, static_cast<int>(n), 1.0, 0.7, da.data(), g.data());
    EXPECT_GT(da[1], 1e-3);
}

// ─── Dense propagator scalar contract ───────────────────────────────────────

TEST(DensePropagator, ScalarCaseIsTheExactExponentialIntegrator) {
    std::vector<double> M = {0.9};
    const double s = 1.25, dt = 7.0, g = 0.3;
    double a = 0.2;
    propagateDense(M, 1, s, dt, &a, &g);
    const double rate = s * 0.9, steady = g / rate;
    EXPECT_NEAR(a, (0.2 - steady) * std::exp(-rate * dt) + steady, 1e-13);
}

namespace {

// Deterministic pseudo-random dense k×k (non-symmetric) and vector.
void fillRandom(std::vector<double>& A, std::vector<double>& c, int k, unsigned seed, double scale) {
    unsigned x = seed;
    auto rnd = [&]() { x = x * 1664525u + 1013904223u; return (static_cast<double>(x >> 8) / 16777216.0) * 2.0 - 1.0; };
    A.assign(static_cast<std::size_t>(k) * static_cast<std::size_t>(k), 0.0);
    c.assign(static_cast<std::size_t>(k), 0.0);
    for (auto& v : A) v = scale * rnd();
    for (auto& v : c) v = rnd();
}

} // namespace

TEST(BatchPropagator, MatchesPerMemberPadeOnSmallNormOperator) {
    const int k = 9;
    std::vector<double> A, c;
    fillRandom(A, c, k, 7u, 0.08);          // ‖A‖ well below θ₁₃: no scaling
    const double s_vals[] = {0.83, 0.91, 1.0, 1.07, 1.19};
    BatchPropagator batch;
    batch.prepare(A, k, 1.19);
    EXPECT_EQ(batch.scalingExponent(), 0);
    for (double s : s_vals) {
        std::vector<double> da_ref(static_cast<std::size_t>(k)), da_bat(static_cast<std::size_t>(k));
        for (int i = 0; i < k; ++i) da_ref[static_cast<std::size_t>(i)] = da_bat[static_cast<std::size_t>(i)] = 0.3 * std::sin(1.3 * i + s);
        // Reference: propagateDense integrates exp(−s·dt·M); with M = −A/dt and dt = 1 that is exp(s·A).
        std::vector<double> M(A.size());
        for (std::size_t i = 0; i < A.size(); ++i) M[i] = -A[i];
        propagateDense(M, k, s, 1.0, da_ref.data(), c.data());
        batch.apply(s, c.data(), da_bat.data());
        double mx = 0.0;
        for (double v : da_ref) mx = std::max(mx, std::fabs(v));
        EXPECT_LT(maxAbsDiff(da_ref, da_bat), 1e-13 * std::max(mx, 1.0)) << "s=" << s;
    }
}

TEST(BatchPropagator, MatchesPerMemberPadeWhenScalingIsNeeded) {
    const int k = 12;
    std::vector<double> A, c;
    fillRandom(A, c, k, 11u, 1.5);          // ‖A‖₁ ≫ θ₁₃: several squarings
    BatchPropagator batch;
    batch.prepare(A, k, 1.25);
    EXPECT_GT(batch.scalingExponent(), 0);
    for (double s : {0.8, 1.0, 1.25}) {
        std::vector<double> da_ref(static_cast<std::size_t>(k), 0.1), da_bat(static_cast<std::size_t>(k), 0.1);
        std::vector<double> M(A.size());
        for (std::size_t i = 0; i < A.size(); ++i) M[i] = -A[i];
        propagateDense(M, k, s, 1.0, da_ref.data(), c.data());
        batch.apply(s, c.data(), da_bat.data());
        double mx = 0.0;
        for (double v : da_ref) mx = std::max(mx, std::fabs(v));
        ASSERT_GT(mx, 0.0);
        EXPECT_LT(maxAbsDiff(da_ref, da_bat), 1e-11 * mx) << "s=" << s;
    }
}

TEST(BatchPropagator, ScalarAndSingularLimits) {
    // k = 1: exact exponential integrator; A = 0: Euler limit (φ₁(0) = I).
    std::vector<double> A = {-0.9 * 7.0};   // −dt·rate, dt = 7
    BatchPropagator b1;
    b1.prepare(A, 1, 1.25);
    double da = 0.2;
    const double c = 0.3 * 7.0;   // dt·g
    b1.apply(1.25, &c, &da);
    const double rate = 1.25 * 0.9, steady = 0.3 / rate;
    EXPECT_NEAR(da, (0.2 - steady) * std::exp(-rate * 7.0) + steady, 1e-13);

    std::vector<double> Z(9, 0.0), cz = {1.0, -2.0, 0.5};
    BatchPropagator b0;
    b0.prepare(Z, 3, 1.1);
    std::vector<double> d = {0.1, 0.2, 0.3};
    b0.apply(1.1, cz.data(), d.data());
    EXPECT_NEAR(d[0], 0.1 + 1.0, 1e-15);
    EXPECT_NEAR(d[1], 0.2 - 2.0, 1e-15);
    EXPECT_NEAR(d[2], 0.3 + 0.5, 1e-15);
}

// ─── SpectralROM1D reduced path ─────────────────────────────────────────────

namespace {

void makeRom(SpectralROM1D& rom, const GraphEigenBasis& basis, int M, const std::vector<double>& h0) {
    rom.basis = &basis;
    rom.n_ensemble = M;
    rom.mannings_pert = 0.20;
    rom.runoff_pert = 0.15;
    rom.initialize();
    rom.seed(h0.data());
}

} // namespace

TEST(SpectralROM1DReduced, DiagLambdaOperatorReproducesDiagonalPath) {
    Chain ch(10, 5, 1.0, 0.1);
    const int k = ch.basis.num_kept;
    const auto n = static_cast<std::size_t>(ch.n_active);
    std::vector<double> h(n), sens(n), runoff(n);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = 100.0 - 2.0 * static_cast<double>(i) + 0.3 * std::sin(0.9 * static_cast<double>(i));
        sens[i] = 0.4 + 0.2 * std::cos(0.7 * static_cast<double>(i));
        runoff[i] = 1e-5 * (1.0 + 0.5 * std::sin(0.4 * static_cast<double>(i)));
    }
    SpectralROM1D diag, red;
    makeRom(diag, ch.basis, 12, h);
    makeRom(red,  ch.basis, 12, h);
    std::vector<double> Mdiag(static_cast<std::size_t>(k) * static_cast<std::size_t>(k), 0.0);
    for (int j = 0; j < k; ++j) Mdiag[static_cast<std::size_t>(j * k + j)] = ch.basis.eigenvalues[static_cast<std::size_t>(j)];
    red.setReducedOperator(Mdiag);
    ASSERT_TRUE(red.hasReducedOperator());

    for (int step = 0; step < 6; ++step) {
        diag.advance(30.0, 0.02, h.data(), runoff.data(), sens.data(), nullptr);
        red.advance(30.0, 0.02, h.data(), runoff.data(), sens.data(), nullptr);
    }
    double max_a = 0.0;
    for (double a : diag.a_ensemble) max_a = std::max(max_a, std::fabs(a));
    ASSERT_GT(max_a, 1e-9) << "vacuous: no deviation developed";
    // The diagonal path skips inactive modes; with runoff + Manning
    // perturbation every mode is active here, so the two must agree.
    EXPECT_LT(maxAbsDiff(diag.a_ensemble, red.a_ensemble), 1e-10 * max_a);
}

TEST(SpectralROM1DReduced, ClearRestoresDiagonalPathExactly) {
    Chain ch(8, 3);
    const auto n = static_cast<std::size_t>(ch.n_active);
    const int k = ch.basis.num_kept;
    std::vector<double> h(n, 1.0), sens(n);
    for (std::size_t i = 0; i < n; ++i) sens[i] = 0.5 + 0.1 * static_cast<double>(i);
    SpectralROM1D a, b;
    makeRom(a, ch.basis, 8, h);
    makeRom(b, ch.basis, 8, h);
    std::vector<double> Mz(static_cast<std::size_t>(k) * static_cast<std::size_t>(k), 0.0);
    b.setReducedOperator(Mz);
    b.clearReducedOperator();
    EXPECT_FALSE(b.hasReducedOperator());
    for (int step = 0; step < 3; ++step) {
        a.advance(20.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
        b.advance(20.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
    }
    EXPECT_EQ(maxAbsDiff(a.a_ensemble, b.a_ensemble), 0.0);
}

TEST(SpectralROM1DReduced, FrozenSourceFixedPointIsOperatorIndependent) {
    // The deviation-form invariant (W3 §"fixed point (mm−1)·b preserved,
    // M-independent"): with a FROZEN sensitivity field both operators settle
    // on δa* = (mm−1)·b exactly. The directional operator therefore cannot
    // change a steady band — it acts on transients only, which is where the
    // H5b leak lives (a pool filling through the window).
    Chain ch(8, 4);
    const auto n = static_cast<std::size_t>(ch.n_active);
    const int k = ch.basis.num_kept;
    std::vector<double> h(n, 50.0), sens(n, 0.0);
    sens[n - 1] = 3.0;
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr_super(static_cast<std::size_t>(ch.n_conduits), 2.0);
    DirectionalOperatorConfig cfg;
    std::vector<double> Mdir;
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr_super.data(), cfg, Mdir));
    SpectralROM1D dir, diag;
    makeRom(dir, ch.basis, 16, h);
    makeRom(diag, ch.basis, 16, h);
    dir.setReducedOperator(Mdir);
    for (int step = 0; step < 60; ++step) {
        dir.advance(60.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
        diag.advance(60.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
    }
    const auto uk = static_cast<std::size_t>(k);
    std::vector<double> b(uk, 0.0);
    for (std::size_t j = 0; j < uk; ++j)
        for (std::size_t i = 0; i < n; ++i) b[j] += ch.basis.P[j * n + i] * sens[i];
    for (std::size_t m = 0; m < 16; ++m)
        for (std::size_t j = 0; j < uk; ++j) {
            const double fixed = (dir.mannings_mult[m] - 1.0) * b[j];
            EXPECT_NEAR(dir.a_ensemble[m * uk + j],  fixed, 1e-9 + 1e-6 * std::fabs(fixed)) << "dir m=" << m << " j=" << j;
            EXPECT_NEAR(diag.a_ensemble[m * uk + j], fixed, 1e-9 + 2e-2 * std::fabs(fixed)) << "diag m=" << m << " j=" << j;
        }
}

TEST(SpectralROM1DReduced, GatedChainCutsTheTransientLeakUpstreamOfAFillingPool) {
    // Sensitivity ramps up at the downstream-most active node only (a pool
    // filling). The symmetric operator relays each increment upstream
    // through its coupling rows; the fully gated operator's upstream rows do
    // not see the pool at all. Compared mid-ramp, upstream of the pool,
    // against the local fixed point both share (a frozen reference at the
    // same instantaneous source, converged).
    Chain ch(7, 3);   // 6 active, k = min(5, 6) = 5 retained: small truncation floor
    const auto n = static_cast<std::size_t>(ch.n_active);
    const int k = ch.basis.num_kept;
    std::vector<double> h(n, 50.0), sens(n, 0.0);
    std::vector<int> sign(static_cast<std::size_t>(ch.n_conduits), 1);
    std::vector<double> fr_super(static_cast<std::size_t>(ch.n_conduits), 2.0);
    std::vector<double> fr_sub(static_cast<std::size_t>(ch.n_conduits), 0.0);
    DirectionalOperatorConfig cfg;
    std::vector<double> Mdir, Msym;
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr_super.data(), cfg, Mdir));
    ASSERT_TRUE(assembleDirectionalReduced(ch.basis.P.data(), k, ch.n_active, ch.n_conduits,
                                           ch.a1.data(), ch.a2.data(), ch.w.data(), sign.data(), fr_sub.data(), cfg, Msym));
    SpectralROM1D dir, sym;
    makeRom(dir, ch.basis, 16, h);
    makeRom(sym, ch.basis, 16, h);
    dir.setReducedOperator(Mdir);
    sym.setReducedOperator(Msym);
    // Slow ramp relative to the modal time constants (K1d·λ_min ≈ 0.05·0.1):
    // 40 steps × 60 s, pool source 0 → 4.
    const int n_steps = 40;
    for (int step = 0; step < n_steps; ++step) {
        sens[n - 1] = 4.0 * (step + 1) / n_steps;
        dir.advance(60.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
        sym.advance(60.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
    }
    std::vector<double> inv(n, 0.0);
    dir.computeQuantiles(h.data(), inv.data());
    sym.computeQuantiles(h.data(), inv.data());
    // Truncation floor: a frozen run at the final source converged to its
    // fixed point (operator-independent, see the previous test).
    SpectralROM1D frozen;
    makeRom(frozen, ch.basis, 16, h);
    for (int step = 0; step < 200; ++step) frozen.advance(60.0, 0.05, h.data(), nullptr, sens.data(), nullptr);
    frozen.computeQuantiles(h.data(), inv.data());

    const double w_dir = dir.q95[0] - dir.q05[0];
    const double w_sym = sym.q95[0] - sym.q05[0];
    const double w_floor = frozen.q95[0] - frozen.q05[0];
    std::printf("[H14 unit] upstream band mid-ramp: directional %.5f  symmetric %.5f  local fixed point %.5f\n",
                w_dir, w_sym, w_floor);
    // Measured while writing this test (not assumed): the symmetric operator
    // does NOT make the upstream band wider mid-ramp here -- its lag makes it
    // NARROWER (0.144 vs the 0.155 local fixed point), because the upstream
    // row is dragged by the pool's still-rising deviation. Either way it is
    // the pool's transient showing up where it cannot physically reach. The
    // property that is actually true, and what H14 delivers, is that the
    // gated upstream node tracks ITS OWN local fixed point regardless of what
    // the pool is doing: measured 0.0005 off vs 0.011 for symmetric (20x).
    ASSERT_GT(std::fabs(w_sym - w_floor), 1e-4) << "symmetric operator must be perturbed by the pool transient (else vacuous)";
    EXPECT_LT(std::fabs(w_dir - w_floor), 0.5 * std::fabs(w_sym - w_floor))
        << "gated upstream node must sit at least twice as close to its local fixed point as the symmetric one";
    EXPECT_GE(w_dir, 0.0);
    // The pool node keeps a band of the same order on both operators.
    const double p_dir = dir.q95[n - 1] - dir.q05[n - 1];
    const double p_sym = sym.q95[n - 1] - sym.q05[n - 1];
    EXPECT_GT(p_dir, 0.5 * p_sym);
    EXPECT_LT(p_dir, 2.0 * p_sym);
}

TEST(SpectralROM1DReduced, RejectsWrongSizeAndUninitialised) {
    Chain ch(7, 3);
    const int k = ch.basis.num_kept;
    SpectralROM1D rom;
    rom.basis = &ch.basis; rom.n_ensemble = 6;
    EXPECT_THROW(rom.setReducedOperator(std::vector<double>(static_cast<std::size_t>(k * k), 0.0)), std::invalid_argument);
    rom.initialize();
    EXPECT_THROW(rom.setReducedOperator(std::vector<double>(static_cast<std::size_t>(k * k + 1), 0.0)), std::invalid_argument);
    EXPECT_NO_THROW(rom.setReducedOperator(std::vector<double>(static_cast<std::size_t>(k * k), 0.0)));
}
