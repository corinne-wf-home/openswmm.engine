/**
 * @file test_soft_forcing_rom.cpp
 * @brief Soft forcing on the 2D ROM; the 1D ROM's member coefficient columns.
 *
 * @ingroup unit_tests
 */

#include <gtest/gtest.h>

#define private public
#include "uncertainty/GraphEigenBasis.hpp"
#include "uncertainty/NetworkLaplacian1D.hpp"
#include "uncertainty/SpectralROM1D.hpp"
#include "uncertainty/LhsShuffle.hpp"
#include "2d/mesh/MeshBuilder.hpp"
#include "2d/uncertainty/MeshEigenBasis.hpp"
#include "2d/uncertainty/SpectralROM.hpp"
#undef private

#include <algorithm>
#include <cmath>
#include <vector>

using namespace openswmm::uncertainty;
using namespace openswmm::twoD;

namespace {

static CsrGraph make_chain_laplacian(int n) {
    int nc = n - 1;
    std::vector<int> n1(static_cast<std::size_t>(nc));
    std::vector<int> n2(static_cast<std::size_t>(nc));
    for (int k = 0; k < nc; ++k) {
        n1[static_cast<std::size_t>(k)] = k;
        n2[static_cast<std::size_t>(k)] = k + 1;
    }
    std::vector<int> is_outfall(static_cast<std::size_t>(n), 0);
    std::vector<int> active_map, full_to_active;
    return NetworkLaplacian1D::buildUniform(n, nc, n1.data(), n2.data(),
                                            is_outfall.data(),
                                            active_map, full_to_active);
}

static MeshData makeStructuredMesh(int N = 3, double domain_m = 3.0) {
    MeshData mesh;
    int nv = (N + 1) * (N + 1);
    int nt = 2 * N * N;
    mesh.resize_vertices(nv);
    mesh.resize_triangles(nt);

    double dx = domain_m / N;
    double dy = domain_m / N;

    for (int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j) {
            int vi = i * (N + 1) + j;
            mesh.vx[static_cast<std::size_t>(vi)] = j * dx;
            mesh.vy[static_cast<std::size_t>(vi)] = i * dy;
            mesh.vz[static_cast<std::size_t>(vi)] = 0.0;
        }

    int t = 0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            int v00 = i * (N + 1) + j,    v01 = i * (N + 1) + j + 1;
            int v10 = (i+1) * (N + 1) + j, v11 = (i+1) * (N + 1) + j + 1;
            mesh.set_triangle(t, v00, v01, v11); ++t;
            mesh.set_triangle(t, v00, v11, v10); ++t;
        }

    buildMeshTopology(mesh);
    return mesh;
}

TEST(SoftForcingROM2D, ProjectsLocAndSpreadAndActivatesModes) {
    MeshData mesh = makeStructuredMesh();
    MeshEigenBasis basis;
    ASSERT_TRUE(basis.build(mesh, 6));

    SpectralROM rom;
    rom.basis = &basis;
    rom.n_ensemble = 11;
    rom.mannings_pert = 0.0;
    rom.rainfall_pert = 0.0;
    rom.initialize();

    std::vector<double> h_det(static_cast<std::size_t>(basis.n_triangles), 0.0);
    std::vector<double> loc(static_cast<std::size_t>(basis.n_triangles));
    std::vector<double> spread(static_cast<std::size_t>(basis.n_triangles));
    const double* mode0 = &basis.P[0];
    for (int i = 0; i < basis.n_triangles; ++i) {
        loc[static_cast<std::size_t>(i)] = mode0[i];
        spread[static_cast<std::size_t>(i)] = 0.25 * mode0[i];
    }

    rom.seed(h_det.data());
    rom.setSoftForcing(loc.data(), spread.data());
    rom.advance(1.0, 0.0, nullptr, nullptr, h_det.data());

    EXPECT_NEAR(rom.r_coarse[0], 1.0, 1.0e-12);
    EXPECT_NEAR(rom.soft_r_spread_[0], 0.25, 1.0e-12);
    EXPECT_GT(rom.n_modes_active, 0);
}

TEST(SoftForcingROM2D, ZeroSpreadEquivalentToNoSoftForcing) {
    MeshData mesh = makeStructuredMesh();
    MeshEigenBasis basis;
    ASSERT_TRUE(basis.build(mesh, 6));

    SpectralROM a;
    a.basis = &basis;
    a.n_ensemble = 11;
    a.mannings_pert = 0.0;
    a.rainfall_pert = 0.0;
    a.initialize();

    SpectralROM c;
    c.basis = &basis;
    c.n_ensemble = 11;
    c.mannings_pert = 0.0;
    c.rainfall_pert = 0.0;
    c.initialize();

    std::vector<double> h_det(static_cast<std::size_t>(basis.n_triangles), 0.0);
    std::vector<double> loc(static_cast<std::size_t>(basis.n_triangles), 0.0);
    std::vector<double> spread(static_cast<std::size_t>(basis.n_triangles), 0.0);

    a.seed(h_det.data());
    a.advance(1.0, 0.0, nullptr, nullptr, h_det.data());
    a.computeQuantiles(h_det.data());

    c.seed(h_det.data());
    c.setSoftForcing(loc.data(), spread.data());
    c.advance(1.0, 0.0, nullptr, nullptr, h_det.data());
    c.computeQuantiles(h_det.data());

    EXPECT_EQ(a.q05, c.q05);
    EXPECT_EQ(a.q50, c.q50);
    EXPECT_EQ(a.q95, c.q95);
}

// ---------------------------------------------------------------------------
// SR-1b family-aware member coefficient (2u-1 for UNIFORM, z=probit(u) else)
// ---------------------------------------------------------------------------

// The 1D ROM no longer carries additive forcings in its modal ODE (H19: they
// ride the full-order additive channel); it only supplies each member's
// coefficient column. These tests pin that column.
static SpectralROM1D makeCoefficientRom() {
    SpectralROM1D rom;
    rom.n_ensemble = 21;            // odd -> one member sits exactly at u = 0.5
    rom.mannings_pert = 0.0;
    rom.runoff_pert = 0.0;
    return rom;
}

TEST(SoftForcingFamily, NormalCoefficientIsProbit) {
    CsrGraph L = make_chain_laplacian(8);
    GraphEigenBasis basis;
    ASSERT_TRUE(basis.build(L, 4));
    SpectralROM1D rom = makeCoefficientRom();
    rom.basis = &basis;
    rom.initialize();

    const auto c = rom.softCoefficients(DistType::NORMAL);
    ASSERT_EQ(static_cast<int>(c.size()), rom.n_ensemble);
    for (int i = 0; i < rom.n_ensemble; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        EXPECT_DOUBLE_EQ(c[ui], probit(rom.soft_u_[ui]));
    }
}

TEST(SoftForcingFamily, LognormalUsesSameProbitCoefficient) {
    CsrGraph L = make_chain_laplacian(8);
    GraphEigenBasis basis;
    ASSERT_TRUE(basis.build(L, 4));
    SpectralROM1D rom = makeCoefficientRom();
    rom.basis = &basis;
    rom.initialize();

    EXPECT_EQ(rom.softCoefficients(DistType::LOGNORMAL),
              rom.softCoefficients(DistType::NORMAL));
}

TEST(SoftForcingFamily, UniformCoefficientIsCenteredBand) {
    CsrGraph L = make_chain_laplacian(8);
    GraphEigenBasis basis;
    ASSERT_TRUE(basis.build(L, 4));
    SpectralROM1D rom = makeCoefficientRom();
    rom.basis = &basis;
    rom.initialize();

    const auto c = rom.softCoefficients(DistType::UNIFORM);
    for (int i = 0; i < rom.n_ensemble; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        EXPECT_DOUBLE_EQ(c[ui], 2.0 * rom.soft_u_[ui] - 1.0);
        EXPECT_GE(c[ui], -1.0);
        EXPECT_LE(c[ui], 1.0);
    }
}

TEST(SoftForcingFamily, ProbitAccuracyOverRange) {
    // Acklam's stated accuracy: |probit(Phi(z)) - z| < 1.2e-8 on z in [-6, 6].
    for (int k = -600; k <= 600; ++k) {
        const double z = 0.01 * k;
        const double phi = 0.5 * std::erfc(-z / std::sqrt(2.0));
        const double zback = probit(phi);
        EXPECT_LT(std::abs(zback - z), 1.2e-8) << "z=" << z;
    }
}

// ---------------------------------------------------------------------------
// SR-3c — lognormal delta-linearization validity (exact vs. delta member field)
//
// The ROM forces the lognormal family with the delta linearization
//   member_i = loc * (1 + z_i * sigma_log)
// while the exact median-preserving field is
//   member_i = loc * exp(z_i * sigma_log).
// The error grows with sigma_log (CV); the SR-3c guard warns once past CV 0.5.
// The one-shot 2D grid warning itself is covered by test_soft_rain_grid_2d
// (LognormalHighCvWarnsOnce); here we bound the numeric accuracy.
// ---------------------------------------------------------------------------

static double lognormalDeltaRelError(double sigma_log, double z) {
    const double f_delta = 1.0 + z * sigma_log;      // loc factored out (loc = 1)
    const double f_exact = std::exp(z * sigma_log);
    return std::abs(f_delta - f_exact) / f_exact;
}

TEST(SoftForcingLognormal, DeltaAccurateAtLowCv) {
    // Interquartile core (|z| <= 1, ~68% of members): CV = 0.2 stays within the
    // ~5% delta bound. exp() convexity makes the low tail worse than the high
    // tail, so the extreme 95th-percentile coefficient is bounded separately.
    const double z68 = 1.0;
    EXPECT_LT(lognormalDeltaRelError(0.2, z68), 0.05);
    EXPECT_LT(lognormalDeltaRelError(0.2, -z68), 0.05);
    // 95th-percentile coefficient the reported q05/q95 bands ride on: still
    // small (< 7.5%) at CV = 0.2.
    const double z95 = 1.6448536269514722;
    EXPECT_LT(lognormalDeltaRelError(0.2, z95), 0.075);
    EXPECT_LT(lognormalDeltaRelError(0.2, -z95), 0.075);
}

TEST(SoftForcingLognormal, DeltaDegradesAtHighCv) {
    // CV = 0.8 is well past the 0.5 guard: the linearization error at z95 is
    // large, documenting why the SR-3c warning fires.
    const double z95 = 1.6448536269514722;
    EXPECT_GT(lognormalDeltaRelError(0.8, z95), 0.15);
}

// ---------------------------------------------------------------------------
// CL-1b — correlated-coherence spatial soft-forcing field (COHERENCE CORR_LEN)
//
// A per-member, per-cell coefficient field W_i[t] replaces the scalar c_i.
// ---------------------------------------------------------------------------

// 2D ROM: a constant-row spatial field (W_i[t] = c_i ∀t) must reproduce the
// comonotone scalar soft path exactly, confirming the shared spatial machinery
// on the 2D SpectralROM.
TEST(SoftForcingCorrelated2D, ComonotoneFieldMatchesScalarPath) {
    MeshData mesh = makeStructuredMesh();
    MeshEigenBasis basis;
    ASSERT_TRUE(basis.build(mesh, 6));

    const int nt = basis.n_triangles;
    std::vector<double> h_det(static_cast<std::size_t>(nt), 0.0);
    std::vector<double> loc(static_cast<std::size_t>(nt));
    std::vector<double> spread(static_cast<std::size_t>(nt));
    const double* mode0 = &basis.P[0];
    const double* mode1 = &basis.P[static_cast<std::size_t>(nt)];
    for (int i = 0; i < nt; ++i) {
        loc[static_cast<std::size_t>(i)]    = mode0[i];
        spread[static_cast<std::size_t>(i)] = 0.3 * mode0[i] + 0.2 * mode1[i];
    }

    SpectralROM a;
    a.basis = &basis;
    a.n_ensemble = 11;
    a.mannings_pert = 0.0;
    a.rainfall_pert = 0.0;
    a.initialize();
    a.seed(h_det.data());
    a.setSoftForcing(loc.data(), spread.data(), DistType::NORMAL);
    a.advance(1.0, 0.0, nullptr, nullptr, h_det.data());
    a.computeQuantiles(h_det.data());

    SpectralROM c;
    c.basis = &basis;
    c.n_ensemble = 11;
    c.mannings_pert = 0.0;
    c.rainfall_pert = 0.0;
    c.initialize();
    c.seed(h_det.data());
    c.setSoftForcing(loc.data(), spread.data(), DistType::NORMAL);
    SpatialUncertaintyField field;
    field.allocate(c.n_ensemble, nt);
    for (int i = 0; i < c.n_ensemble; ++i)
        for (int t = 0; t < nt; ++t)
            field.at(i, t) = c.soft_coeff_[static_cast<std::size_t>(i)];
    c.setSoftForcing(loc.data(), spread.data(), DistType::NORMAL, &field);
    c.advance(1.0, 0.0, nullptr, nullptr, h_det.data());
    c.computeQuantiles(h_det.data());

    ASSERT_EQ(a.n_modes_active, c.n_modes_active);
    for (std::size_t j = 0; j < a.a_ensemble.size(); ++j)
        EXPECT_NEAR(a.a_ensemble[j], c.a_ensemble[j], 1.0e-14) << "flat index " << j;
    for (int t = 0; t < nt; ++t) {
        const auto ut = static_cast<std::size_t>(t);
        EXPECT_NEAR(a.q05[ut], c.q05[ut], 1.0e-14);
        EXPECT_NEAR(a.q50[ut], c.q50[ut], 1.0e-14);
        EXPECT_NEAR(a.q95[ut], c.q95[ut], 1.0e-14);
    }
}

} // anonymous namespace


// ===========================================================================
// PR H7 — true per-family coefficient planes (MIXED + multi-gage)
//
// Before H7 the ROM carried ONE (spread, family) channel, so a MIXED source
// had to borrow a single family's coefficient column for every cell/gage and
// range-match the others by pre-scaling their spread (the SR-4b ~3.0x hack) or
// fall back to the first family outright (SR-1b). H7 adds a second, nullable
// (spread_b, family_b) channel: two projections per advance, each carrying its
// OWN family's per-member coefficient,
//
//     g_ij  +=  c_i * (P^T spread)_j  +  c_i^B * (P^T spread_b)_j
//
// with both coefficient columns drawn from the SAME u_i stream
// (shuffledStrata(M, sample_seed+4)) so comonotone coherence -- one rank per
// member everywhere, the COHERENCE FULL contract -- is preserved across
// families by construction.
//
// NOTE ON WHAT IS AND IS NOT TESTABLE PER CELL. The spec's phrasing ("UNIFORM
// cells' member values uniform ... not z-shaped") cannot be asserted on a
// reconstructed per-CELL value: the ROM's eigenmodes are GLOBAL, so every mode
// sees both planes and every cell's reconstruction is a linear combination of
// both coefficients. That mixing is a property of the spectral ROM itself, not
// something per-family planes can or should remove. The claim that IS exactly
// true -- and is what H7 actually fixes -- is per SOURCE PLANE: each plane's
// modal forcing carries its own family's coefficient column. The marginal-shape
// tests below therefore use spread planes aligned with ORTHOGONAL eigenmodes,
// which isolates one plane per mode and makes the marginal claim exact rather
// than approximate, and the checkerboard fixture then verifies the full
// two-plane decomposition end to end against an independent analytic
// per-member reference.
// ===========================================================================

// Analytic single-step modal solution from a zero initial deviation:
//   d(da)/dt = -rate*da + g,  da(0) = 0  =>  da(dt) = (g/rate)*(1-exp(-rate*dt))
// (Euler below the rate floor). Used as an independent reference so the tests
// do not merely re-derive the ROM's own arithmetic.
static double analyticStepFromZero(double g, double rate, double dt) {
    if (rate > 1.0e-12)
        return (g / rate) * (1.0 - std::exp(-rate * dt));
    return g * dt;
}

// Spearman-style check: two coefficient columns are comonotone when sorting
// members by one also sorts them by the other (identical rank permutation).
static bool sameRankOrder(const std::vector<double>& a,
                          const std::vector<double>& b) {
    const std::size_t M = a.size();
    if (b.size() != M) return false;
    std::vector<std::size_t> ia(M), ib(M);
    for (std::size_t i = 0; i < M; ++i) { ia[i] = i; ib[i] = i; }
    std::sort(ia.begin(), ia.end(), [&](std::size_t p, std::size_t q) { return a[p] < a[q]; });
    std::sort(ib.begin(), ib.end(), [&](std::size_t p, std::size_t q) { return b[p] < b[q]; });
    return ia == ib;
}

// ---------------------------------------------------------------------------
// Coefficient-column contract
// ---------------------------------------------------------------------------

TEST(SoftForcingFamilyPlanes, BothFamiliesShareOneRankPerMember) {
    CsrGraph L = make_chain_laplacian(8);
    GraphEigenBasis basis;
    ASSERT_TRUE(basis.build(L, 4));
    SpectralROM1D rom = makeCoefficientRom();
    rom.basis = &basis;
    rom.initialize();

    const auto cA = rom.softCoefficients(DistType::NORMAL);
    const auto cB = rom.softCoefficients(DistType::UNIFORM);

    // COHERENCE FULL contract: one rank per member across families.
    EXPECT_TRUE(sameRankOrder(cA, cB));

    // Both columns are zero-mean, so neither family shifts the median.
    double mean_a = 0.0, mean_b = 0.0;
    for (int i = 0; i < rom.n_ensemble; ++i) {
        mean_a += cA[static_cast<std::size_t>(i)];
        mean_b += cB[static_cast<std::size_t>(i)];
    }
    EXPECT_NEAR(mean_a / rom.n_ensemble, 0.0, 1.0e-12);
    EXPECT_NEAR(mean_b / rom.n_ensemble, 0.0, 1.0e-12);

    // The nominal member is nominal in BOTH families simultaneously.
    int nominal = -1;
    for (int i = 0; i < rom.n_ensemble; ++i)
        if (std::abs(cB[static_cast<std::size_t>(i)]) < 1.0e-15)
            nominal = i;
    ASSERT_GE(nominal, 0) << "odd M must place one member exactly at u = 0.5";
    EXPECT_NEAR(cA[static_cast<std::size_t>(nominal)], 0.0, 1.0e-12);
}

// ---------------------------------------------------------------------------
// 2D: the checkerboard MIXED grid the SR-4b hack was written for. Cells
// alternate NORMAL/UNIFORM family; each family's cells go into their own plane
// (zero elsewhere), and the result is checked against an independent analytic
// per-member reference built from the two coefficient columns.
// ---------------------------------------------------------------------------

TEST(SoftForcingFamilyPlanes2D, CheckerboardMixedGridUsesBothFamilies) {
    MeshData mesh = makeStructuredMesh(4, 4.0);
    MeshEigenBasis basis;
    ASSERT_TRUE(basis.build(mesh, 6));

    const int nt = basis.n_triangles;
    const auto unt = static_cast<std::size_t>(nt);
    std::vector<double> h_det(unt, 0.0);
    std::vector<double> spread_norm(unt, 0.0), spread_unif(unt, 0.0);

    // /family_code checkerboard: even cells NORMAL, odd cells UNIFORM. Each
    // plane holds only its own cells' spread; the two are disjoint and sum to
    // the full spread field, which is exactly the split H7 introduces.
    const double* mode0 = &basis.P[0];
    for (int t = 0; t < nt; ++t) {
        const auto ut = static_cast<std::size_t>(t);
        const double sp = 0.05 + 0.2 * std::abs(mode0[t]);
        if (t % 2 == 0) spread_norm[ut] = sp;
        else            spread_unif[ut] = sp;
    }

    SpectralROM rom;
    rom.basis = &basis;
    rom.n_ensemble = 21;
    rom.mannings_pert = 0.0;
    rom.rainfall_pert = 0.0;
    rom.mode_drop_threshold = 0.0;
    rom.initialize();
    rom.seed(h_det.data());
    rom.setSoftForcing(nullptr, spread_norm.data(), DistType::NORMAL, nullptr,
                       spread_unif.data(), DistType::UNIFORM);

    ASSERT_TRUE(rom.softPlanesError().empty()) << rom.softPlanesError();
    ASSERT_EQ(static_cast<int>(rom.softCoeffB().size()), rom.n_ensemble);
    EXPECT_TRUE(sameRankOrder(rom.softCoeff(), rom.softCoeffB()));

    const double dt = 1.0;
    rom.advance(dt, 0.0, nullptr, nullptr, h_det.data());
    rom.computeQuantiles(h_det.data());

    // Independent reference: project each plane by hand off the basis and
    // integrate the modal ODE in closed form, per member.
    const auto nk = static_cast<std::size_t>(rom.n_kept);
    std::vector<double> RA(nk, 0.0), RB(nk, 0.0);
    for (std::size_t j = 0; j < nk; ++j) {
        const double* Pj = &basis.P[j * unt];
        for (std::size_t t = 0; t < unt; ++t) {
            RA[j] += Pj[t] * spread_norm[t];
            RB[j] += Pj[t] * spread_unif[t];
        }
    }
    // Both families genuinely reach the modal forcing.
    double maxA = 0.0, maxB = 0.0;
    for (std::size_t j = 0; j < nk; ++j) {
        maxA = std::max(maxA, std::abs(RA[j]));
        maxB = std::max(maxB, std::abs(RB[j]));
    }
    EXPECT_GT(maxA, 1.0e-6);
    EXPECT_GT(maxB, 1.0e-6);

    for (int i = 0; i < rom.n_ensemble; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        const double cA = rom.softCoeff()[ui];
        const double cB = rom.softCoeffB()[ui];
        for (std::size_t j = 0; j < nk; ++j) {
            // K_eff = 0 here (no Manning channel), so rate = 0 -> Euler branch.
            const double expected = analyticStepFromZero(cA * RA[j] + cB * RB[j],
                                                         0.0, dt);
            EXPECT_NEAR(rom.a_ensemble[ui * nk + j], expected, 1.0e-13)
                << "member " << i << " mode " << j;
        }
    }

    // A real band, and the median still tracks the deterministic field.
    double widest = 0.0;
    for (std::size_t t = 0; t < unt; ++t)
        widest = std::max(widest, rom.q95[t] - rom.q05[t]);
    EXPECT_GT(widest, 1.0e-6);
    for (std::size_t t = 0; t < unt; ++t)
        EXPECT_NEAR(rom.q50[t], 0.0, 1.0e-12) << "cell " << t;
}

// 2D regression lock, mirroring the 1D one.
TEST(SoftForcingFamilyPlanes2D, NullAndZeroSecondPlaneAreBitIdentical) {
    MeshData mesh = makeStructuredMesh();
    MeshEigenBasis basis;
    ASSERT_TRUE(basis.build(mesh, 6));

    const int nt = basis.n_triangles;
    const auto unt = static_cast<std::size_t>(nt);
    std::vector<double> h_det(unt, 0.0), loc(unt), spread(unt);
    std::vector<double> zero_plane(unt, 0.0);
    const double* mode0 = &basis.P[0];
    const double* mode1 = &basis.P[unt];
    for (std::size_t t = 0; t < unt; ++t) {
        loc[t]    = mode0[t];
        spread[t] = 0.3 * mode0[t] + 0.2 * mode1[t];
    }

    auto run = [&](const double* second) {
        SpectralROM rom;
        rom.basis = &basis;
        rom.n_ensemble = 11;
        rom.mannings_pert = 0.0;
        rom.rainfall_pert = 0.0;
        rom.initialize();
        rom.seed(h_det.data());
        rom.setSoftForcing(loc.data(), spread.data(), DistType::NORMAL,
                           nullptr, second, DistType::UNIFORM);
        rom.advance(1.0, 0.0, nullptr, nullptr, h_det.data());
        rom.computeQuantiles(h_det.data());
        return rom;
    };

    SpectralROM none = run(nullptr);
    SpectralROM zero = run(zero_plane.data());

    ASSERT_EQ(none.n_modes_active, zero.n_modes_active);
    for (std::size_t i = 0; i < none.a_ensemble.size(); ++i)
        EXPECT_EQ(none.a_ensemble[i], zero.a_ensemble[i]) << "flat index " << i;
    for (std::size_t t = 0; t < unt; ++t) {
        EXPECT_EQ(none.q05[t], zero.q05[t]);
        EXPECT_EQ(none.q50[t], zero.q50[t]);
        EXPECT_EQ(none.q95[t], zero.q95[t]);
    }
}
