/**
 * @file test_rom_surcharge_attenuation.cpp
 * @brief PR H5 — unit tests for the pure computeSurchargeAlpha() function.
 *
 * @details Hand-built fixtures, independent of the real hydraulics: these
 *          pin the exact ramp formula and its degenerate/guard behavior.
 *          Engine-level plumbing coverage (crown_elev/invert_elev/head
 *          reaching the function through the real engine) lives in
 *          test_engine_1d_rom_lifecycle.cpp; the effect on the ROM's
 *          projected sensitivity lives in the SurchargeAttenuation suite of
 *          test_spectral_rom1d.cpp.
 *
 * @ingroup engine_uncertainty
 */

#include <gtest/gtest.h>

#include <vector>

#include "uncertainty/RomSurchargeAttenuation.hpp"

using openswmm::uncertainty::computeSurchargeAlpha;
using openswmm::uncertainty::SurchargeAttenuationConfig;

namespace {

// Single-node convenience wrapper: builds full-node-space arrays for one
// active node at full index 0, with default ramp band, and returns alpha[0].
double alphaFor(double crown_elev, double invert_elev, double head,
                 SurchargeAttenuationConfig cfg = {}) {
    const std::vector<double> crown  = {crown_elev};
    const std::vector<double> invert = {invert_elev};
    const std::vector<double> heads  = {head};
    const std::vector<int> active_to_full = {0};
    std::vector<double> alpha(1, -1.0);
    computeSurchargeAlpha(1, active_to_full.data(), crown.data(),
                           invert.data(), heads.data(), 1, cfg,
                           alpha.data());
    return alpha[0];
}

}  // namespace

TEST(RomSurchargeAttenuation, WellBelowRampLoIsAlphaFree) {
    // Free surface: the normal-depth elasticity 0.6 (H5b), not 1.
    EXPECT_EQ(alphaFor(/*crown=*/2.0, /*invert=*/0.0, /*head=*/1.0), 0.6);
}

TEST(RomSurchargeAttenuation, WellAboveRampHiReachesAlphaSurcharged) {
    // Deep surcharge: the pressurised-friction elasticity 2.0 (H5b). Note the
    // direction: H5 shipped a floor of 0.05 here; the per-node breakdown
    // showed that attenuated exactly the nodes that were under-predicted.
    EXPECT_NEAR(alphaFor(/*crown=*/2.0, /*invert=*/0.0, /*head=*/4.0), 2.0, 1e-12);
}

TEST(RomSurchargeAttenuation, ExactlyAtRampLoIsAlphaFree) {
    EXPECT_EQ(alphaFor(/*crown=*/1.0, /*invert=*/0.0, /*head=*/0.9), 0.6);
}

TEST(RomSurchargeAttenuation, ExactlyAtRampHiReachesAlphaSurcharged) {
    EXPECT_NEAR(alphaFor(/*crown=*/1.0, /*invert=*/0.0, /*head=*/1.1), 2.0, 1e-12);
}

TEST(RomSurchargeAttenuation, MidpointRatioGivesExactMidpoint) {
    // ratio 1.0 is the midpoint of [0.9, 1.1] -> halfway from 0.6 to 2.0 = 1.3.
    EXPECT_NEAR(alphaFor(/*crown=*/1.0, /*invert=*/0.0, /*head=*/1.0), 1.3, 1e-12);
}

TEST(RomSurchargeAttenuation, DegenerateCrownDepthFailsToFreeSurfaceLaw) {
    // No classifiable crown (outfalls): the free-surface elasticity, the
    // regime such nodes are in. (H5 failed open to 1.0.)
    EXPECT_EQ(alphaFor(/*crown=*/0.0, /*invert=*/0.0, /*head=*/50.0), 0.6);
    EXPECT_EQ(alphaFor(/*crown=*/1.0, /*invert=*/2.0, /*head=*/50.0), 0.6);
}

TEST(RomSurchargeAttenuation, OutOfRangeActiveIndexFailsOpenToOne) {
    // A full-node index outside [0, n_nodes_full) must not read out of
    // bounds and must fail open (alpha = 1.0), same convention as
    // RomDiagTrust's out-of-range defensiveness.
    const std::vector<double> crown  = {2.0};
    const std::vector<double> invert = {0.0};
    const std::vector<double> heads  = {10.0};  // would be ratio 5.0 if read
    const std::vector<int> active_to_full = {99};
    std::vector<double> alpha(1, -1.0);
    computeSurchargeAlpha(1, active_to_full.data(), crown.data(),
                           invert.data(), heads.data(), 1,
                           SurchargeAttenuationConfig{}, alpha.data());
    EXPECT_EQ(alpha[0], 0.6);   // out of range -> free-surface law (H5b)
}

TEST(RomSurchargeAttenuation, CustomRampBandChangesOutput) {
    // Same physical state (ratio = 1.0) under two different configured
    // bands must give two different results -- proves the formula reads
    // cfg.ramp_lo/ramp_hi, not hardcoded 0.9/1.1.
    SurchargeAttenuationConfig narrow{/*ramp_lo=*/0.95, /*ramp_hi=*/1.05};
    SurchargeAttenuationConfig wide{/*ramp_lo=*/0.5, /*ramp_hi=*/1.5};
    const double a_narrow = alphaFor(1.0, 0.0, 1.0, narrow);
    const double a_wide   = alphaFor(1.0, 0.0, 1.0, wide);
    EXPECT_NEAR(a_narrow, 1.3, 1e-12);  // still the midpoint of its own band (0.6 .. 2.0)
    EXPECT_NEAR(a_wide, 1.3, 1e-12);
    // Push ratio to 1.1: narrow band has already reached its floor, wide
    // band's ramp has not (still strictly above the floor).
    const double a_narrow_at_1p1 = alphaFor(1.0, 0.0, 1.1, narrow);
    const double a_wide_at_1p1   = alphaFor(1.0, 0.0, 1.1, wide);
    EXPECT_NEAR(a_narrow_at_1p1, narrow.alpha_surcharged, 1e-12);
    EXPECT_LT(a_wide_at_1p1, wide.alpha_surcharged);   // still on the ramp in the wide band
}

TEST(RomSurchargeAttenuation, ElasticitiesAreConfigurable) {
    // Both plateaus read the config, not hardcoded constants.
    SurchargeAttenuationConfig cfg;
    cfg.alpha_free = 1.0;        // the pre-H5b "elasticity one" fixed point
    cfg.alpha_surcharged = 0.05; // H5's old attenuation floor, reproducible on request
    EXPECT_NEAR(alphaFor(/*crown=*/2.0, /*invert=*/0.0, /*head=*/1.0, cfg), 1.0, 1e-12);
    EXPECT_NEAR(alphaFor(/*crown=*/2.0, /*invert=*/0.0, /*head=*/4.0, cfg), 0.05, 1e-12);
}

TEST(RomSurchargeAttenuation, MultipleNodesAreIndependent) {
    // Three active nodes in one call, mixed regimes: free-surface,
    // surcharged, and mid-ramp. Each slot must reflect only its own node's
    // state -- no cross-contamination between array slots.
    const std::vector<double> crown  = {2.0, 2.0, 2.0};
    const std::vector<double> invert = {0.0, 0.0, 0.0};
    // ratios: 0.5 (free surface), 2.0 (surcharged), 1.0 (midpoint).
    const std::vector<double> heads  = {1.0, 4.0, 2.0};
    const std::vector<int> active_to_full = {0, 1, 2};
    std::vector<double> alpha(3, -1.0);
    computeSurchargeAlpha(3, active_to_full.data(), crown.data(),
                           invert.data(), heads.data(), 3,
                           SurchargeAttenuationConfig{}, alpha.data());
    EXPECT_EQ(alpha[0], 0.6);
    EXPECT_NEAR(alpha[1], 2.0, 1e-12);   // deep surcharge -> the pressurised elasticity
    EXPECT_NEAR(alpha[2], 1.3, 1e-12);   // ramp midpoint
}
