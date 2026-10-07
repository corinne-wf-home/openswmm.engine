/**
 * @file RomSurchargeAttenuation.hpp
 * @brief PR H5 / H5b — per-active-node Manning-sensitivity ELASTICITY factor
 *        `alpha_n` for the 1D ROM: 0.6 in free-surface flow, 2.0 under
 *        surcharge, ramped across the crown.
 *
 * @details The modal Manning-sensitivity source term
 *          `-lam_j*K1d*(1/mm-1)*b_j` (SpectralROM1D::advance) has the fixed
 *          point `delta_a = (mm-1)*b_j`, i.e. it makes the local head respond
 *          to Manning's n with elasticity ONE: dh/h = dn/n. Neither flow
 *          regime does that.
 *
 *          - Free surface: normal depth y ~ (n*Q)^(3/5) for a wide section,
 *            elasticity 0.6. Elasticity 1 over-predicts by 1/0.6 = 1.67x,
 *            which is the ~1.6-2.1x over-width measured on every
 *            free-surface node (P8: 1.96; the H5 chain's J1-J2: 2.1).
 *          - Pressurised: the friction head loss h_f ~ n^2 (Manning full-pipe
 *            or Darcy-Weisbach with Manning n), elasticity 2. Elasticity 1
 *            UNDER-predicts by 2x, which is what the H5 chain's chokepoint
 *            nodes showed once measured per node with no attenuation at all
 *            (J5: 5.0 ft vs MC 10.0 ft).
 *
 *          `alpha_n` folds elementwise into the sensitivity reference BEFORE
 *          projection (`b_j = P[:,j]^T*(alpha ⊙ bref)`), scaling only the
 *          Manning-sensitivity channel; the forcing-sensitivity channel
 *          (runoff / soft rain) is projected separately and untouched.
 *
 *          History (keep, it is the point of the per-node metric): PR H5
 *          (2026-08) shipped this as an ATTENUATION ramping 1 -> 0.05 above
 *          crown, calibrated to make the surcharged chain's MEDIAN width
 *          ratio ~1.0. The P7 audit (2026-10-04) showed that cell's member
 *          coverage was 0.51, and the H5b per-node breakdown (2026-10-06)
 *          showed why: the median averaged free-surface nodes at ~2x
 *          over-wide with surcharged nodes crushed to ~0.01x -- the floor was
 *          attenuating exactly the nodes that were already under-predicted.
 *          With alpha == 1 everywhere coverage was 0.74; with the chokepoint
 *          nodes at 2 it was 0.86. The pre-H5 "50-190x over-prediction" the
 *          attenuation was built against came from the pre-PR-10
 *          absolute-head projection and did not exist on this fixture once
 *          depth was the sensitivity reference. See VALIDATION.md, "H5b".
 *
 *          Node-level formulation (unchanged from H5): the engine exposes
 *          crown data as a per-NODE scalar (`ctx_.nodes.crown_elev`), the same
 *          quantity DWSolver uses for `HSnapshot::node_surcharged`. The factor
 *          is a smooth ramp in `ratio_n = (head-invert)/(crown_elev-invert)`
 *          between `ramp_lo` and `ramp_hi`.
 *
 * @ingroup engine_uncertainty
 */
#ifndef OPENSWMM_ENGINE_UNCERTAINTY_ROM_SURCHARGE_ATTENUATION_HPP
#define OPENSWMM_ENGINE_UNCERTAINTY_ROM_SURCHARGE_ATTENUATION_HPP

#include <algorithm>

namespace openswmm::uncertainty {

/// Ramp band for the surcharge attenuation factor, as config fields (not
/// hardcoded constants) so fixtures/tests can exercise non-default bands.
struct SurchargeAttenuationConfig {
    /// depth/crown_depth ratio at and below which alpha == alpha_free.
    double ramp_lo = 0.9;
    /// depth/crown_depth ratio at and above which alpha == alpha_surcharged.
    double ramp_hi = 1.1;
    /// Elasticity of local head to Manning's n in FREE-SURFACE flow. Normal
    /// depth goes as (n*Q)^(3/5): 0.6. The ROM's fixed point assumes 1, so this
    /// is the factor that removes the ~1.67x free-surface over-prediction
    /// (H13). Physical constant, not a fitted dial.
    double alpha_free = 0.6;
    /// Elasticity of local head to Manning's n under SURCHARGE. Pressurised
    /// friction head goes as n^2: 2.0. Replaces H5's attenuation floor (0.05),
    /// which had the wrong sign; see the file comment. Physical constant.
    double alpha_surcharged = 2.0;
};

/**
 * @brief Per-active-node Manning-sensitivity elasticity factor.
 *
 * `alpha_n = alpha_free + (alpha_surcharged - alpha_free) * clamp((ratio_n -
 * ramp_lo) / (ramp_hi - ramp_lo), 0, 1)`, with `ratio_n = (head_n -
 * invert_n) / (crown_elev_n - invert_n)`. alpha_free well below crown,
 * alpha_surcharged well above, linear in between.
 *
 * Degenerate guard: a node whose crown_depth (`crown_elev - invert_elev`) is
 * <= 1e-6 (outfalls, or any node without a classifiable crown), and any
 * out-of-range index in active_to_full, cannot be ramped and gets alpha ==
 * alpha_free -- the free-surface law, the regime such nodes are in.
 *
 * @param n_active       Length of active_to_full and alpha_out.
 * @param active_to_full Active-index -> full-node-index map, length n_active.
 * @param crown_elev     Full-node-space crown elevation, length n_nodes_full.
 * @param invert_elev    Full-node-space invert elevation, length n_nodes_full.
 * @param head           Full-node-space current head, length n_nodes_full.
 * @param n_nodes_full   Length of crown_elev/invert_elev/head.
 * @param cfg            Ramp band + the two elasticities.
 * @param alpha_out      Output buffer, length n_active.
 */
inline void computeSurchargeAlpha(int n_active, const int* active_to_full,
                                   const double* crown_elev,
                                   const double* invert_elev,
                                   const double* head, int n_nodes_full,
                                   const SurchargeAttenuationConfig& cfg,
                                   double* alpha_out) noexcept {
    const double band = cfg.ramp_hi - cfg.ramp_lo;
    for (int ai = 0; ai < n_active; ++ai) {
        const int ui = active_to_full[ai];
        if (ui < 0 || ui >= n_nodes_full) {
            alpha_out[ai] = cfg.alpha_free;
            continue;
        }
        const double crown_depth = crown_elev[ui] - invert_elev[ui];
        if (crown_depth <= 1.0e-6) {
            alpha_out[ai] = cfg.alpha_free;
            continue;
        }
        const double depth = head[ui] - invert_elev[ui];
        const double ratio = depth / crown_depth;
        // s = 0 at/below ramp_lo (free surface), 1 at/above ramp_hi (surcharged)
        double s = (band > 0.0) ? (ratio - cfg.ramp_lo) / band : (ratio < cfg.ramp_lo ? 0.0 : 1.0);
        s = std::clamp(s, 0.0, 1.0);
        alpha_out[ai] = cfg.alpha_free + (cfg.alpha_surcharged - cfg.alpha_free) * s;
    }
}

}  // namespace openswmm::uncertainty

#endif  // OPENSWMM_ENGINE_UNCERTAINTY_ROM_SURCHARGE_ATTENUATION_HPP
