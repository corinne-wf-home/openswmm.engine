/**
 * @file mc_quantiles.hpp
 * @brief Monte-Carlo quantile / interval-coverage estimators for the ROM-vs-MC
 *        regression gates.
 *
 * Added 2026-10-04 after a cross-vendor review of PR SP3 found two metric
 * defects in the soft-rain gates (shared with the PR-10 harness idiom):
 *
 *  1. "q95 - q05" of the 21-member MC reference was computed as h[19] - h[1],
 *     which spans plotting positions 0.071..0.929, not 0.05..0.95. For a normal
 *     spread that is ~11% narrower than the true 5-95 span, biasing every
 *     reported ROM/MC width ratio HIGH by ~12%.
 *  2. "coverage" counted only whether the MC MEDIAN fell inside the ROM band.
 *     That is median containment, not interval coverage: it says nothing about
 *     whether ~90% of outcomes fall inside the ROM's q05-q95 interval.
 *
 * Convention. The ROMs take their q05/q95 at the member whose LHS stratum
 * midpoint u_i = (i+0.5)/M equals 0.05 / 0.95 exactly (idx = round(p*(M-1))).
 * The consistent estimator for an N-member LHS reference is therefore the
 * midpoint plotting position p_i = (i+0.5)/N with linear interpolation between
 * neighbours (Hazen / Hyndman-Fan type 5). With N = 21 and p = 0.05 that is
 * 0.55 of the way from h[0] to h[1]; the old h[1] was the 0.071 point.
 */
#ifndef OPENSWMM_TESTS_MC_QUANTILES_HPP
#define OPENSWMM_TESTS_MC_QUANTILES_HPP

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <vector>

namespace mcq {

/// Quantile at probability p of an ASCENDING-sorted sample, midpoint plotting
/// positions p_i = (i+0.5)/N, linear interpolation, clamped to the extremes.
inline double quantileMidpoint(const std::vector<double>& sorted, double p) {
    const std::size_t n = sorted.size();
    if (n == 0) return 0.0;
    if (n == 1) return sorted[0];
    const double pos = p * static_cast<double>(n) - 0.5;   // fractional index
    if (pos <= 0.0) return sorted.front();
    if (pos >= static_cast<double>(n - 1)) return sorted.back();
    const auto i = static_cast<std::size_t>(pos);
    const double w = pos - static_cast<double>(i);
    return sorted[i] + w * (sorted[i + 1] - sorted[i]);
}

/// Fraction of sample values inside the closed interval [lo, hi]. For a
/// correct 5-95 band this should be ~0.90 (19 of 21 strata midpoints lie in
/// u in [0.05, 0.95]: 0.905).
inline double intervalCoverage(const std::vector<double>& values, double lo, double hi) {
    if (values.empty()) return 0.0;
    std::size_t inside = 0;
    for (double v : values) if (v >= lo && v <= hi) ++inside;
    return static_cast<double>(inside) / static_cast<double>(values.size());
}

/// C1 (owner, 2026-10-05): a cell meets the FLOOR when its empirical member
/// coverage is >= kCalibratedFloor; otherwise it is RANKING ONLY -- still
/// registered, still printing, documented with its fix PR, never loosened.
/// C2 (owner, 2026-10-06, option 1): calibration is two-sided. A cell that
/// meets the floor is CALIBRATED only if its median ROM/MC width ratio is
/// also <= kCalibratedWidthCeiling; above that it is CONSERVATIVE (over-wide):
/// enough outcomes inside, but by being too wide. Conservative is a label,
/// not a failure. Cells that claim "validated" assert the floor; the ceiling
/// is reported and documented, uniformly, with no per-test variation.
constexpr double kCalibratedFloor        = 0.80;
constexpr double kCalibratedWidthCeiling = 1.50;

enum class Calibration { RankingOnly, Conservative, Calibrated };

inline Calibration classifyCalibration(double member_cov, double ratio_med) {
    if (member_cov < kCalibratedFloor) return Calibration::RankingOnly;
    if (ratio_med > kCalibratedWidthCeiling) return Calibration::Conservative;
    return Calibration::Calibrated;
}

/// Prints the C1/C2 verdict line for a cell. Returns whether the FLOOR is met
/// (the gated part); a Conservative cell returns true.
inline bool reportCalibration(const char* cell, double member_cov, double ratio_med,
                              const char* fix_pr) {
    const Calibration c = classifyCalibration(member_cov, ratio_med);
    const char* verdict =
        c == Calibration::Calibrated   ? "CALIBRATED (coverage >= 0.80, width <= 1.5)" :
        c == Calibration::Conservative ? "CONSERVATIVE (coverage >= 0.80 but over-wide, width > 1.5)" :
                                         "RANKING ONLY (coverage < 0.80)";
    std::printf("[C1/C2] %-32s member-coverage=%.3f width-med=%.3f -> %s%s%s\n",
                cell, member_cov, ratio_med, verdict,
                c == Calibration::Calibrated ? "" : "; fix: ",
                c == Calibration::Calibrated ? "" : fix_pr);
    return c != Calibration::RankingOnly;
}

}  // namespace mcq

#endif
