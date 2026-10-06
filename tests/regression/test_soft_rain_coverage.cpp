/**
 * @file test_soft_rain_coverage.cpp
 * @brief SR-5 — soft-rainfall ROM bands vs brute-force Monte Carlo.
 *
 * Reference: 21 deterministic engine runs (no ROM), each with the rain-gage
 * input scaled by the materialized member `1 + z(u_stratum)·CV` — the same
 * NORMAL location-scale prior the soft-rain ROM propagates. ROM run: identical
 * network with `[SOFT_RAINGAGES] RG1 NORMAL CV 0.20`, M = 50.
 *
 * Assertions (checklist floors, tightened toward the first measured run whose
 * actuals are documented inline):
 *   (a) coverage — ROM [q05, q95] contains the MC median at report times
 *       t > 60 s for at least the stated fraction of (node, time) samples;
 *   (b) width ratio — ROM (q95−q05) vs MC (q95−q05) stays within [0.3x, 3x]
 *       for at least the stated fraction of samples where the MC width is
 *       resolvable.
 *
 * The network uses storage nodes with a large constant surface area so heads
 * rise gradually through the whole run (transient dh/dt), keeping the
 * rainfall-rate-driven soft spread active at every report boundary and giving
 * the MC members resolvable head spread.
 *
 * Runtime: 22 engine runs on a 6-node, 1-hour network — a few seconds total.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "uncertainty/SpectralROM1D.hpp"
#include "uncertainty/LhsShuffle.hpp"
#include "mc_quantiles.hpp"

#if defined(_WIN32)
#  include <process.h>   // _getpid
#  define getpid _getpid
#else
#  include <unistd.h>    // getpid
#endif

namespace {

constexpr double kBaseRainInHr = 1.0;     // location rain rate (in/hr)
// Sweep hooks (2026-10-04): SR5_CV overrides the CV, SR5_MEMBERS the ROM
// ensemble size (via a [2D_ROM] MEMBERS line, which buildROM1D inherits).
// SR5_TRACE=1 prints every (node, time) width ratio. Defaults reproduce the
// gate exactly; the hooks exist for the calibration sweeps recorded in
// VALIDATION.md ("Interval-coverage audit") and are not used by ctest.
inline double envDouble(const char* name, double dflt) {
    const char* v = std::getenv(name); return (v && *v) ? std::atof(v) : dflt;
}
inline int envInt(const char* name, int dflt) {
    const char* v = std::getenv(name); return (v && *v) ? std::atoi(v) : dflt;
}
static const double kCV        = envDouble("SR5_CV", 0.20);    // NORMAL coefficient of variation
static const int    kMembers   = envInt("SR5_MEMBERS", 0);      // 0 = engine default (50)
static const double kK1dScale  = envDouble("SR5_K1D_SCALE", 1.0); // SR-6 falsification knob
static const double kPipeDiam  = envDouble("SR5_PIPE_DIAM", 0.5);  // SR-6: surcharge-elasticity probe (ft)
static const double kSubWidth  = envDouble("SR5_SUBCATCH_WIDTH", 5.0); // SR-6: runoff-elasticity probe (ft)
static const double kSubArea   = envDouble("SR5_SUBCATCH_AREA", 5.0);  // SR-6: equilibrium probe (acres)
static const int    kElasOn    = envInt("SR5_RUNOFF_ELASTICITY", 1);    // SR-6 A/B: 0 = pre-SR-6 mapping
constexpr int    kMcRuns       = 21;
constexpr double kReportStep   = 300.0;   // s (5 min)
constexpr double kEndTime      = 3600.0;  // s (1 h sampling window)

// Build the fixture .inp. `rain_mult` scales the gage timeseries (the
// materialized MC member); `with_soft` adds the [SOFT_RAINGAGES] entry that
// activates the soft-rain 1D ROM.
std::string fixtureInp(double rain_mult, bool with_soft) {
    const double rain = kBaseRainInHr * rain_mult;
    std::ostringstream f;
    f << std::setprecision(10);
    f << "[TITLE]\nSR-5 soft-rain ROM vs MC coverage\n\n";
    f << "[OPTIONS]\n";
    f << "FLOW_UNITS CFS\nINFILTRATION HORTON\nFLOW_ROUTING DYNWAVE\n";
    f << "START_DATE 01/01/2025\nSTART_TIME 00:00:00\n";
    f << "REPORT_START_DATE 01/01/2025\nREPORT_START_TIME 00:00:00\n";
    f << "END_DATE 01/01/2025\nEND_TIME 01:05:00\n";
    f << "REPORT_STEP 00:05:00\nWET_STEP 00:01:00\nDRY_STEP 00:05:00\nROUTING_STEP 00:00:30\n";
    f << "MINIMUM_STEP 0.5\nTHREADS 1\n\n";
    f << "[EVAPORATION]\nCONSTANT 0.0\n\n";
    f << "[RAINGAGES]\nRG1 INTENSITY 0:05 1.0 TIMESERIES TS0\n\n";
    // Constant rain at 5-min recording interval so the rate is flat over the run.
    f << "[TIMESERIES]\n";
    for (int m = 0; m <= 65; m += 5) {
        const int hh = m / 60, mm = m % 60;
        f << "TS0 01/01/2025 " << std::setfill('0') << std::setw(2) << hh
          << ":" << std::setfill('0') << std::setw(2) << mm << " " << rain << "\n";
    }
    f << std::setfill(' ') << "\n";
    f << "[SUBCATCHMENTS]\n";
    for (int i = 1; i <= 5; ++i)
        f << "S" << i << " RG1 J" << i << " " << kSubArea << " " << kSubWidth << " 100 1.0 0\n";
    f << "\n[SUBAREAS]\n";
    for (int i = 1; i <= 5; ++i)
        f << "S" << i << " 0.015 0.20 0.00 0.00 100 OUTLET\n";
    f << "\n[INFILTRATION]\n";
    for (int i = 1; i <= 5; ++i)
        f << "S" << i << " 0.0 0.0 0.0 0.0 0.0\n";
    // Large-area storage nodes → heads rise gradually (transient dh/dt).
    f << "\n[STORAGE]\n";
    for (int i = 1; i <= 5; ++i)
        f << "J" << i << " " << (100 - i) << " 60 0 FUNCTIONAL 0 0 20000\n";
    f << "\n[OUTFALLS]\nO1 90 FREE NO\n\n";
    f << "[CONDUITS]\n";
    for (int i = 1; i <= 4; ++i)
        f << "C" << i << " J" << i << " J" << (i + 1) << " 200 0.02 0 0 0 0\n";
    f << "C5 J5 O1 200 0.02 0 0 0 0\n\n";
    f << "[XSECTIONS]\n";
    for (int i = 1; i <= 5; ++i)
        f << "C" << i << " CIRCULAR " << kPipeDiam << " 0 0 0 1\n";
    f << "\n[REPORT]\nINPUT NO\nCONTINUITY YES\nNODES ALL\nLINKS ALL\n\n";
    if (with_soft) {
        f << "[SOFT_RAINGAGES]\nRG1 NORMAL CV " << kCV << "\n";
        if (kMembers >= 2) f << "\n[2D_ROM]\nMEMBERS " << kMembers << "\n";
    }
    return f.str();
}

struct RunResult {
    std::map<std::string, std::vector<double>> heads;
    std::map<std::string, std::vector<double>> q05, q50, q95;
    std::map<std::string, std::vector<double>> latq;   // lateral inflow (cfs) at report times
    std::vector<double> times;
    bool ok = false;
};

RunResult runCase(const std::string& inp_text, const char* tag, bool with_soft) {
    RunResult out;
    // Per-process prefix to avoid collisions under parallel ctest -j.
    const std::string pfx = "/tmp/sr5_" + std::to_string(getpid()) + "_" + tag;
    const std::string inp_path = pfx + ".inp";
    const std::string rpt_path = pfx + ".rpt";
    { std::ofstream f(inp_path); f << inp_text; }

    SWMM_Engine handle = swmm_engine_create();
    if (!handle) return out;
    auto cleanup = [&]() {
        swmm_engine_end(handle);
        swmm_engine_close(handle);
        swmm_engine_destroy(handle);
        std::remove(inp_path.c_str());
        std::remove(rpt_path.c_str());
        std::remove((rpt_path.substr(0, rpt_path.size() - 4) + ".uncertainty.csv").c_str());
    };

    if (swmm_engine_open(handle, inp_path.c_str(), rpt_path.c_str(), nullptr, nullptr) != 0 ||
        swmm_engine_initialize(handle) != 0 ||
        swmm_engine_start(handle, 0) != 0) {
        cleanup();
        return out;
    }

    auto* eng = static_cast<openswmm::SWMMEngine*>(handle);
    // SR-6 falsification knob: scales the K1d handed to the ROM every step
    // (K1d is recomputed per step, so setting it before the first step is
    // sufficient). 1.0 = bit-identical to the gate.
    if (with_soft) eng->rom1dK1dScale() = kK1dScale;
    if (with_soft) eng->rom1dRunoffElasticityEnabled() = (kElasOn != 0);
    const auto& ctx = eng->context();
    const int n_nodes = static_cast<int>(ctx.nodes.head.size());

    std::map<std::string, int> node_idx;
    for (int i = 0; i < n_nodes; ++i) {
        std::string nm = ctx.node_names.name_of(i);
        if (!nm.empty() && nm[0] == 'J') node_idx[nm] = i;
    }

    const auto* rom = with_soft ? eng->rom1d() : nullptr;

    double elapsed = 1.0;
    double t_prev = 0.0;
    int guard = 0, stalled = 0;
    while (elapsed != 0.0 && guard < 20000) {
        if (swmm_engine_step(handle, &elapsed) != 0) break;
        ++guard;
        const double t_now = (elapsed > 0.0) ? elapsed * 86400.0 : kEndTime;
        stalled = (t_now - t_prev < 1e-9) ? stalled + 1 : 0;
        if (stalled > 200) break;

        for (double b = std::floor(t_prev / kReportStep + 1.0) * kReportStep;
             b <= t_now + 1e-6; b += kReportStep) {
            if (b > kEndTime + 1e-6) break;
            out.times.push_back(b);
            for (const auto& [nm, ui] : node_idx) {
                out.heads[nm].push_back(ctx.nodes.head[static_cast<std::size_t>(ui)]);
                out.latq[nm].push_back(ctx.nodes.lat_flow[static_cast<std::size_t>(ui)]);
                if (rom && rom->is_ready() &&
                    ui < static_cast<int>(rom->full_to_active.size())) {
                    const int ai = rom->full_to_active[static_cast<std::size_t>(ui)];
                    if (ai >= 0) {
                        out.q05[nm].push_back(rom->q05[static_cast<std::size_t>(ai)]);
                        out.q50[nm].push_back(rom->q50[static_cast<std::size_t>(ai)]);
                        out.q95[nm].push_back(rom->q95[static_cast<std::size_t>(ai)]);
                    }
                }
            }
        }
        t_prev = t_now;
        if (t_now >= kEndTime) break;
    }
    // Require a non-trivial number of report samples so later indexing is safe.
    out.ok = !out.times.empty();
    if (out.ok) {
        const std::size_t n = out.times.size();
        for (const auto& [nm, _] : out.heads) {
            if ((out.q05.count(nm) && out.q05[nm].size() != n)
                || (out.q50.count(nm) && out.q50[nm].size() != n)
                || (out.q95.count(nm) && out.q95[nm].size() != n))
                out.ok = false;
        }
    }
    cleanup();
    return out;
}

}  // namespace

TEST(SoftRainCoverage, BandsBracketBruteForceMonteCarlo) {
    using openswmm::uncertainty::probit;

    // --- Reference: 21 deterministic runs at NORMAL strata midpoints --------
    // member_i rain scale = 1 + z_i·CV, z_i = probit((i+0.5)/21).
    std::vector<RunResult> mc(kMcRuns);
    for (int i = 0; i < kMcRuns; ++i) {
        const double u = (static_cast<double>(i) + 0.5) / kMcRuns;
        const double z = probit(u);
        const double mult = 1.0 + z * kCV;
        const std::string tag = "mc" + std::to_string(i);
        mc[static_cast<std::size_t>(i)] =
            runCase(fixtureInp(mult, /*with_soft=*/false), tag.c_str(), false);
        ASSERT_TRUE(mc[static_cast<std::size_t>(i)].ok) << "MC run " << i << " failed";
    }

    // --- ROM run -------------------------------------------------------------
    RunResult rom = runCase(fixtureInp(1.0, /*with_soft=*/true), "rom", true);
    ASSERT_TRUE(rom.ok) << "ROM run failed";
    ASSERT_FALSE(rom.q05.empty()) << "ROM produced no soft-rain quantiles";

    // --- Verify all runs share the same number of report samples -------------
    // Use the minimum across all runs as the safe comparison window.
    std::size_t n_samples = rom.times.size();
    for (int i = 0; i < kMcRuns; ++i) {
        n_samples = std::min(n_samples, mc[static_cast<std::size_t>(i)].times.size());
    }
    ASSERT_GT(n_samples, 0u) << "No report samples collected";
    for (const auto& [nm, _] : rom.q05) {
        ASSERT_GE(rom.q05.at(nm).size(), n_samples)
            << "ROM quantile vector too short for node " << nm;
    }

    // --- Compare at every (junction, report time) with t > 60 s -------------
    int n_total = 0, n_covered = 0;
    int n_width = 0, n_width_ok = 0;
    int n_late = 0;
    double member_cov_sum = 0.0, member_cov_late_sum = 0.0;
    double ratio_min = 1e300, ratio_max = 0.0;
    std::vector<double> ratios;

    for (const auto& [nm, rom_q05] : rom.q05) {
        const auto& rom_q95 = rom.q95.at(nm);
        for (std::size_t k = 0; k < n_samples && k < rom_q05.size(); ++k) {
            if (rom.times[k] <= 60.0) continue;
            const bool late = rom.times[k] >= 0.5 * kEndTime;

            std::vector<double> h(kMcRuns);
            for (int i = 0; i < kMcRuns; ++i)
                h[static_cast<std::size_t>(i)] = mc[static_cast<std::size_t>(i)].heads.at(nm)[k];
            std::sort(h.begin(), h.end());
            // Midpoint-plotting-position quantiles (see mc_quantiles.hpp): the
            // ROM's own convention. The former h[19]-h[1] was the 0.071-0.929
            // span and biased the width ratio high by ~12%.
            const double mc_q50   = mcq::quantileMidpoint(h, 0.50);
            const double mc_q05   = mcq::quantileMidpoint(h, 0.05);
            const double mc_q95   = mcq::quantileMidpoint(h, 0.95);
            const double mc_width = mc_q95 - mc_q05;

            ++n_total;
            // Median containment (the pre-review "coverage") and TRUE empirical
            // interval coverage: the fraction of MC members inside the ROM band.
            if (rom_q05[k] <= mc_q50 && mc_q50 <= rom_q95[k]) ++n_covered;
            member_cov_sum += mcq::intervalCoverage(h, rom_q05[k], rom_q95[k]);
            if (late) { member_cov_late_sum += mcq::intervalCoverage(h, rom_q05[k], rom_q95[k]); ++n_late; }

            // Width comparison only in the saturated regime (second half):
            // the deviation-form spread starts at zero and grows toward its
            // parametric steady state, so early-window under-prediction is
            // expected (DEVIATION_FORM.md §4.3).
            if (mc_width > 1e-6 && std::getenv("SR5_TRACE"))
                std::printf("  trace t=%6.0f %s ratio=%.3f rom_w=%.4e mc_w=%.4e\n",
                            rom.times[k], nm.c_str(),
                            (rom_q95[k] - rom_q05[k]) / mc_width,
                            rom_q95[k] - rom_q05[k], mc_width);
            if (late && mc_width > 1e-6) {
                const double rom_width = rom_q95[k] - rom_q05[k];
                const double ratio = rom_width / mc_width;
                ++n_width;
                ratios.push_back(ratio);
                ratio_min = std::min(ratio_min, ratio);
                ratio_max = std::max(ratio_max, ratio);
                if (ratio >= 0.3 && ratio <= 3.0) ++n_width_ok;
            }
        }
    }
    ASSERT_GT(n_total, 0);
    ASSERT_GT(n_width, 0) << "MC produced no resolvable spread — fixture too static";

    const double coverage = static_cast<double>(n_covered) / n_total;   // median containment
    const double member_cov = member_cov_sum / n_total;                  // empirical interval coverage
    const double member_cov_late = (n_late > 0) ? member_cov_late_sum / n_late : 0.0;
    const double width_frac = static_cast<double>(n_width_ok) / n_width;
    std::sort(ratios.begin(), ratios.end());
    const double ratio_med = ratios[ratios.size() / 2];

    {
        // Proportional-scaling check: if head responded exactly in proportion
        // to rain, the MC 5-95 band would be ~3.40*CV*dh (21 midpoint strata)
        // and the ROM band ~3.29*CV*dh (50 strata). Print both ratios at J1 and
        // J5 at the last report time so the shortfall can be attributed to the
        // MC side (elasticity != 1) or the ROM side.
        const std::size_t kk = n_samples - 1;
        for (const char* nm : {"J1", "J5"}) {
            const auto& hs = rom.heads.at(nm);
            const double dh = hs[kk] - hs[0];
            std::vector<double> h(kMcRuns);
            for (int i = 0; i < kMcRuns; ++i) h[static_cast<std::size_t>(i)] = mc[static_cast<std::size_t>(i)].heads.at(nm)[kk];
            std::sort(h.begin(), h.end());
            const double mc_w = mcq::quantileMidpoint(h, 0.95) - mcq::quantileMidpoint(h, 0.05);
            const double rom_w = rom.q95.at(nm)[kk] - rom.q05.at(nm)[kk];
            // Inflow volume into the node over the window vs the stored volume
            // A*dh (A = 20000 ft^2): the ratio is inflow/(inflow - outflow).
            double v_in = 0.0;
            const auto& lq = rom.latq.at(nm);
            for (std::size_t j = 1; j <= kk && j < lq.size(); ++j)
                v_in += 0.5 * (lq[j] + lq[j - 1]) * (rom.times[j] - rom.times[j - 1]);
            // Runoff elasticity to rain, straight from the extreme MC members:
            // ln(V_hi/V_lo) / ln(mult_hi/mult_lo). 1.0 = runoff proportional to rain.
            auto vol = [&](int i) {
                const auto& q = mc[static_cast<std::size_t>(i)].latq.at(nm);
                double v = 0.0;
                for (std::size_t j = 1; j <= kk && j < q.size(); ++j)
                    v += 0.5 * (q[j] + q[j - 1]) * (rom.times[j] - rom.times[j - 1]);
                return v;
            };
            const double z_lo = openswmm::uncertainty::probit(0.5 / kMcRuns);
            const double z_hi = openswmm::uncertainty::probit((kMcRuns - 0.5) / kMcRuns);
            const double elas = std::log(vol(kMcRuns - 1) / vol(0))
                              / std::log((1.0 + z_hi * kCV) / (1.0 + z_lo * kCV));
            const double rain_equiv_cfs = kSubArea * 43560.0 * (kBaseRainInHr / 12.0) / 3600.0;
            std::printf("[SoftRain-elasticity] %s t=%.0f dh_det=%.4f  MC/(3.40*CV*dh)=%.3f  ROM/(3.29*CV*dh)=%.3f  "
                        "V_in/(A*dh)=%.3f  runoff-elasticity=%.3f  mean_runoff/rain=%.2f\n",
                        nm, rom.times[kk], dh, mc_w / (3.40 * kCV * dh), rom_w / (3.29 * kCV * dh),
                        v_in / (20000.0 * dh), elas, (v_in / rom.times[kk]) / rain_equiv_cfs);
        }
        double max_fill = 0.0;
        for (const auto& [nm, hs] : rom.heads) {
            if (nm.size() < 2 || nm[0] != 'J') continue;
            const double invert = 100.0 - (nm[1] - '0');
            for (double h : hs) max_fill = std::max(max_fill, (h - invert) / kPipeDiam);
        }
        std::printf("[SoftRain-fixture] pipe_diam=%.2f ft  max(depth/diam)=%.2f  (>1 = pressurised)\n",
                    kPipeDiam, max_fill);
    }
    std::printf("[SoftRain-vs-MC] samples=%d median-containment=%.3f "
                "member-coverage=%.3f (saturated %.3f, nominal 0.905) width-ratio "
                "min/med/max = %.3f / %.3f / %.3f (in-band frac %.3f of %d)\n",
                n_total, coverage, member_cov, member_cov_late,
                ratio_min, ratio_med, ratio_max, width_frac, n_width);

    // Checklist floors: coverage >= 0.90, width-ratio in [0.3, 3.0] at >= 0.80.
    // Actuals are printed above and recorded in VALIDATION.md on first run.
    //
    // Measured on the first full run (2026-07-16, this fixture):
    //   samples=60  coverage=1.000  width-ratio min/med/max = 0.754/0.821/0.872
    //   in-band fraction 1.000 of 35.
    // The ROM band is slightly narrower than MC (median ratio ~0.82) — the
    // expected mild under-prediction of the delta-linearized soft forcing —
    // but it brackets every MC median and stays well inside [0.3x, 3x].
    // NOTE (2026-10-04 review): `coverage` is MEDIAN CONTAINMENT, not interval
    // coverage. C1 (owner, 2026-10-05): calibrated means member coverage >= 0.80.
    // SR-6 (2026-10-05) raised this cell from 0.761 to 0.881 by propagating the
    // gage spread through the rain->runoff elasticity; it is now a validated
    // cell and asserts the floor (SR5_RUNOFF_ELASTICITY=0 reproduces 0.761).
    EXPECT_TRUE(mcq::reportCalibration("SR-5 soft rain, FULL", member_cov, ratio_med, "SR-6"))
        << "SR-5 is a validated cell since SR-6: member coverage must stay >= 0.80 (measured 0.881)";
    EXPECT_GE(coverage, 0.90)
        << "ROM [q05,q95] must contain the MC median at >=90% of samples";
    EXPECT_GE(width_frac, 0.80)
        << "ROM band width must stay within [0.3x, 3x] of MC width at >=80% "
        << "of saturated-regime samples";
}
