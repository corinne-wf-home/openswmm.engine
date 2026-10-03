/**
 * @file test_rom_coverage_bellinge.cpp
 * @brief PR 10 MC coverage validation — Bellinge model variant.
 *
 * @details The original test_rom_coverage.cpp uses a 5-junction Phase-9 chain
 *          where the dominant ROM mode has τ_0 ≈ 15 hours, so a 1-hour window
 *          only reaches 6.4% saturation (root-cause: network geometry).
 *          Bellinge is a real 1020-node sewer network where K1d ≈ 0.123 is 540×
 *          larger than the Phase-9 fixture, collapsing τ_0 to ~100 seconds. The
 *          1-hour window reaches ~99% saturation — the correct measurement regime.
 *
 *          This test drives the public C API on the real Bellinge SWMM model with
 *          ROM enabled via [UNCERTAINTY] section. It measures the band widths
 *          (coverage, width ratio) when the ROM is provably in saturation, so the
 *          numbers can serve as a regression lock (NOT an MC coverage
 *          validation — it runs the ROM once; see the assertions block).
 *
 *          P6 (2026-10-02): the original loop stopped at 20,000 engine steps,
 *          but 24 h on this model needs 65,589, so the "24h" numbers it
 *          reported (coverage 0.990, max width ratio 0.007) described a
 *          truncated window. The loop now terminates on simulated time and
 *          the test asserts the clock reached 24 h.
 *
 *          Data: set OPENSWMM_BELLINGE_INP / OPENSWMM_BELLINGE_RAIN to point at
 *          the model and rainfall record; the test skips if they are absent.
 *          Runtime 5-10 min (labelled "slow").
 *
 * @ingroup engine_uncertainty
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <map>
#include <string>
#include <vector>

#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "uncertainty/SpectralROM1D.hpp"

namespace fs = std::filesystem;

namespace {

// ============================================================================
// Bellinge model setup — real 1020-node network, run under ROM
// ============================================================================

constexpr double kMcRuns        = 21;     // LHS strata midpoints
constexpr double kPert          = 0.20;   // ±20% Manning's n
constexpr double kRoutingStep   = 30.0;   // s
constexpr double kReportStep    = 300.0; // s (5 min — match longer time window)
constexpr double kEndTime       = 86400.0; // s (24 hours — ensures rainfall activity + saturation)

// Data locations. The model and its rainfall record live OUTSIDE the repo and
// are NOT in the same directory on this machine: the .inp is under
// SWMM_inp/Bellinge/7_SWMM/ but the rainfall .dat is not (P6, 2026-10-02).
// Override with OPENSWMM_BELLINGE_INP / OPENSWMM_BELLINGE_RAIN; the test SKIPS
// (does not fail) when either file is absent.
std::string envOr(const char* name, const std::string& fallback) {
    const char* v = std::getenv(name);
    return (v && *v) ? std::string(v) : fallback;
}

std::string homeDir() {
    const char* h = std::getenv("HOME");
    return h ? std::string(h) : std::string();
}

std::string bellingeInpPath() {
    return envOr("OPENSWMM_BELLINGE_INP",
        homeDir() + "/Projects/SWMM_inp/Bellinge/7_SWMM/BellingeSWMM_v021_nopervious.inp");
}

std::string bellingeRainPath() {
    return envOr("OPENSWMM_BELLINGE_RAIN",
        homeDir() + "/Downloads/article/runs/bellinge/rg_bellinge_Jun2010_Aug2021.dat");
}

// Safety cap on engine steps. It is a SAFETY CAP, never the measurement window:
// the window is kEndTime. Reaching 24 h on this model takes ~65,600 steps
// (measured), so the cap is ~7.6x that. The run result records whether the
// simulated clock actually reached kEndTime, and the test asserts it.
constexpr int kMaxEngineSteps = 500000;

// Modify the .inp to inject [UNCERTAINTY] section, fix reporting, and use absolute rainfall path
// (Bellinge as-shipped has no uncertainty spec and has REPORT_START_TIME=04:00:00)
std::string bellingeWithRomSpec(double rough_mult) {
    // Read the original file
    std::ifstream f(bellingeInpPath());
    std::string content((std::istreambuf_iterator<char>(f)),
                        std::istreambuf_iterator<char>());

    // Replace REPORT_START_TIME to start from beginning
    {
        size_t pos = content.find("REPORT_START_TIME");
        if (pos != std::string::npos) {
            size_t eol = content.find('\n', pos);
            if (eol != std::string::npos) {
                content.replace(pos, eol - pos, "REPORT_START_TIME           00:00:00");
            }
        }
    }

    // Replace relative rainfall path with absolute path
    // Original: "rg_bellinge_Jun2010_Aug2021.dat"
    // Replace with: absolute path to the file
    {
        const std::string abs_path = bellingeRainPath();
        size_t pos = content.find("rg_bellinge_Jun2010_Aug2021.dat");
        while (pos != std::string::npos) {
            // Find the quotes around it
            size_t quote_start = content.rfind('"', pos);
            size_t quote_end = content.find('"', pos + 30);
            if (quote_start != std::string::npos && quote_end != std::string::npos) {
                content.replace(quote_start + 1, quote_end - quote_start - 1, abs_path);
                pos = content.find("rg_bellinge_Jun2010_Aug2021.dat", quote_start + abs_path.length());
            } else {
                break;
            }
        }
    }

    // Inject [UNCERTAINTY] before the first [REPORT] or at EOF
    std::string spec = "[UNCERTAINTY]\n1D MANNINGS_N 0.20\n\n";
    size_t pos = content.find("[REPORT]");
    if (pos != std::string::npos) {
        content.insert(pos, spec);
    } else {
        content += "\n" + spec;
    }

    return content;
}

// ============================================================================
// Engine run helper
// ============================================================================

struct RunResult {
    std::map<std::string, std::vector<double>> heads;  // deterministic
    std::map<std::string, std::vector<double>> q05, q50, q95;  // ROM quantiles
    std::vector<double> times;
    std::map<std::string, double> invert;  // internal feet, junctions only
    bool ok = false;
    int    steps = 0;          // engine steps taken
    double last_t = 0.0;       // simulated seconds reached
    bool   reached_end = false;  // simulated clock reached kEndTime
};

RunResult runCase(const std::string& tag, bool with_rom) {
    RunResult out;
    const std::string pfx = "/tmp/bellinge_rom_cov_" + tag;
    const std::string inp_path = pfx + ".inp";
    const std::string rpt_path = pfx + ".rpt";

    {
        std::ofstream f(inp_path);
        if (with_rom) {
            f << bellingeWithRomSpec(1.0);
        } else {
            std::ifstream src(bellingeInpPath());
            f << std::string((std::istreambuf_iterator<char>(src)),
                            std::istreambuf_iterator<char>());
        }
    }

    SWMM_Engine handle = swmm_engine_create();
    if (!handle) return out;

    auto cleanup = [&]() {
        swmm_engine_end(handle);
        swmm_engine_close(handle);
        swmm_engine_destroy(handle);
        // Keep .inp/.rpt for debugging:
        // std::remove(inp_path.c_str());
        // std::remove(rpt_path.c_str());
        // std::remove((rpt_path.substr(0, rpt_path.size() - 4) + ".uncertainty.csv").c_str());
    };

    if (swmm_engine_open(handle, inp_path.c_str(), rpt_path.c_str(), nullptr, nullptr) != 0 ||
        swmm_engine_initialize(handle) != 0 ||
        swmm_engine_start(handle, 0) != 0) {
        cleanup();
        return out;
    }

    auto* eng = static_cast<openswmm::SWMMEngine*>(handle);
    const auto& ctx = eng->context();
    const int n_nodes = static_cast<int>(ctx.nodes.head.size());

    // Map all junctions (only measure junctions, not outfalls)
    std::map<std::string, int> node_idx;
    for (int i = 0; i < n_nodes; ++i) {
        std::string nm = ctx.node_names.name_of(i);
        if (!nm.empty()) {
            auto ui = static_cast<std::size_t>(i);
            // Include node if it's not an outfall
            if (ctx.nodes.type[ui] != openswmm::NodeType::OUTFALL) {
                node_idx[nm] = i;
                out.invert[nm] = ctx.nodes.invert_elev[ui];
            }
        }
    }

    const auto* rom = with_rom ? eng->rom1d() : nullptr;

    // Step the engine, sampling at report boundaries
    double elapsed = 1.0;
    double t_prev = 0.0;
    int guard = 0;
    std::printf("DEBUG runCase: with_rom=%s, rom=%p\n", with_rom ? "true" : "false", (void*)rom);
    if (rom) {
        std::printf("  rom->is_ready()=%s\n", rom->is_ready() ? "true" : "false");
        std::printf("  rom->n_kept=%d\n", rom->n_kept);
    }

    // Terminate on SIMULATED time (kEndTime), not on a step count. The previous
    // `guard < 20000` stopped at roughly a third of the way through the 24 h
    // window (65,589 steps needed), before or partway through the storm.
    while (elapsed != 0.0 && guard < kMaxEngineSteps) {
        if (swmm_engine_step(handle, &elapsed) != 0) break;
        ++guard;
        const double t_now = (elapsed > 0.0) ? elapsed * 86400.0 : kEndTime;

        // Report boundary crossing
        for (double b = std::floor(t_prev / kReportStep + 1.0) * kReportStep;
             b <= t_now; b += kReportStep) {
            out.times.push_back(b);

            for (const auto& [nm, idx] : node_idx) {
                auto ui = static_cast<std::size_t>(idx);
                out.heads[nm].push_back(ctx.nodes.head[ui]);
                if (rom && rom->is_ready()) {
                    out.q05[nm].push_back(rom->q05[ui]);
                    out.q50[nm].push_back(rom->q50[ui]);
                    out.q95[nm].push_back(rom->q95[ui]);
                }
            }
        }

        t_prev = t_now;
        out.last_t = t_now;
        if (t_now >= kEndTime) { out.reached_end = true; break; }
    }
    // The engine can also finish on its own (elapsed == 0) exactly at END_TIME.
    if (elapsed == 0.0 && t_prev >= kEndTime - kReportStep) out.reached_end = true;
    out.steps = guard;

    out.ok = !out.times.empty();
    cleanup();
    return out;
}

}  // namespace

// ============================================================================
// Test: Bellinge ROM coverage at saturation
// ============================================================================

TEST(RomCoverageBellinge, BandsBracketBruteForceMonteCarlo) {
    // This test is PLACEHOLDER/DEMONSTRATOR only.
    // Once Bellinge-based ROM coverage numbers are measured and confirmed,
    // this will be the baseline for H5/H11 work (replacing or supplementing
    // the original Phase-9 test).

    if (!fs::exists(bellingeInpPath()) || !fs::exists(bellingeRainPath())) {
        GTEST_SKIP() << "Bellinge data not found. inp=" << bellingeInpPath()
                     << " rain=" << bellingeRainPath()
                     << " (override with OPENSWMM_BELLINGE_INP / OPENSWMM_BELLINGE_RAIN)";
    }

    // Run the ROM once to measure baseline band widths
    const auto rom = runCase("bellinge_rom", /*with_rom=*/true);
    ASSERT_TRUE(rom.ok) << "Bellinge ROM run failed";

    // THE P6 GUARD: the measurement window is only valid if the simulated
    // clock really reached the end. A silently truncated window is exactly how
    // the old baseline went wrong.
    std::printf("RUN: %d engine steps, simulated clock reached %.2f h (target %.2f h)\n",
                rom.steps, rom.last_t / 3600.0, kEndTime / 3600.0);
    ASSERT_TRUE(rom.reached_end)
        << "Run stopped at " << rom.last_t / 3600.0 << " h after " << rom.steps
        << " steps — the 24 h window was NOT completed (step cap "
        << kMaxEngineSteps << "). Numbers below would describe a truncated window.";

    // Debug output
    std::printf("DEBUG: rom.q05.size() = %zu, rom.q50.size() = %zu, rom.q95.size() = %zu\n",
                rom.q05.size(), rom.q50.size(), rom.q95.size());
    std::printf("DEBUG: rom.times.size() = %zu, rom.heads.size() = %zu\n",
                rom.times.size(), rom.heads.size());

    ASSERT_FALSE(rom.q05.empty()) << "ROM produced no q05 quantiles";

    // Placeholder assertion: verify ROM ran and produced spreads
    // (Real assertions would compare against MC runs, once they're generated)
    int total_q95_q05 = 0, nonzero_spreads = 0;
    int sample_count = 0;
    for (const auto& [nm, q05_vec] : rom.q05) {
        const auto it_q95 = rom.q95.find(nm);
        const auto it_q50 = rom.q50.find(nm);
        if (it_q95 == rom.q95.end()) continue;
        const auto& q95_vec = it_q95->second;
        const auto& q50_vec = (it_q50 != rom.q50.end()) ? it_q50->second : q95_vec;

        for (std::size_t i = 0; i < q95_vec.size() && i < q05_vec.size(); ++i) {
            const double spread = q95_vec[i] - q05_vec[i];
            ++total_q95_q05;
            if (spread > 1e-8) {
                ++nonzero_spreads;
            }
            if (sample_count < 5) {
                std::printf("  %s[%zu]: q05=%e, q50=%e, q95=%e, spread=%e\n",
                            nm.c_str(), i, q05_vec[i], q50_vec[i], q95_vec[i], spread);
                ++sample_count;
            }
        }
    }

    std::printf("SUMMARY: %d spreads > 1e-8 out of %d (%.1f%%) measurements\n",
                nonzero_spreads, total_q95_q05, 100.0 * nonzero_spreads / std::max(1, total_q95_q05));

    // Compute statistics for the regression baseline
    double min_ratio = 1e100, med_ratio = 0, max_ratio = 0;
    std::vector<double> all_ratios;
    int in_band_count = 0; // spreads in [0.3*median, 3*median] range (will refine after median)

    for (const auto& [nm, q05_vec] : rom.q05) {
        const auto it_q95 = rom.q95.find(nm);
        const auto it_q50 = rom.q50.find(nm);
        if (it_q95 == rom.q95.end()) continue;
        const auto& q95_vec = it_q95->second;
        const auto& q50_vec = (it_q50 != rom.q50.end()) ? it_q50->second : q95_vec;

        for (std::size_t i = 0; i < q95_vec.size() && i < q05_vec.size(); ++i) {
            const double spread = q95_vec[i] - q05_vec[i];
            if (spread < 1e-8) continue; // ignore machine epsilon
            const double median = q50_vec[i];
            if (median < 1e-6) continue; // avoid division by zero / numerical noise
            const double ratio = spread / median;
            all_ratios.push_back(ratio);
            min_ratio = std::min(min_ratio, ratio);
            max_ratio = std::max(max_ratio, ratio);
        }
    }

    // Compute median ratio
    if (!all_ratios.empty()) {
        std::sort(all_ratios.begin(), all_ratios.end());
        med_ratio = all_ratios[all_ratios.size() / 2];

        // Count in-band: [0.3*median, 3*median]
        const double lower = 0.3 * med_ratio;
        const double upper = 3.0 * med_ratio;
        for (double r : all_ratios) {
            if (r >= lower && r <= upper) ++in_band_count;
        }
    }

    // Band relative to LOCAL WATER DEPTH (q50 - invert). The ratio above divides
    // by absolute head (~328 ft on Bellinge, internal feet), which mostly
    // measures how high the network sits above datum, not how uncertain the
    // water level is. Depth is the quantity comparable between a metres-deep
    // trunk and centimetre-deep laterals. Nodes shallower than kMinDepthFt are
    // excluded (depth ~ 0 makes the ratio meaningless).
    constexpr double kMinDepthFt = 0.05;
    std::vector<double> depth_ratios, abs_bands_ft;
    for (const auto& [nm, q05_vec] : rom.q05) {
        const auto it_q95 = rom.q95.find(nm);
        const auto it_q50 = rom.q50.find(nm);
        const auto it_inv = rom.invert.find(nm);
        if (it_q95 == rom.q95.end() || it_q50 == rom.q50.end() ||
            it_inv == rom.invert.end()) continue;
        const std::size_t n = std::min({q05_vec.size(), it_q95->second.size(),
                                        it_q50->second.size()});
        for (std::size_t i = 0; i < n; ++i) {
            const double band  = it_q95->second[i] - q05_vec[i];
            const double depth = it_q50->second[i] - it_inv->second;
            abs_bands_ft.push_back(band);
            if (depth >= kMinDepthFt) depth_ratios.push_back(band / depth);
        }
    }
    auto pct = [](std::vector<double> v, double q) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[static_cast<std::size_t>(q * (v.size() - 1))];
    };
    std::printf("Band / absolute head      : p50=%.4f p90=%.4f p99=%.4f max=%.4f (n=%zu)\n",
                pct(all_ratios, 0.5), pct(all_ratios, 0.9), pct(all_ratios, 0.99),
                max_ratio, all_ratios.size());
    std::printf("Band / local depth (>=%.2f ft): p50=%.4f p90=%.4f p99=%.4f max=%.4f (n=%zu)\n",
                kMinDepthFt, pct(depth_ratios, 0.5), pct(depth_ratios, 0.9),
                pct(depth_ratios, 0.99), pct(depth_ratios, 1.0), depth_ratios.size());
    std::printf("Absolute band (ft)        : p50=%.4f p90=%.4f p99=%.4f max=%.4f\n",
                pct(abs_bands_ft, 0.5), pct(abs_bands_ft, 0.9), pct(abs_bands_ft, 0.99),
                pct(abs_bands_ft, 1.0));

    const double in_band_pct = (all_ratios.empty()) ? 0.0 :
        (100.0 * in_band_count / static_cast<int>(all_ratios.size()));

    std::printf("=== BELLINGE ROM BAND BASELINE (full 24 h window) ===\n");
    std::printf("Coverage: %.3f (spreads on %.1f%% of measurements)\n",
                static_cast<double>(nonzero_spreads) / total_q95_q05,
                100.0 * nonzero_spreads / total_q95_q05);
    std::printf("Width ratio:\n");
    std::printf("  min  = %.3f\n", min_ratio);
    std::printf("  med  = %.3f\n", med_ratio);
    std::printf("  max  = %.3f\n", max_ratio);
    std::printf("In-band [0.3*med, 3*med]: %.1f%%\n", in_band_pct);

    // ---- Assertions: a REGRESSION LOCK on the full 24 h window ----------------
    // This test runs the ROM ONCE; there is no brute-force MC here, so these
    // are change detectors around measured values, NOT a coverage validation
    // ("coverage" above is the fraction of samples with nonzero spread). Real
    // MC coverage on Bellinge at storm peak does not exist; the band magnitude
    // in the surcharging regime is UNVALIDATED (see VALIDATION.md, Bellinge).
    //
    // Measured on the COMPLETE window (65,589 steps, 24.00 h, 2026-10-02):
    //   nonzero-spread fraction 0.988
    //   band / absolute head : p50 0.0002  p90 0.0013  p99 0.0118  max 0.277
    //   band / local depth   : p50 0.153   p90 0.73    p99 2.00    max 10.7
    // The previous baseline (max 0.007, coverage 0.990) came from a run that
    // stopped at 20,000 steps and never reached the storm; it was ~40x too low.
    EXPECT_GT(nonzero_spreads, 0)
        << "ROM produced zero spread everywhere — not exercising uncertainty";
    EXPECT_GE(static_cast<double>(nonzero_spreads) / total_q95_q05, 0.90)
        << "Should have spreads on >=90% of measurements (measured 0.988)";
    EXPECT_LE(max_ratio, 0.5)
        << "Max band/absolute-head should stay bounded (measured 0.277)";
    // The band relative to local depth is the comparable quantity; pin its
    // median to a wide bracket around the measured 0.153 so a collapse to ~0
    // (band vanishes) or a blow-up (>0.5 median) is caught.
    const double depth_p50 = pct(depth_ratios, 0.5);
    EXPECT_GT(depth_p50, 0.02) << "Median band/depth collapsed (measured 0.153)";
    EXPECT_LT(depth_p50, 0.5)  << "Median band/depth blew up (measured 0.153)";
    // No assertion on "in-band": ratio-vs-its-own-median is not a meaningful band.

    std::printf("Bellinge ROM full-window regression lock passed (change detector, not MC validation).\n");
}
