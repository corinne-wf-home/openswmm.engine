/**
 * @file test_engine_wq_uncertainty.cpp
 * @brief Integration test: Water-quality uncertainty layer engine lifecycle.
 *
 * @details Tests the full pipeline:
 *   parse .inp with [UNCERTAINTY] QUALITY TSS 0.30
 *   → swmm_engine_open → initialize → start → step loop → end → close
 *
 * Verifications:
 *   1. wq_unc_active() returns true after initialize().
 *   2. wq_conc_prev_ buffer is properly sized.
 *   3. WQ uncertainty CSV file is created.
 *   4. Pollutant name resolution works correctly.
 *   5. WQ bounds are computed and written to CSV.
 *
 * Network: single node with TSS pollutant.
 * Forced inflow at node to ensure non-trivial concentration field.
 *
 * @ingroup engine_integration
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "openswmm/engine/openswmm_engine.h"
#define private public
#include "core/SWMMEngine.hpp"
#undef private

namespace {

// Path to the test input file
static const char* k_inp_path = "wq_uncertainty_test.inp";
static const char* k_rpt_path = "/tmp/wq_uncertainty_test.rpt";

} // anonymous namespace


// ============================================================================
// Test: WQ uncertainty layer integration
// ============================================================================

TEST(EngineWQUncertainty, LifecycleIntegration) {
    // Create engine
    SWMM_Engine handle = swmm_engine_create();
    ASSERT_NE(handle, nullptr) << "Failed to create engine";

    // Open engine
    int error = swmm_engine_open(handle, k_inp_path, k_rpt_path, nullptr, nullptr);
    ASSERT_EQ(error, 0) << "Failed to open engine";

    // Get engine reference
    auto* engine = static_cast<openswmm::SWMMEngine*>(handle);

    // Initialize
    error = swmm_engine_initialize(handle);
    ASSERT_EQ(error, 0) << "Failed to initialize engine: " << error;

    // Verify WQ uncertainty is active
    EXPECT_TRUE(engine->wq_unc_active()) << "WQ uncertainty should be active";

    // Verify the report-boundary seed matches the deterministic start state.
    EXPECT_EQ(engine->wq_conc_prev_.size(), engine->ctx_.nodes.conc.size());
    EXPECT_THAT(engine->wq_conc_prev_, ::testing::ElementsAreArray(engine->ctx_.nodes.conc));

    // Verify pollutant indices are resolved
    // Note: This would require access to private members or additional API

    // Start simulation
    error = swmm_engine_start(handle, 0);
    ASSERT_EQ(error, 0) << "Failed to start engine: " << error;

    // Run the simulation to completion so report boundaries are reached.
    double elapsed_time = 1.0;
    do {
        error = swmm_engine_step(handle, &elapsed_time);
        ASSERT_EQ(error, 0) << "Failed to step engine: " << error;
    } while (elapsed_time > 0.0);

    // End simulation
    error = swmm_engine_end(handle);
    ASSERT_EQ(error, 0) << "Failed to end engine: " << error;

    // Close engine
    error = swmm_engine_close(handle);
    ASSERT_EQ(error, 0) << "Failed to close engine: " << error;
}

// ============================================================================
// Test: WQ uncertainty CSV output
// ============================================================================

TEST(EngineWQUncertainty, CSVOutput) {
    // Create engine
    SWMM_Engine handle = swmm_engine_create();
    ASSERT_NE(handle, nullptr) << "Failed to create engine";

    // Open engine
    int error = swmm_engine_open(handle, k_inp_path, k_rpt_path, nullptr, nullptr);
    ASSERT_EQ(error, 0) << "Failed to open engine";

    // Initialize
    error = swmm_engine_initialize(handle);
    ASSERT_EQ(error, 0) << "Failed to initialize engine: " << error;

    // Start simulation
    error = swmm_engine_start(handle, 0);
    ASSERT_EQ(error, 0) << "Failed to start engine: " << error;

    // Run the simulation to completion so report boundaries are reached.
    double elapsed_time = 1.0;
    do {
        error = swmm_engine_step(handle, &elapsed_time);
        ASSERT_EQ(error, 0) << "Failed to step engine: " << error;
    } while (elapsed_time > 0.0);

    // End simulation
    error = swmm_engine_end(handle);
    ASSERT_EQ(error, 0) << "Failed to end engine: " << error;

    // Close engine
    error = swmm_engine_close(handle);
    ASSERT_EQ(error, 0) << "Failed to close engine: " << error;

    // Check that WQ uncertainty CSV was created
    std::string csv_path = k_rpt_path;
    const std::string rpt_ext = ".rpt";
    if (csv_path.size() >= rpt_ext.size() &&
        csv_path.compare(csv_path.size() - rpt_ext.size(), rpt_ext.size(), rpt_ext) == 0) {
        csv_path.replace(csv_path.size() - rpt_ext.size(), rpt_ext.size(), ".wq_uncertainty.csv");
    } else {
        csv_path += ".wq_uncertainty.csv";
    }
    std::ifstream csv_file(csv_path);
    EXPECT_TRUE(csv_file.good()) << "WQ uncertainty CSV file should be created";

    std::string header;
    ASSERT_TRUE(std::getline(csv_file, header));
    EXPECT_EQ(header, "time_s,node_name,pollutant,q05,q50,q95");

    std::string row;
    ASSERT_TRUE(std::getline(csv_file, row))
        << "WQ uncertainty CSV must contain at least one data row";
    std::stringstream row_stream(row);
    std::vector<std::string> cols;
    std::string col;
    while (std::getline(row_stream, col, ',')) cols.push_back(col);
    ASSERT_EQ(cols.size(), 6u);
    EXPECT_EQ(cols[2], "TSS");

    // Clean up
    std::remove(k_rpt_path);
    std::remove(csv_path.c_str());
}
// ============================================================================
// H6: content, not shape. The lifecycle fixture above has no TSS source, so
// every concentration (and every band) is zero -- it proves the file exists,
// not that the numbers are right. These cases put a real TSS load on J1, run
// each distribution family through the whole engine, and check the CSV rows
// against the closed forms the unit suite pins.
// ============================================================================

#include <cmath>
#include <filesystem>
#include <map>

namespace {

struct WqRow { double t, q05, q50, q95; };

// Fixture: J1 receives 0.1 CMS dry-weather flow carrying 100 mg/L TSS with a
// first-order decay of k_day; 10-minute run, 1-minute report step.
//
// The TSS load rides the pollutant-wide DWF concentration (the Cdwf column of
// [POLLUTANTS]), NOT a per-node "[DWF] J1 TSS 100" row: on this line's
// engine (merge-base with swmm6_rel 2026-07-06) a pollutant [DWF] row is
// summed into the node's dry-weather FLOW (100 mg/L booked as 100 cfs) and
// its load is never applied. Upstream swmm6_rel fixed that (FLOW rows only
// in the water pass); the defect is out of this test's scope and flagged in
// VALIDATION.md's H6 section.
std::string wqContentInp(const std::string& qual_line, double k_day) {
    std::ostringstream o;
    o << "[TITLE]\nH6 WQ content fixture\n\n"
         "[OPTIONS]\nFLOW_UNITS CMS\nFLOW_ROUTING DYNWAVE\n"
         "START_DATE 01/01/2025\nSTART_TIME 00:00:00\n"
         "REPORT_START_DATE 01/01/2025\nREPORT_START_TIME 00:00:00\n"
         "END_DATE 01/01/2025\nEND_TIME 00:10:00\nREPORT_STEP 00:01:00\n"
         "ROUTING_STEP 0:00:15\nMIN_SURFAREA 1.0\nTHREADS 1\n\n"
         "[POLLUTANTS]\nTSS mg/L 0.0 0.0 0.0 " << k_day << " NO * 0.0 100 0.0\n\n"
         "[JUNCTIONS]\nJ1 100 5 1.0 0 0\n\n"
         "[OUTFALLS]\nO1 95 FREE NO\n\n"
         "[CONDUITS]\nC1 J1 O1 100 0.013 0 0 0 0\n\n"
         "[XSECTIONS]\nC1 CIRCULAR 1.0 0 0 0 1\n\n"
         "[DWF]\nJ1 FLOW 0.1\n\n"
         "[UNCERTAINTY]\n" << qual_line << "\n";
    return o.str();
}

// Runs the fixture and returns J1's rows that carry a non-zero concentration,
// plus the number of warnings the engine emitted.
std::vector<WqRow> runWqContent(const std::string& qual_line, double k_day,
                                int* n_warnings = nullptr,
                                bool* csv_has_rows = nullptr) {
    namespace fs = std::filesystem;
    static int serial = 0;
    const fs::path dir = fs::temp_directory_path();
    const std::string stem = "h6_wq_content_" + std::to_string(serial++);
    const std::string inp = (dir / (stem + ".inp")).string();
    const std::string rpt = (dir / (stem + ".rpt")).string();
    const std::string csv = (dir / (stem + ".wq_uncertainty.csv")).string();
    { std::ofstream f(inp); f << wqContentInp(qual_line, k_day); }

    std::vector<WqRow> rows;
    SWMM_Engine h = swmm_engine_create();
    EXPECT_EQ(swmm_engine_open(h, inp.c_str(), rpt.c_str(), nullptr, nullptr), 0);
    EXPECT_EQ(swmm_engine_initialize(h), 0);
    EXPECT_EQ(swmm_engine_start(h, 0), 0);
    double elapsed = 1.0;
    int guard = 0;
    do { if (swmm_engine_step(h, &elapsed) != 0) break; }
    while (elapsed > 0.0 && ++guard < 100000);
    if (n_warnings) *n_warnings = swmm_get_warning_count(h);
    swmm_engine_end(h);
    swmm_engine_close(h);
    swmm_engine_destroy(h);

    std::ifstream f(csv);
    std::string line;
    std::getline(f, line);                       // header
    bool any = false;
    while (std::getline(f, line)) {
        any = true;
        std::stringstream ss(line);
        std::string t, node, pol, a, b, c;
        std::getline(ss, t, ','); std::getline(ss, node, ',');
        std::getline(ss, pol, ','); std::getline(ss, a, ',');
        std::getline(ss, b, ','); std::getline(ss, c, ',');
        if (node != "J1" || pol != "TSS") continue;
        const WqRow r{std::stod(t), std::stod(a), std::stod(b), std::stod(c)};
        if (r.q50 > 1.0) rows.push_back(r);
    }
    if (csv_has_rows) *csv_has_rows = any;
    std::remove(inp.c_str()); std::remove(rpt.c_str()); std::remove(csv.c_str());
    return rows;
}

// k = 0.02/day over a 60 s report step: k*dt = 1.39e-5, so the UNIFORM band's
// relative width is k*dt*(0.9*2p) = 7.5e-6 at p = 0.30. Deliberately small:
// this line's routing applies Kdecay per SECOND (stored per day, never
// converted -- upstream KD1 fixed it), so any k near 1/day zeroes J1's
// concentration in routing before the band layer ever sees it. The band
// layer reads the stored 1/day value, so its width is the physical one.
constexpr double kDay = 0.02;
constexpr double kUniformRelWidth = 0.02 * (60.0 / 86400.0) * 0.9 * 2.0 * 0.30;

} // namespace

TEST(EngineWQUncertaintyContent, RealLoadGivesOrderedNonZeroBands) {
    const auto rows = runWqContent("QUALITY TSS 0.30", kDay);
    ASSERT_GE(rows.size(), 5u) << "J1 should carry TSS at most report boundaries";
    for (const auto& r : rows) {
        EXPECT_LT(r.q05, r.q50) << "t=" << r.t;
        EXPECT_LT(r.q50, r.q95) << "t=" << r.t;
        // The physical width: k [1/day] x report step [days] x 1.8p. Off by
        // 86,400x if the per-day/per-second convention drifts on a rebase.
        EXPECT_NEAR((r.q95 - r.q05) / r.q50, kUniformRelWidth,
                    0.05 * kUniformRelWidth) << "t=" << r.t;
    }
}

// Same run, three families: NORMAL ≈ 0.61× UNIFORM's width, LOGNORMAL skewed
// toward more decay -- the unit suite's closed forms, now through the parser,
// the engine's late name resolution and the report-boundary writer.
TEST(EngineWQUncertaintyContent, DistributionFamilyReachesTheBand) {
    const auto u = runWqContent("QUALITY TSS 0.30 UNIFORM",   kDay);
    const auto n = runWqContent("QUALITY TSS 0.30 NORMAL",    kDay);
    const auto l = runWqContent("QUALITY TSS 0.30 LOGNORMAL", kDay);
    ASSERT_EQ(u.size(), n.size());
    ASSERT_EQ(u.size(), l.size());
    ASSERT_GE(u.size(), 5u);
    for (std::size_t i = 0; i < u.size(); ++i) {
        // Deterministic routing is identical across the three runs.
        EXPECT_DOUBLE_EQ(u[i].q50, n[i].q50) << "row " << i;
        const double wu = u[i].q95 - u[i].q05, wn = n[i].q95 - n[i].q05;
        EXPECT_NEAR(wn / wu, 0.609, 0.04) << "row " << i;
        const double lo = l[i].q50 - l[i].q05, hi = l[i].q95 - l[i].q50;
        EXPECT_NEAR(lo / hi, 1.30, 0.10) << "row " << i;
    }
}

// No decay -> no decay-rate uncertainty: q05 == q50 == q95, exactly.
TEST(EngineWQUncertaintyContent, ZeroDecayGivesZeroWidth) {
    const auto rows = runWqContent("QUALITY TSS 0.30 LOGNORMAL", 0.0);
    ASSERT_GE(rows.size(), 5u);
    for (const auto& r : rows) {
        EXPECT_EQ(r.q05, r.q50);
        EXPECT_EQ(r.q50, r.q95);
    }
}

// An unknown pollutant name is warned about at init and reported nowhere;
// the run itself is unaffected.
TEST(EngineWQUncertaintyContent, UnknownPollutantWarnsAndWritesNoRows) {
    int warnings = 0;
    bool any_rows = true;
    const auto rows = runWqContent("QUALITY NOT_A_POLLUTANT 0.30", kDay,
                                   &warnings, &any_rows);
    EXPECT_TRUE(rows.empty());
    EXPECT_FALSE(any_rows);
    EXPECT_GE(warnings, 1);
}
