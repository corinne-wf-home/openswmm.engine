/**
 * @file test_rom_coverage_structures.cpp
 * @brief H17 — 1D ROM bands vs brute-force Monte Carlo on a network with
 *        hydraulic structures (EPA Example 8: one pump, one orifice, four
 *        weirs, a storage well and a fixed-stage outfall).
 *
 * Reference: 21 deterministic runs of the unmodified Example 8 deck with
 * every CONDUIT's Manning n scaled by the LHS strata midpoints of +/-20%
 * (structures are not conduits and are left alone -- the same parameter
 * the ROM's `1D MANNINGS_N` perturbs). ROM run: same deck plus
 * `[UNCERTAINTY] 1D MANNINGS_N 0.20`.
 *
 * Every junction and the storage node are scored per (node, report time)
 * with the MC-side quantiles at midpoint plotting positions (mc_quantiles.hpp).
 * The per-node table is printed so the effect of the structure edges can
 * be read node by node (the interceptor JI1..JI5 behind the orifice, the
 * force-main JI10..JI13 behind the pump, and the weir-fed surface nodes).
 *
 * MEASURED 2026-10-12 (member coverage / width median, all samples):
 *   conduit-only ROM graph, depth x elasticity source   0.761 / 1.636
 *   + structure edges                                    0.773 / 1.966
 *   + flow-response source                               0.757 / 1.046
 *   + quasi-static residual (shipped)                    0.810 / 1.102  CALIBRATED
 * The widths are right node by node (JI4 15.9 -> 1.03, JI5 13.6 -> 0.63,
 * Well 4.6 -> 0.99, force main while pumping 7 ft -> 0.26 ft vs MC 0.21 ft).
 * The quasi-static residual (SpectralROM1D::quasi_static_residual) lifts the
 * cell over the floor: the kept modes no longer smear a localized response.
 * The remaining misses are centimetre-scale spreads at nodes with no steady
 * flow: the
 * force main after the pump stops (JI10..JI13) and the weir-fed surface
 * channel (J5, J8, J9, J11). There the members differ by stored VOLUME --
 * when the pump stopped, how much spilled -- which a steady head response
 * cannot carry. Recorded in VALIDATION.md "H17".
 *
 * Asserted: the run completes, the width median stays inside the standing
 * [0.3, 3.0] sanity band (measured 1.046), and the samples the output CSV
 * flags `steady` are calibrated (step 4; measured 0.944 / 1.350). The pooled
 * C1 verdict (0.810, margin 0.01) is printed, not asserted; the misses sit
 * under the `no_flow` (0.694 / 0.213) and `transient` (0.465 / 0.393) flags,
 * which is what those flags tell the user.
 *
 * Env knobs: H17_WINDOW_H (sampling window, default the full 12 h),
 * H17_NODE=a,b prints a time trace for those nodes, H17_FLOW_SOURCE=0
 * restores the pre-H17 source for A/B.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "uncertainty/SpectralROM1D.hpp"
#include "mc_quantiles.hpp"

#if defined(_WIN32)
#  include <process.h>
#  define getpid _getpid
#else
#  include <unistd.h>
#endif

#ifndef REGRESSION_DATA_DIR
#  define REGRESSION_DATA_DIR "."
#endif

namespace {

constexpr double kPert       = 0.20;
constexpr int    kMcRuns     = 21;
constexpr double kReportStep = 60.0;

std::string readFile(const std::string& path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// Scale column 5 (Manning N) of every [CONDUITS] data line by `mult`, and
// optionally append an [UNCERTAINTY] block.
std::string deck(const std::string& base, double mult, bool with_rom) {
    std::istringstream in(base);
    std::ostringstream out;
    std::string line, section;
    while (std::getline(in, line)) {
        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        if (!trimmed.empty() && trimmed[0] == '[') {
            section = trimmed.substr(0, trimmed.find(']') + 1);
            out << line << '\n';
            continue;
        }
        if (section == "[CONDUITS]" && !trimmed.empty() && trimmed[0] != ';') {
            std::istringstream ls(trimmed);
            std::vector<std::string> tok;
            std::string t;
            while (ls >> t) tok.push_back(t);
            if (tok.size() >= 5) {
                char buf[64];
                std::snprintf(buf, sizeof buf, "%.8g", std::atof(tok[4].c_str()) * mult);
                tok[4] = buf;
                std::string joined;
                for (const auto& s : tok) joined += s + "  ";
                out << joined << '\n';
                continue;
            }
        }
        out << line << '\n';
    }
    if (with_rom) out << "\n[UNCERTAINTY]\n1D MANNINGS_N 0.20\n";
    return out.str();
}

std::set<std::string> scoredNodes(const std::string& base) {
    std::set<std::string> names;
    std::istringstream in(base);
    std::string line, section;
    while (std::getline(in, line)) {
        std::string tr = line;
        tr.erase(0, tr.find_first_not_of(" \t"));
        if (!tr.empty() && tr[0] == '[') { section = tr.substr(0, tr.find(']') + 1); continue; }
        if ((section == "[JUNCTIONS]" || section == "[STORAGE]") && !tr.empty() && tr[0] != ';') {
            std::istringstream ls(tr);
            std::string nm;
            ls >> nm;
            names.insert(nm);
        }
    }
    return names;
}

struct RunResult {
    std::map<std::string, std::vector<double>> heads, q05, q95;
    std::map<std::string, std::vector<int>> flag;   // H17 step 4 (ROM run)
    std::vector<double> times;
    bool ok = false;
};

RunResult runCase(const std::string& inp_text, const std::string& tag, bool with_rom,
                  const std::set<std::string>& names, double end_time_s) {
    RunResult out;
    const std::string pfx = "/tmp/rom_struct_" + std::to_string(getpid()) + "_" + tag;
    const std::string inp_path = pfx + ".inp", rpt_path = pfx + ".rpt";
    { std::ofstream f(inp_path); f << inp_text; }

    SWMM_Engine h = swmm_engine_create();
    auto cleanup = [&]() {
        swmm_engine_end(h);
        swmm_engine_close(h);
        swmm_engine_destroy(h);
        std::remove(inp_path.c_str());
        std::remove(rpt_path.c_str());
        std::remove((pfx + ".uncertainty.csv").c_str());
        std::remove((pfx + ".rom_diag.csv").c_str());
        std::remove((pfx + ".rom_threshold.csv").c_str());
    };
    bool opened = swmm_engine_open(h, inp_path.c_str(), rpt_path.c_str(), nullptr, nullptr) == 0;
    if (opened) {
        auto* e0 = static_cast<openswmm::SWMMEngine*>(h);
        if (const char* v = std::getenv("H17_K_MODES")) e0->rom1dModeCountOverride() = std::atoi(v);
        if (const char* v = std::getenv("H17_K1D_SCALE")) e0->rom1dK1dScale() = std::atof(v);
        if (const char* v = std::getenv("H17_ADV")) e0->rom1dK1dAdvective() = std::atoi(v) != 0;
        if (const char* q = std::getenv("H17_QSR")) e0->rom1dQuasiStatic() = std::atoi(q) != 0;
    }
    if (!opened || swmm_engine_initialize(h) != 0 || swmm_engine_start(h, 0) != 0) {
        std::printf("[H17] %s failed to start: %s\n", tag.c_str(), swmm_get_last_error_msg(h));
        cleanup();
        return out;
    }
    auto* eng = static_cast<openswmm::SWMMEngine*>(h);
    const auto& ctx = eng->context();
    std::map<std::string, int> idx;
    for (int i = 0; i < static_cast<int>(ctx.nodes.head.size()); ++i) {
        const std::string nm = ctx.node_names.name_of(i);
        if (names.count(nm)) idx[nm] = i;
    }
    const auto* rom = with_rom ? eng->rom1d() : nullptr;
    if (with_rom && !rom) { cleanup(); return out; }
    if (const char* v = std::getenv("H17_FLOW_SOURCE")) eng->rom1dFlowSource() = std::atoi(v) != 0;

    double elapsed = 1.0, t_prev = 0.0;
    int stalled = 0;
    std::ofstream phi_out;
    std::map<std::string, int> phi_nodes;
    if (with_rom && std::getenv("H17_PHI_OUT") && std::getenv("H17_PHI_NODES")) {
        phi_out.open(std::getenv("H17_PHI_OUT"));
        std::stringstream ss(std::getenv("H17_PHI_NODES"));
        std::string nm;
        while (std::getline(ss, nm, ',')) if (idx.count(nm)) phi_nodes[nm] = idx.at(nm);
    }
    // H2 re-measure (2026-10-12): H17_PUMP_TRACE=1 prints every pump on/off
    // switch and every basis rebuild (cold-forced or not) with its time.
    const bool pump_trace = with_rom && std::getenv("H17_PUMP_TRACE");
    std::vector<int> pumps;
    std::vector<char> pump_on;
    int last_att = 0, last_cold = 0;
    if (pump_trace) {
        for (int j = 0; j < static_cast<int>(ctx.links.type.size()); ++j)
            if (ctx.links.type[static_cast<std::size_t>(j)] == openswmm::LinkType::PUMP) {
                pumps.push_back(j); pump_on.push_back(0);
            }
    }
    for (int guard = 0; elapsed != 0.0 && guard < 500000; ++guard) {
        if (swmm_engine_step(h, &elapsed) != 0) break;
        const double t_now = (elapsed > 0.0) ? elapsed * 86400.0 : end_time_s;
        if (pump_trace && rom) {
            for (std::size_t p = 0; p < pumps.size(); ++p) {
                const char on = std::fabs(ctx.links.flow[static_cast<std::size_t>(pumps[p])]) > 1e-6;
                if (on != pump_on[p]) {
                    std::printf("[H2-trace] t=%7.0f pump %s %s\n", t_now,
                                ctx.link_names.name_of(pumps[p]).c_str(), on ? "ON" : "OFF");
                    pump_on[p] = on;
                }
            }
            if (rom->basis_updates_attempted_ != last_att || rom->basis_rebuilds_cold_forced_ != last_cold) {
                std::printf("[H2-trace] t=%7.0f rebuild%s (attempts %d, cold-forced %d)\n", t_now,
                            rom->basis_rebuilds_cold_forced_ != last_cold ? " COLD" : "",
                            rom->basis_updates_attempted_, rom->basis_rebuilds_cold_forced_);
                last_att = rom->basis_updates_attempted_; last_cold = rom->basis_rebuilds_cold_forced_;
            }
        }
        stalled = (t_now - t_prev < 1e-9) ? stalled + 1 : 0;
        if (stalled > 200) break;
        for (double b = std::floor(t_prev / kReportStep + 1.0) * kReportStep;
             b <= t_now + 1e-6 && b <= end_time_s + 1e-6; b += kReportStep) {
            out.times.push_back(b);
            // H17 diagnostic: phi and its truncated reconstruction P P^T phi at
            // the nodes named in H17_PHI_NODES, one line per node per report.
            if (rom && phi_out.is_open()) {
                const auto& phi = eng->rom1dFlowResponse();
                const auto& P = rom->basis->P;
                const std::size_t nn = static_cast<std::size_t>(rom->n_nodes);
                const std::size_t nk = static_cast<std::size_t>(rom->n_kept);
                if (phi.size() == nn) {
                    std::vector<double> c(nk, 0.0);
                    for (std::size_t j = 0; j < nk; ++j)
                        for (std::size_t i = 0; i < nn; ++i) c[j] += P[j * nn + i] * phi[i];
                    for (const auto& [nm, ui] : phi_nodes) {
                        const int ai = rom->full_to_active[static_cast<std::size_t>(ui)];
                        if (ai < 0) continue;
                        double rec = 0.0;
                        for (std::size_t j = 0; j < nk; ++j) rec += P[j * nn + static_cast<std::size_t>(ai)] * c[j];
                        phi_out << b << ' ' << nm << ' ' << phi[static_cast<std::size_t>(ai)] << ' ' << rec << ' '
                                << (ctx.nodes.head[static_cast<std::size_t>(ui)] - ctx.nodes.invert_elev[static_cast<std::size_t>(ui)])
                                << ' ' << (rom->q95[static_cast<std::size_t>(ai)] - rom->q05[static_cast<std::size_t>(ai)]) << '\n';
                    }
                }
            }
            for (const auto& [nm, ui] : idx) {
                out.heads[nm].push_back(ctx.nodes.head[static_cast<std::size_t>(ui)]);
                if (rom && rom->is_ready()) {
                    const int ai = rom->full_to_active[static_cast<std::size_t>(ui)];
                    out.q05[nm].push_back(ai >= 0 ? rom->q05[static_cast<std::size_t>(ai)] : NAN);
                    out.q95[nm].push_back(ai >= 0 ? rom->q95[static_cast<std::size_t>(ai)] : NAN);
                    const auto& fl = eng->rom1dBandFlags();
                    out.flag[nm].push_back(ai >= 0 && static_cast<std::size_t>(ai) < fl.size()
                                           ? fl[static_cast<std::size_t>(ai)] : -1);
                }
            }
        }
        t_prev = t_now;
        if (t_now >= end_time_s) break;
    }
    out.ok = !out.times.empty();
    cleanup();
    return out;
}

}  // namespace

// H17: one-member dump mode for large networks (Bellinge), so members can run
// as parallel processes. Skipped unless H17_ONE_INP is set. Env:
//   H17_ONE_INP   deck path      H17_ONE_MULT  conduit n multiplier
//   H17_ONE_ROM   1 = add the ROM spec and dump q05/q95/flag
//   H17_ONE_OUT   dump path      H17_WINDOW_H  window (h)
// Dump: line 1 = report times; then per node "name kind v1 v2 ..." with
// kind in {h, q05, q95, flag}. Scored by tests/regression/score_h17_mc.py.
TEST(RomCoverageStructures, OneMemberDump) {
    const char* inp = std::getenv("H17_ONE_INP");
    if (!inp) GTEST_SKIP() << "H17_ONE_INP not set";
    const std::string base = readFile(inp);
    ASSERT_FALSE(base.empty()) << inp;
    const double mult = std::getenv("H17_ONE_MULT") ? std::atof(std::getenv("H17_ONE_MULT")) : 1.0;
    const bool with_rom = std::getenv("H17_ONE_ROM") && std::atoi(std::getenv("H17_ONE_ROM")) != 0;
    const char* outp = std::getenv("H17_ONE_OUT");
    ASSERT_NE(outp, nullptr) << "H17_ONE_OUT not set";
    double window_h = 24.0;
    if (const char* v = std::getenv("H17_WINDOW_H")) window_h = std::atof(v);
    const auto names = scoredNodes(base);
    const std::string tag = std::string("one_") + (with_rom ? "rom_" : "mc_") + std::to_string(mult);
    RunResult r = runCase(deck(base, mult, with_rom), tag, with_rom, names, window_h * 3600.0);
    ASSERT_TRUE(r.ok);
    std::ofstream f(outp);
    f.precision(10);
    for (double t : r.times) f << t << ' ';
    f << '\n';
    for (const auto& [nm, hv] : r.heads) {
        f << nm << " h";
        for (double v : hv) f << ' ' << v;
        f << '\n';
        if (with_rom && r.q05.count(nm)) {
            f << nm << " q05"; for (double v : r.q05.at(nm)) f << ' ' << v; f << '\n';
            f << nm << " q95"; for (double v : r.q95.at(nm)) f << ' ' << v; f << '\n';
            f << nm << " flag"; for (int v : r.flag.at(nm)) f << ' ' << v; f << '\n';
        }
    }
    std::printf("[H17-one] %s times=%zu nodes=%zu -> %s\n", tag.c_str(), r.times.size(), r.heads.size(), outp);
}

TEST(RomCoverageStructures, Example8ConduitManningBandsVsMonteCarlo) {
    if (std::getenv("H17_ONE_INP")) GTEST_SKIP() << "one-member dump mode";
    const std::string base = readFile(std::string(REGRESSION_DATA_DIR) + "/Example8.inp");
    ASSERT_FALSE(base.empty()) << "Example8.inp not found under " << REGRESSION_DATA_DIR;
    const auto names = scoredNodes(base);
    ASSERT_GE(names.size(), 20u);

    double window_h = 12.0;
    if (const char* v = std::getenv("H17_WINDOW_H")) window_h = std::atof(v);
    const double end_s = window_h * 3600.0;

    std::vector<RunResult> mc(kMcRuns);
    for (int i = 0; i < kMcRuns; ++i) {
        const double mult = (1.0 - kPert) + (i + 0.5) / kMcRuns * 2.0 * kPert;
        mc[i] = runCase(deck(base, mult, false), "mc" + std::to_string(i), false, names, end_s);
        ASSERT_TRUE(mc[i].ok) << "MC run " << i;
    }
    RunResult rom = runCase(deck(base, 1.0, true), "rom", true, names, end_s);
    ASSERT_TRUE(rom.ok) << "ROM run";
    ASSERT_FALSE(rom.q05.empty()) << "ROM produced no quantiles";

    std::size_t n_s = rom.times.size();
    for (const auto& r : mc) n_s = std::min(n_s, r.times.size());

    struct NodeStat { int n = 0; double cov = 0; std::vector<double> ratios; };
    std::map<std::string, NodeStat> per;
    double cov_sum = 0; int n_tot = 0;
    double mat_cov = 0; int n_mat = 0; std::vector<double> mat_ratios;   // MC width >= 0.05 ft
    std::map<int, NodeStat> by_flag;   // H17 step 4: score per band validity flag
    std::vector<double> all_ratios;
    for (const auto& [nm, q05] : rom.q05) {
        const auto& q95 = rom.q95.at(nm);
        auto& st = per[nm];
        for (std::size_t k = 0; k < n_s; ++k) {
            if (rom.times[k] <= 60.0 || std::isnan(q05[k])) continue;
            std::vector<double> hv(kMcRuns);
            for (int i = 0; i < kMcRuns; ++i) hv[i] = mc[i].heads.at(nm)[k];
            std::sort(hv.begin(), hv.end());
            const double mcw = mcq::quantileMidpoint(hv, 0.95) - mcq::quantileMidpoint(hv, 0.05);
            if (mcw < 1e-3) continue;   // 1e-3 ft: MC has no resolvable spread here
            const double c = mcq::intervalCoverage(hv, q05[k], q95[k]);
            st.cov += c; ++st.n;
            cov_sum += c; ++n_tot;
            const double r = (q95[k] - q05[k]) / mcw;
            st.ratios.push_back(r);
            all_ratios.push_back(r);
            if (mcw >= 0.05) { mat_cov += c; ++n_mat; mat_ratios.push_back(r); }
            const int fg = (rom.flag.count(nm) && k < rom.flag.at(nm).size()) ? rom.flag.at(nm)[k] : -1;
            auto& fs = by_flag[fg];
            fs.cov += c; ++fs.n; fs.ratios.push_back(r);
        }
    }
    ASSERT_GT(n_tot, 0) << "MC produced no resolvable spread";
    if (const char* tn = std::getenv("H17_NODE")) {
      std::stringstream nl(tn);
      std::string nm;
      while (std::getline(nl, nm, ',')) {
        if (rom.q05.count(nm)) {
            for (std::size_t k = 0; k < n_s; k += 15) {
                std::vector<double> hv(kMcRuns);
                for (int i = 0; i < kMcRuns; ++i) hv[i] = mc[i].heads.at(nm)[k];
                std::sort(hv.begin(), hv.end());
                std::printf("[H17-trace] %s t=%6.0f det=%.3f mc[q05,q95]=[%.3f,%.3f] rom[q05,q95]=[%.3f,%.3f]\n",
                            nm.c_str(), rom.times[k], rom.heads.at(nm)[k],
                            mcq::quantileMidpoint(hv, 0.05), mcq::quantileMidpoint(hv, 0.95),
                            rom.q05.at(nm)[k], rom.q95.at(nm)[k]);
            }
        }
      }
    }
    auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; };

    std::printf("[H17-node] %-6s %5s %7s %8s\n", "node", "n", "cov", "w-med");
    for (const auto& [nm, st] : per) {
        if (st.n == 0) { std::printf("[H17-node] %-6s %5d   (no resolvable MC spread)\n", nm.c_str(), 0); continue; }
        std::printf("[H17-node] %-6s %5d %7.3f %8.3f\n", nm.c_str(), st.n, st.cov / st.n, med(st.ratios));
    }
    for (const auto& [fg, st] : by_flag)
        std::printf("[H17-flag] %-13s samples=%5d member-coverage=%.3f width-med=%.3f\n",
                    fg < 0 ? "?" : openswmm::SWMMEngine::rom1dBandFlagName(static_cast<std::uint8_t>(fg)),
                    st.n, st.cov / st.n, med(st.ratios));
    // H17 step 4: on this fixture the `steady` samples are calibrated; gate it.
    // (Not a network-independent promise: Bellinge steady is 0.714 / 1.575,
    // VALIDATION.md "H17 §8".)
    ASSERT_TRUE(by_flag.count(0) && by_flag.at(0).n > 0) << "no steady-flagged samples";
    EXPECT_TRUE(mcq::reportCalibration("Example 8, steady-flagged bands", by_flag.at(0).cov / by_flag.at(0).n,
                                       med(by_flag.at(0).ratios), "-"))
        << "bands flagged `steady` must stay calibrated (measured 0.944 / 1.350)";
    const double cov = cov_sum / n_tot, wmed = med(all_ratios);
    std::printf("[H17] Example 8 samples=%d member-coverage=%.3f width-med=%.3f\n", n_tot, cov, wmed);
    mcq::reportCalibration("Example 8 conduit Manning (structures)", cov, wmed,
                           "-");
    EXPECT_GE(wmed, 0.3) << "Example 8 band collapsed (measured 1.102)";
    EXPECT_LE(wmed, 3.0) << "Example 8 band blew up (measured 1.102)";
    if (n_mat > 0)
        std::printf("[H17] material spreads only (MC width >= 0.05 ft): samples=%d member-coverage=%.3f width-med=%.3f\n",
                    n_mat, mat_cov / n_mat, med(mat_ratios));
}
