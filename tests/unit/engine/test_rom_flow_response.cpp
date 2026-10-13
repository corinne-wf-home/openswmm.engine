/**
 * @file test_rom_flow_response.cpp
 * @brief PR H17 — the 1D ROM's Manning sensitivity field phi = dh/d(ln n),
 *        solved upstream from the outfalls (RomFlowResponse.hpp).
 *
 * Each case builds the oriented edge list by hand and checks phi against the
 * closed form the linearization implies.
 */
#include <gtest/gtest.h>

#include <vector>

#include "uncertainty/RomFlowResponse.hpp"

using openswmm::uncertainty::FlowEdge;
using openswmm::uncertainty::solveFlowResponse;

namespace {

// Conduit edge from the linearization: Ju = q(1/(2dh) + kfac), Jd = -q/(2dh).
FlowEdge conduit(int up, int dn, double q, double dh, double kfac) {
    FlowEdge e;
    e.up = up; e.dn = dn; e.conduit = true;
    e.q = q; e.ju = q * (0.5 / dh + kfac); e.jd = -q * 0.5 / dh;
    return e;
}
FlowEdge structure(int up, int dn, double d) {
    FlowEdge e;
    e.up = up; e.dn = dn; e.conduit = false; e.d = d;
    return e;
}

}  // namespace

// Full pipes (kfac = 0): friction loss ~ n^2, so each pipe adds 2*dh and the
// responses accumulate from the outfall upstream.
TEST(RomFlowResponse, FullPipeChainAccumulatesTwiceTheFrictionLoss) {
    const std::vector<FlowEdge> e = {conduit(0, 1, 3.0, 1.5, 0.0),
                                     conduit(1, 2, 3.0, 2.0, 0.0),
                                     conduit(2, -1, 3.0, 0.5, 0.0)};
    std::vector<double> phi;
    EXPECT_EQ(solveFlowResponse(3, e, phi), 1);
    EXPECT_NEAR(phi[2], 2.0 * 0.5, 1e-12);
    EXPECT_NEAR(phi[1], 2.0 * (0.5 + 2.0), 1e-12);
    EXPECT_NEAR(phi[0], 2.0 * (0.5 + 2.0 + 1.5), 1e-12);
}

// Steep free surface (|dh| >> depth, wide section K ~ y^(5/3)): phi tends to
// 1/kfac = 0.6*y at every node, with no relay from downstream.
TEST(RomFlowResponse, SteepFreeSurfaceGivesSixTenthsOfDepth) {
    const double y = 0.8, kfac = (5.0 / 3.0) / y, dh = 1.0e6;
    const std::vector<FlowEdge> e = {conduit(0, 1, 2.0, dh, kfac),
                                     conduit(1, -1, 2.0, dh, kfac)};
    std::vector<double> phi;
    solveFlowResponse(2, e, phi);
    EXPECT_NEAR(phi[0], 0.6 * y, 1e-5);
    EXPECT_NEAR(phi[1], 0.6 * y, 1e-5);
}

// A node with no outflow (still water, or the tip of a filling front) has no
// steady state for roughness to shift: phi = 0, and it pins its feeders.
TEST(RomFlowResponse, SinkAndStillNodesAreZero) {
    const std::vector<FlowEdge> e = {conduit(0, 1, 1.0, 1.0, 0.0)};
    std::vector<double> phi;
    solveFlowResponse(3, e, phi);
    EXPECT_EQ(phi[1], 0.0);   // sink
    EXPECT_EQ(phi[2], 0.0);   // no link at all
    EXPECT_NEAR(phi[0], 2.0, 1e-12);
}

// A free structure passes a flow set by its upstream head alone and carries
// no Manning source: draining only through it gives phi = 0, and adding it
// beside a conduit absorbs part of the rise (overflow takes the extra head).
TEST(RomFlowResponse, StructuresAbsorbRatherThanRelay) {
    std::vector<double> phi;
    solveFlowResponse(2, {structure(0, 1, 5.0), conduit(1, -1, 1.0, 1.0, 0.0)}, phi);
    EXPECT_EQ(phi[0], 0.0);
    EXPECT_NEAR(phi[1], 2.0, 1e-12);

    const auto c = conduit(0, -1, 1.0, 1.0, 0.0);
    solveFlowResponse(1, {c}, phi);
    const double alone = phi[0];
    solveFlowResponse(1, {c, structure(0, -1, 0.5)}, phi);
    EXPECT_NEAR(phi[0], c.q / (c.ju + 0.5), 1e-12);
    EXPECT_LT(phi[0], alone);
}

// A flow split: the node's single row balances both outgoing pipes.
TEST(RomFlowResponse, SplitBalancesBothOutflows) {
    const auto a = conduit(0, -1, 1.0, 1.0, 0.0);
    const auto b = conduit(0, -1, 3.0, 1.0, 0.0);
    std::vector<double> phi;
    solveFlowResponse(1, {a, b}, phi);
    EXPECT_NEAR(phi[0], (a.q + b.q) / (a.ju + b.ju), 1e-12);
}

// A flow cycle cannot be ordered; Gauss-Seidel settles it to the same answer
// the node equations give.
TEST(RomFlowResponse, CycleConvergesToTheNodeEquations) {
    const std::vector<FlowEdge> e = {conduit(0, 1, 1.0, 1.0, 0.2),
                                     conduit(1, 0, 0.5, 1.0, 0.2),
                                     conduit(1, -1, 1.0, 1.0, 0.2)};
    std::vector<double> phi;
    EXPECT_GT(solveFlowResponse(2, e, phi), 1);
    // Node 0: phi0 = (q0 - jd0*phi1)/ju0; node 1: phi1 = (q1+q2 - jd1*phi0)/(ju1+ju2).
    const double r0 = (e[0].q - e[0].jd * phi[1]) / e[0].ju;
    const double r1 = (e[1].q + e[2].q - e[1].jd * phi[0]) / (e[1].ju + e[2].ju);
    EXPECT_NEAR(phi[0], r0, 1e-8);
    EXPECT_NEAR(phi[1], r1, 1e-8);
}

// A junction overflowing at its rim has its head pinned by the spill: its
// response is zero and it relays nothing upstream (Bellinge G70F143,
// 2026-10-13: at the rim with zero Monte Carlo spread, unpinned it read twice
// the drop below it).
TEST(RomFlowResponse, OverflowingNodeIsPinned) {
    const std::vector<FlowEdge> e = {conduit(0, 1, 3.0, 1.0, 0.0),
                                     conduit(1, -1, 3.0, 20.0, 0.0)};
    std::vector<double> phi;
    solveFlowResponse(2, e, phi);
    EXPECT_NEAR(phi[1], 40.0, 1e-9);
    EXPECT_NEAR(phi[0], 42.0, 1e-9);
    const std::vector<char> pinned = {0, 1};
    solveFlowResponse(2, e, phi, 50, &pinned);
    EXPECT_EQ(phi[1], 0.0);
    EXPECT_NEAR(phi[0], 2.0, 1e-9);
}
