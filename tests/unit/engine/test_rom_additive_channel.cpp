/**
 * @file test_rom_additive_channel.cpp
 * @brief PR H19 — full-order linearized additive channel (RomAdditiveChannel.hpp).
 *
 * S dPsi/dt = -J Psi + dq, backward Euler, Gauss-Seidel. Each case checks a
 * property the modal additive forcing could not hold: the steady response
 * carries the accumulated extra flow, the early response stays local, and
 * volume is conserved.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "uncertainty/RomAdditiveChannel.hpp"

using openswmm::uncertainty::AdditiveSystem;
using openswmm::uncertainty::FlowEdge;

namespace {

FlowEdge conduit(int up, int dn, double ju, double jd, double q = 1.0) {
    FlowEdge e;
    e.up = up; e.dn = dn; e.conduit = true; e.q = q; e.ju = ju; e.jd = jd;
    return e;
}

// Chain 0 -> 1 -> ... -> n-1 -> outfall, one conduit per node.
std::vector<FlowEdge> chain(int n, double ju, double jd) {
    std::vector<FlowEdge> e;
    for (int i = 0; i < n; ++i) e.push_back(conduit(i, i + 1 < n ? i + 1 : -1, ju, jd));
    return e;
}

}  // namespace

TEST(AdditiveChannel, SteadyStateCarriesTheAccumulatedExtraFlow) {
    // Kinematic chain (Jd = 0): at steady state conduit i carries the extra
    // inflow of every node above it, so Psi_i = (i + 1) * dq / Ju. This is the
    // transfer the symmetric modal operator dropped (SR-5 J5).
    const int n = 5;
    const double ju = 0.4, dq = 1.0;
    const auto edges = chain(n, ju, 0.0);
    std::vector<double> S(n, 100.0), q(n, dq), psi(n, 0.0);
    AdditiveSystem sys;
    sys.build(n, edges, S.data(), 1.0e9);
    for (int k = 0; k < 5; ++k) sys.step(q.data(), psi.data());
    for (int i = 0; i < n; ++i)
        EXPECT_NEAR(psi[static_cast<std::size_t>(i)], (i + 1) * dq / ju, 1.0e-6) << "node " << i;
}

TEST(AdditiveChannel, EarlyResponseStaysLocal) {
    // Storage-dominated (the SR-5 regime): one second after extra water enters
    // node 0, node 3 has barely moved. Relaxing toward the steady target moved
    // it at once (A.4b, VALIDATION "H18" s10b).
    const int n = 5;
    const auto edges = chain(n, 0.4, 0.0);
    std::vector<double> S(n, 2.0e4), q(n, 0.0), psi(n, 0.0);
    q[0] = 1.0;
    AdditiveSystem sys;
    sys.build(n, edges, S.data(), 1.0);
    sys.step(q.data(), psi.data());
    EXPECT_NEAR(psi[0], 1.0 / (2.0e4 + 0.4), 1.0e-12);   // essentially dq*dt/S
    EXPECT_LT(psi[3], 1.0e-12 * psi[0]);
}

TEST(AdditiveChannel, VolumeIsConserved) {
    // Each step: change in stored volume + volume leaving at the outfall =
    // extra inflow * dt (backward Euler, converged).
    const int n = 4;
    const double ju = 0.7, jd = -0.2, dt = 30.0;
    const auto edges = chain(n, ju, jd);
    std::vector<double> S = {50.0, 80.0, 120.0, 60.0}, q = {0.3, 0.1, 0.0, 0.2}, psi(n, 0.0);
    AdditiveSystem sys;
    sys.build(n, edges, S.data(), dt);
    for (int k = 0; k < 20; ++k) {
        std::vector<double> old = psi;
        sys.step(q.data(), psi.data(), 500, 1.0e-14);
        double dV = 0.0, in = 0.0;
        for (int i = 0; i < n; ++i) {
            dV += S[static_cast<std::size_t>(i)] * (psi[static_cast<std::size_t>(i)] - old[static_cast<std::size_t>(i)]);
            in += q[static_cast<std::size_t>(i)] * dt;
        }
        const double out = ju * psi[n - 1] * dt;   // last conduit to the outfall (Psi_outfall = 0)
        EXPECT_NEAR(dV + out, in, 1.0e-9 * in) << "step " << k;
    }
}

TEST(AdditiveChannel, CoupledSystemMatchesDenseSolve) {
    // Subcritical chain with downstream head coupling plus a structure into
    // the last node: Gauss-Seidel equals Gaussian elimination of the same
    // backward-Euler system.
    const int n = 4;
    std::vector<FlowEdge> edges = {conduit(0, 1, 0.9, -0.5), conduit(1, 2, 0.8, -0.4),
                                   conduit(2, -1, 0.6, -0.1)};
    FlowEdge st; st.up = 3; st.dn = 2; st.conduit = false; st.d = 0.3; st.q = 0.5;
    edges.push_back(st);
    std::vector<double> S = {10.0, 20.0, 15.0, 5.0}, q = {0.1, 0.0, 0.2, 0.05}, psi(n, 0.0);
    const double dt = 5.0;
    AdditiveSystem sys;
    sys.build(n, edges, S.data(), dt);
    sys.step(q.data(), psi.data(), 1000, 1.0e-15);
    // Dense matrix: row i = S_i/dt x_i + sum_out dQ - sum_in dQ.
    double A[4][4] = {}, b[4];
    for (int i = 0; i < n; ++i) { A[i][i] = S[static_cast<std::size_t>(i)] / dt; b[i] = q[static_cast<std::size_t>(i)]; }
    for (const auto& e : edges) {
        auto add = [&](int row, int c, double v) { if (row >= 0 && c >= 0) A[row][c] += v; };
        if (e.conduit) {
            add(e.up, e.up, e.ju); add(e.up, e.dn, e.jd);
            add(e.dn, e.up, -e.ju); add(e.dn, e.dn, -e.jd);
        } else {
            add(e.up, e.up, e.d); add(e.dn, e.up, -e.d);
        }
    }
    for (int k = 0; k < n; ++k)
        for (int r = k + 1; r < n; ++r) {
            const double f = A[r][k] / A[k][k];
            for (int c = k; c < n; ++c) A[r][c] -= f * A[k][c];
            b[r] -= f * b[k];
        }
    double x[4];
    for (int r = n - 1; r >= 0; --r) {
        double s = b[r];
        for (int c = r + 1; c < n; ++c) s -= A[r][c] * x[c];
        x[r] = s / A[r][r];
    }
    for (int i = 0; i < n; ++i) EXPECT_NEAR(psi[static_cast<std::size_t>(i)], x[i], 1.0e-12) << "node " << i;
}

TEST(AdditiveChannel, PinnedNodeStaysAtZeroAndPassesNothingOn) {
    const int n = 3;
    const auto edges = chain(n, 0.5, 0.0);
    std::vector<double> S(n, 10.0), q = {1.0, 0.0, 0.0}, psi(n, 0.0);
    std::vector<char> pinned = {0, 1, 0};
    AdditiveSystem sys;
    sys.build(n, edges, S.data(), 1.0e6, &pinned);
    for (int k = 0; k < 3; ++k) sys.step(q.data(), psi.data());
    EXPECT_EQ(psi[1], 0.0);
    EXPECT_NEAR(psi[2], 0.0, 1.0e-14);
    EXPECT_GT(psi[0], 0.0);
}

TEST(AdditiveChannel, ZeroSourceStaysExactlyZero) {
    const int n = 5;
    const auto edges = chain(n, 0.4, -0.1);
    std::vector<double> S(n, 30.0), q(n, 0.0), psi(n, 0.0);
    AdditiveSystem sys;
    sys.build(n, edges, S.data(), 10.0);
    for (int k = 0; k < 10; ++k) sys.step(q.data(), psi.data());
    for (double v : psi) EXPECT_EQ(v, 0.0);
}
