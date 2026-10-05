/**
 * @file RunoffElasticity.hpp
 * @brief SR-6: per-subcatchment elasticity of runoff to rainfall, from two
 *        persistent perturbed runoff states advanced through the production
 *        kernel.
 *
 * Why. The soft-rain spread was mapped as if runoff scaled one-for-one with
 * rain. On the rising limb of Manning overland flow (Q ~ d^(5/3)) the delivered
 * runoff responds with an elasticity between 1 and 5/3, and the ROM band came
 * out ~0.71 of the brute-force band for exactly that reason (VALIDATION.md,
 * "SR-6 step 1"). A single-step trial cannot see this: most of the
 * sensitivity is in the accumulated ponded depth, not in one step's rain. So
 * two copies of the runoff state are carried persistently at rain*(1-delta)
 * and rain*(1+delta), advanced every runoff step with the same arguments as
 * the deterministic call, and
 *
 *     E_s(t) = ln(q_plus / q_minus) / ln((1+delta)/(1-delta))
 *
 * is the local log-slope of runoff to rain for subcatchment s. The spread
 * mapping multiplies the gage's relative spread by E_s. Exact for a power law
 * at any delta; 1.0 when either member is dry.
 *
 * What it does NOT include: LID surface runoff (added to subcatches.runoff by
 * the engine after execute()) and subcatchment cascading (members read the
 * DETERMINISTIC runon). Both are documented approximations.
 *
 * Cost: two extra RunoffSolver::execute() calls per RUNOFF step (minutes, not
 * routing steps) plus vector copies of the solver state.
 *
 * Bit-identity: the deterministic kernel is never modified; members run on a
 * swapped-in copy of the solver state and the nine ctx.subcatches arrays
 * execute() writes are restored afterwards.
 *
 * @ingroup engine_uncertainty
 */
#ifndef OPENSWMM_ENGINE_UNCERTAINTY_RUNOFF_ELASTICITY_HPP
#define OPENSWMM_ENGINE_UNCERTAINTY_RUNOFF_ELASTICITY_HPP

#include "../hydrology/Runoff.hpp"

#include <vector>

namespace openswmm { struct SimulationContext; }

namespace openswmm::uncertainty {

class RunoffElasticityProbe {
public:
    /// Relative rain perturbation of the two members (default 0.1).
    double delta = 0.10;
    /// Elasticity clamp (a power law on the rising limb is <= 5/3; 3 is a
    /// generous guard against a nearly-dry member making the log-ratio blow up).
    double max_elasticity = 3.0;

    /// Advance both perturbed members one runoff step (call BEFORE the
    /// deterministic execute() with the same arguments) and refresh
    /// elasticity(). Members are seeded from the deterministic state on the
    /// first call. Leaves the solver and ctx exactly as found.
    void step(runoff::RunoffSolver& solver, SimulationContext& ctx, double dt,
              double evap_rate, double infil_factor, double recovery_factor, int month);

    /// Per-subcatchment elasticity, 1.0 until the first step (and wherever a
    /// member is dry). Sized to n_subcatch after the first step.
    const std::vector<double>& elasticity() const noexcept { return elasticity_; }
    bool seeded() const noexcept { return seeded_; }
    int  steps() const noexcept { return n_steps_; }

private:
    bool seeded_ = false;
    int  n_steps_ = 0;
    runoff::RunoffSolver::State minus_, plus_;
    std::vector<double> q_minus_, q_plus_, elasticity_;
};

} // namespace openswmm::uncertainty

#endif
