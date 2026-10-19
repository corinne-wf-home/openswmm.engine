# Full-order additive channel (H19) — design note

Status: design fixed 2026-10-16 (owner chose option 2 after H18 A.4b,
VALIDATION.md "H18" §10–§10c); the default (and only) additive path since A.5,
2026-10-16. It replaces the modal forcing of every additive
uncertainty: soft rain (gage planes and CORR_LEN fields), the built-in runoff
multiplier, and registered `FORCING_VECTOR` sources such as `INFLOW`. The
Manning channel (modal, H17 target + H18 rates) is unchanged.

## 1. Why

Two measured failures of the modal additive forcing:

1. **The downstream transfer.** A member with more rain carries more flow in
   every pipe below the extra water. H18's stable symmetric operator drops that
   transfer (SR-5 J5 at 0.3–2 % of the MC). Relaxing toward the steady
   flow-response target restores it, but delivers it at once (J5 3–25× wide on a
   storage-dominated network). Only dynamics that contain the one-way transfer
   are right both before and after the flows settle, and projecting a one-way
   operator onto a few modes was unstable in October.
2. **The source.** Soft rain was driven by the deterministic head rate dh/dt
   times the relative spread. Under steady rain in steady flow dh/dt = 0, so no
   band forms although the rain uncertainty changes every flow.

## 2. The model

For each additive source p the channel carries one full-order head deviation
field Ψ_p (all active nodes, no mode truncation), shared by all members:

    S dΨ_p/dt = −J Ψ_p + δq_p

- **S:** node storage (storage-curve area + half the top width of attached
  conduits, floored at MIN_SURFAREA; the dynamic-wave rule, from published
  state).
- **J:** the linearized flow Jacobian of the flow-response sweep: per conduit
  δQ = Ju·Ψ_up + Jd·Ψ_dn, and per structure δQ = |dQ/dH|·Ψ_up, oriented by
  flow. This includes the one-way transfer, the Froude relay and the full-pipe
  bound exactly as the Manning target uses them. Outfalls hold zero deviation.
  A node overflowing at its rim is pinned at zero (the extra water spills).
- **δq_p:** the source's extra inflow per unit member coefficient (§4).

Member i's additive deviation is Σ_p κ_ip·Ψ_p, added to its reconstructed head
next to the Manning deviation:
- comonotone soft rain: κ = c_i (plane A) and c_i^B (plane B);
- CORR_LEN reduced: one Ψ_m per spatial mode, κ = a_im;
- CORR_LEN materialized: one Ψ per member, κ = identity;
- runoff multiplier: κ = scale_i − 1;
- forcing vectors: κ = θ_i − 1.

**Time step.** Backward Euler, (S/dt + J)·Ψ_new = (S/dt)·Ψ_old + δq. The
matrix is a strictly column-diagonally-dominant M-matrix: each column's
off-diagonal sum is exactly its diagonal minus S/dt, and the column sums are
volume conservation. So it is nonsingular, the step is unconditionally stable,
and Gauss–Seidel converges. Sweeps alternate upstream-first and
downstream-first, warm-started from Ψ_old, to a relative tolerance of 1e-10,
with at most 200 sweeps (the count is reported).

**Limits of the linear model (documented, not fixed):**
- Members do not scale J by 1/mm, so the roughness × rain interaction, a
  second-order effect, is dropped.
- Conduits with no flow carry no conductance, so a still pool simply stores
  its extra inflow. That is the physics of a pool without outflow.

## 3. What it costs

One Gauss–Seidel solve per plane per routing step, at O(edges) per sweep.
- Gage soft rain: two planes at most.
- CORR_LEN reduced: K_s planes.
- CORR_LEN materialized: M planes.
- Reconstruction adds Σ_p κ_ip Ψ_p at report boundaries: O(M·P·n).

## 4. Sources (extra inflow per unit coefficient, cfs)

- **Soft rain, plane A or B:** Σ over the node's contributing subcatchments of
  runoff_s × relative spread of its gage × rain→runoff elasticity E_s (SR-6),
  restricted to that plane's family. The relative spread is CV, or SD or
  HALFRANGE over the gage rainfall, as today. No area weighting is needed,
  because each subcatchment contributes its own runoff.
- **CORR_LEN:** the plane-A source multiplied pointwise by the spatial field (ψ_m
  or W_i).
- **Runoff multiplier:** the node's subcatchment runoff.
- **Forcing vectors:** the node's lateral inflow, so "INFLOW 0.30" means
  inflow 30 % larger. Before, it meant head changes 30 % stronger.

## 5. What the ROM does

`SpectralROM1D` still generates the per-member
coefficients (`softCoefficients(family)` returns c_i or c_i^B). It applies no
additive modal forcing (the modal path was deleted on 2026-10-19), and its
reconstruction (`computeQuantiles`, `reconstructHead`) adds Σ_p κ_ip Ψ_p from
`setAdditivePlanes(P, Ψ, κ)`. The Manning channel, the phase coordinate and the
flags are unchanged.

## 6. Validation (gate)

| Gate | Requirement |
|---|---|
| SR-5, CL-1e (DW and FV) | calibrated |
| **New:** steady rain on steady flow (SR-7) | calibrated; the old channel gives no band |
| Every Manning gate (chains DW/FV, Example 8, Bellinge MC) | unchanged (no additive source) |
| Unit tests | steady state = flow-response sweep with accumulated throughflow; early response local (no instant transfer); volume conserved; zero coefficient exact |

After this gate passes, H18 A.5 resumes: H18 as the default, then K1d, the
advective hook and the old modal additive forcing are deleted.
