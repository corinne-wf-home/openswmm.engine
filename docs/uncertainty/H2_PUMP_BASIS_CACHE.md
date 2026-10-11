# H2 — Pump-cycle basis cache: design note and premise check

**Status**: premise measured and falsified on this line (2026-10-10). The cache
is **not implemented**; it is deferred behind a prerequisite, **H17**
(non-conduit links in the 1D ROM operator). This note is the design one-pager
the checklist asks for, including the §7.1 forward question, so H2 can be built
directly once H17 lands.

**Branch**: `hsym2/prh2-pump-basis-cache` (stacked on H6).

---

## 1. What H2 assumed

> Every pump on/off transition changes H discontinuously → cold re-solve every
> cycle. A two-state pump cycling every few minutes pays full Lanczos cost
> twice per cycle forever, though the system only ever visits two operator
> configurations.

The fix it proposed: an LRU of `{P, λ, conduit weights}` keyed on the discrete
control state, restored on a state flip instead of re-solving.

## 2. What the code does

1. **The 1D ROM graph is built from conduits only.** `buildROM1D()` and
   `updateBasis()` take their edges from `DWSolver::HSnapshot`, which exports
   `0.5·dt·dqdh` for conduits and nothing else. Pumps, orifices, weirs and
   outlets are not edges of the operator.
2. **So a pump switch is not an operator event.** Turning a pump off does not
   remove an edge; it changes flows, and with them the conduit weights near the
   wet well, smoothly over the following steps. The weighted Laplacian never
   sees the pump's on/off state.
3. **Coefficients already survive any basis change.** A rebuild re-projects the
   ensemble through the full overlap `R = P_newᵀ·P_old`, so sign flips and mode
   reorderings are absorbed. The spec's "2× spread jump per swap without sign
   alignment" cannot happen on this path; a cache hit would reuse the same
   re-projection.

## 3. Measurement: EPA Example 8 (one type-III pump, 12 h)

Probe: the engine driven through the C API with `[UNCERTAINTY] 1D MANNINGS_N
0.20`, logging every pump transition and every forced-cold basis rebuild.

| Sim time (s) | Event |
|---|---|
| 78 – 550 | forced-cold rebuilds #1–#7 (no pump activity; network wetting up from dry) |
| 2,016 | pump ON |
| 2,132 | forced-cold rebuild #8 (one node above crown) |
| 8,423 | pump OFF |
| 24,890 | pump ON |
| 25,919 | pump OFF |

Totals: 4 pump transitions, 125 basis-update attempts, 0 failed, 8 forced cold.
Seven of the eight cold rebuilds are the dry-start wet-up (H1's edge-drift
trigger, as designed). Three of the four pump transitions trigger nothing; the
one that is followed by a cold rebuild is followed 116 s later, coincident with
a node surcharging, which is H1's surcharge trigger.

**Structural consequence, worth more than the cache.** Example 8's conduit-only
graph splits into two components (26 and 5 nodes), joined in the real network
only by the pump. Each has its own outfall, so both stay grounded and the
solve is well-posed. But a deviation upstream of the pump cannot reach the force
main in the ROM, and weir/orifice overflow paths are likewise invisible.

## 4. Conclusion

As specified, H2 would cache operators keyed on a state the operator does not
depend on. It would rarely be consulted, and when consulted would return a basis
no different from the warm-started one. Building it now adds an LRU, a key
extractor and a staleness check with no measurable consumer and no test that
could show a benefit. **Not built.**

## 5. Prerequisite: H17 — non-conduit links in the ROM operator

The data already exists. `HydStructures.cpp` sets `links.dqdh` for pumps
(curve slope), orifices and weirs every iteration, exactly as legacy does. H17
would:

- export non-conduit links in the snapshot alongside conduits (`0.5·dt·dqdh`,
  node pair), with a per-link type tag;
- include them as edges in `NetworkLaplacian1D::buildWeighted`, using the same
  `1e-6` weight floor so topology stays fixed while a pump is off (an "off"
  pump is a near-zero edge, not a deleted one, so the active-node-count guard in
  `updateBasis()` never trips);
- leave the H14 directional operator's Froude gate on conduits only (Froude is
  undefined for a pump).

Effect: a pump transition becomes a discrete weight change, which is the event
H1's edge-drift trigger already detects and H2's cache can amortize. This
**changes ROM numerics on every network with non-conduit links**, so it needs
its own MC validation on a pumped fixture (Example 8 is the natural one, with
member-level brute-force runs), and the Bellinge lock re-run. That is an owner
decision, not a mechanical follow-up.

## 6. Cache design, for when H17 lands

- **Key**: the vector of discrete control states over non-conduit links —
  pump on/off, and orifice/weir/outlet settings quantized to 1e-3 relative —
  hashed. Read from `links.setting` / pump status at snapshot time; read-only.
- **Entry**: a value object `{key, P, λ, weights_at_solve}`. Independent of the
  ROM's live `basis` pointer, so entries can outlive a swap.
- **Lookup on a key change**: hit → staleness check (the H1 edge-drift
  criterion between `weights_at_solve` and the live weights; fail ⇒ treat as a
  miss) → install `P, λ` → re-project the ensemble through `R = P_newᵀ·P_old`
  (§2 point 3). Miss → the existing H1 path (cold Lanczos), then insert.
- **Eviction**: LRU, cap 4. Memory per entry `k·N` doubles (1.6 MB at
  N = 10⁴, k = 20).
- **Counter**: `basis_cache_hits_` (the `rom_diag.csv` column H3 reserved).
- **Tests**: a two-state pumped chain toggling ≥ 3 cycles (second and later
  transitions hit); quantile continuity across a swap; LRU eviction with a
  5-state fixture; the memory bound; and a hit compared with a cold solve on
  the same key (eigenpairs equal to solver tolerance).

## 7. The §7.1 forward question

> Does the cache key and eviction design admit member-indexed lookup later?

**Yes, by keeping two properties.** (1) Entries are value objects keyed only by
control state, never by "the" current basis, so several entries can be live at
once. (2) The advance path already groups members by operator: H14b's
`BatchPropagator` forms the shared powers once per operator and solves per
member. With a per-member state index, members sharing a key share one cache
entry and one propagator, so two branches cost two propagators per step, not
`M`. Nothing in §6 needs to change for that; only the advance loop would gain a
grouping step.

The residual §7.1 already named is unchanged: the deviation form anchors every
member on one deterministic trajectory `h_det`, which has a single control
state. Per-branch anchors are a separate formulation question.

## 8. Reproduction

The probe is a scratch program compiled against the engine dylib (the W3/H4
calibration pattern, not a ctest target): append `[UNCERTAINTY]\n1D MANNINGS_N
0.20` to `tests/regression/data/Example8.inp`, step the engine, and log pump-link
flow sign changes and `rom1d()->basis_rebuilds_cold_forced_` changes per step.
