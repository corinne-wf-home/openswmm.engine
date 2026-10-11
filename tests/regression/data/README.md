# Regression Test Data

EPA SWMM5 Example networks used by `test_engine_regression` (both suites:
`Examples/RegressionTest.*` cross-engine vs the legacy engine, and
`Examples/SelfConsistencyTest.*` bit-identical determinism) and by
`bench_engine_vs_legacy`. All six are public domain (EPA SWMM 5 distribution /
Applications Manual) and are committed here; none exceeds 200 KB, so no LFS.

| File | Source | Features | Run |
|---|---|---|---|
| `Example1.inp` | SWMM 5 distribution | small drainage, DW routing | 6 h |
| `Example2.inp` | SWMM 5 distribution | LID, combined sewer | 6 h |
| `Example3.inp` | SWMM 5 distribution | water-quality routing | 6 h |
| `Example6.inp` | Applications Manual, Example6-Final (PR H8, 2026-10-10) | runoff treatment: storage, orifices, weirs, pollutants | 2 d |
| `Example7.inp` | Applications Manual, Example7-Final (PR H8) | dual drainage, irregular (TRANSECT) channels | 12 h |
| `Example8.inp` | Applications Manual, Example8 (PR H8) | combined sewer: pumps, storage, orifices, weirs, DWF | 12 h |

Examples 6-8 were fetched from the `paleolimbot/swmm` R package
(`inst/swmm_examples/`), which redistributes the EPA files unchanged; the EPA
originals ship with the SWMM 5 installer
(https://www.epa.gov/water-research/storm-water-management-model-swmm).

Cross-engine result (2026-10-11, after the rebase onto swmm6_rel): every
model in the suite, Examples 6-8 included, matches the legacy engine at the
suite's strict tolerance (0.001 abs / 0.1 % rel) at EVERY routing step, and
every model is bit-identical across two runs of the new engine.

History: H8 first registered Examples 6-8 on a line 725 commits behind
swmm6_rel, where the cross-engine engine differences forced the older Phase 7
"5 % of comparisons may exceed 0.06 / 25 %" gate (Example 2 and Example 6 then
showed 0.97 % and 1.11 % exceedances). Upstream's parity work closed those
differences; the rebase kept upstream's strict suite and added the three
networks to it unchanged.
