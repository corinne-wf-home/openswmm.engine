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

Cross-engine result at adoption (abs 0.06 / rel 25 %, 5 % fraction gate), all
runs bit-identical across two runs:

| Example | periods | checks | exceeding | fraction | worst single point |
|---|---|---|---|---|---|
| 1 | 360 | 8,280 | 0 | 0.000 % | none |
| 2 | 216 | 4,968 | 48 | 0.966 % | link 9 at t=0.49 d: legacy 1.11, new 64.6 |
| 3 | 216 | 4,968 | 0 | 0.000 % | none |
| 6 | 2,880 | 118,080 | 1,315 | 1.114 % | link 7 at t=0.04 d: legacy 3e-6, new 0.229 |
| 7 | 720 | 28,080 | 0 | 0.000 % | none |
| 8 | 720 | 47,520 | 0 | 0.000 % | none |

The two outliers (Example2 link 9, Example6 link 7) are short transients and
sit well inside the gate; they are recorded here so a change in them is
noticed, not accepted as understood.
