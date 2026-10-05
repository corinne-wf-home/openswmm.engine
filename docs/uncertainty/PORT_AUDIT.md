# Port Audit: `swmm6_rel` integration

**Date**: 2026-07-20
**Target Branch**: `swmm6_rel`
**Local Lineage**: `feature/uncertainty-sidecar` through `softrain/cl2d-docs` (62 commits since merge-base `58c2f82c`)
**Upstream Progress**: 216 commits since merge-base

This document analyzes the porting surface to integrate the uncertainty sidecar (H_sym, reform PRs, post-reform, and soft rainfall CL-1/CL-2) into the active `swmm6_rel` production line.

## A. Purely-Additive Local Files (Expected Clean)

These files do not exist upstream and will port cleanly. There are 89 such files, primarily entirely new modules for uncertainty quantification.

**Key directories & files:**
- `docs/uncertainty/*` (USER_GUIDE, HOW_IT_WORKS, VALIDATION, etc.)
- `src/engine/uncertainty/*` (GraphEigenBasis, SpectralROM1D, LhsShuffle, SpdeSpatialBasis, etc.)
- `src/engine/2d/uncertainty/*` (SpectralROM, CorrelatedFieldGenerator, etc.)
- `src/engine/hydraulics/SpectralCoarse.*`
- `src/engine/2d/solver/SpectralPrecond2D.*`
- `tests/regression/test_*` (ROM coverage, profiling, etc.)
- `tests/benchmarks/bench_2d_spectral*`, `demo_1d_spectral_rom*`

## B. Genuine Conflict Surfaces (Intersection)

These are files modified both locally and upstream since the merge base.
They require careful reconciliation.

**1. CMake Configurations**
- `src/engine/CMakeLists.txt`
- `tests/unit/engine/CMakeLists.txt`
*(Additions of the `uncertainty/` sources and test executables need to be merged into the current upstream lists.)*

**2. SWMMEngine & Core Lifecycle Hooks**
- `src/engine/core/SWMMEngine.cpp` & `.hpp`
- `src/engine/core/SimulationContext.hpp`
- `src/engine/core/SimulationOptions.hpp`
*(ROM initialization, memory allocation, `buildRom1DSoftField`, and `execute()` timestep injection hooks. Note `NODE_CONTINUITY SEMI_IMPLICIT` was added to `SimulationOptions` upstream.)*

**3. 2D Solver Hooks & Refactoring**
- `src/engine/2d/solver/CvodeSurfaceSolver.cpp` & `.hpp`
- `src/engine/2d/SurfaceRouter2D.cpp` & `.hpp`
- `src/engine/2d/input/SectionHandlers2D.cpp` & `.hpp`
- `src/engine/2d/coupling/NodeCoupling.hpp`
*(Upstream introduced `ISurfaceSolver` and `ArkodeSurfaceSolver`. `CvodeSurfaceSolver` contains our sidecar hooks, which must be hoisted to `SurfaceRouter2D` or defined cleanly to support both CVODE and ARKODE.)*

**4. Input Initialization & Parsing**
- `src/engine/input/handlers/OptionsHandler.cpp`
- `src/engine/input/PostParseResolver.cpp`

**5. `links.*` Readers (Phase 6 Data-Layout Refactor)**
*Note: The upstream Phase 6 refactor deleted wide `LinkData` arrays in favor of side-tables.*
Our local sidecar contains many readers that will break compilation or silently fail:
- `computeK1d`
- `NetworkLaplacian1D`
- Test fixtures touching DW/KW links.
*These must be manually migrated to the side-table pattern in P2.*

## C. Upstream Fixes & Parity Constraints

Local engine fixes or adjustments made during the sidecar development must be reconciled with upstream fixes:
- **Anderson Acceleration**: Upstream added proper convergence and history mixing logic (PRs #97, #98, #100, #101). The ROM's state snapshot MUST occur post-`execute()` to read the canonical accepted node depths.
- **DW wet/dry wetting criterion**: We experienced final-boundary stalls and relaxed tolerances on `dw-ritter-drybed-strip`. Upstream maintains legacy float32 bit-parity. We must ensure our integration does not alter the solver output unless specifically option-gated (like `[OPTIONS] DW_WETTING THIN_FILM`).
- **OADate / Final-Boundary Stall Fixes**: Do not blindly re-apply local hotfixes if upstream has a canonical solution that maintains the parity ladder.

## D. Recommended Porting Strategy

**Strategy**: Topically re-apply the sidecar as a fresh stacked series, rather than doing a 62-to-216-commit rebase.

**Why?**
1. Re-applying avoids 216-commit-deep conflict archaeology on structural changes.
2. It allows us to naturally adapt the `links.*` (Phase 6) and 2D solver (`ISurfaceSolver`) architectural shifts into the new commits.
3. The majority of the sidecar (89 files) is purely additive. The hooks in SWMMEngine and the CMake files can be re-established cleanly.

**Owner Sign-off**: ✅ Target `swmm6_rel`. Proceed via fresh application on the R0 porting branches.

## E. Execution Plan for Wave R0 PR P2 (Mechanical Port)
*Instructions for Claude Code execution of the P2 sidecar port.*

**Target Branch setup:**
```bash
git checkout swmm6_rel
git checkout -b port/p2-sidecar-on-rel
```

**Step-by-step chunks:**

1. **Chunk 1: Additive Files & Build System**
   - Copy over the 89 purely additive files (from `docs/uncertainty`, `src/engine/uncertainty`, `src/engine/2d/uncertainty`, benchmarks, and unit tests).
   - Integrate the `CMakeLists.txt` additions: add the new source files and test executable targets to upstream's `src/engine/CMakeLists.txt` and `tests/unit/engine/CMakeLists.txt`.

2. **Chunk 2: Phase 6 Data-Layout Migration (Crucial)**
   - Find all `links.` and `nodes.` readers in the sidecar code (e.g., `computeK1d` in `SpectralROM1D.cpp`, `NetworkLaplacian1D`, `NodeCoupling`).
   - Upstream removed wide property arrays in Phase 6. Migrate these to use the new side-table data structures (follow the compiler errors, but pay extreme attention to `mod_length` authority and unit conversions).
   - *Validation*: K1d values on test fixtures must exactly match pre-port values (e.g., to 1e-12).

3. **Chunk 3: SWMMEngine & Parser Hooks**
   - Re-insert the initialization, advance, and output hooks into `SWMMEngine.cpp` / `.hpp` and `SimulationContext.hpp`.
   - Port the parser hooks: verify if the sections should still sit in `SectionHandlers2D.cpp` or utilize the new `input/handlers/` pattern.

4. **Chunk 4: 2D Solver Hooks**
   - Restore the sidecar integration points within `CvodeSurfaceSolver.cpp` (and `SurfaceRouter2D`).
   - *Note*: PR P3 will deal with properly re-homing these into the `ISurfaceSolver` architecture to support both CVODE and ARKODE. For P2, the goal is just getting it to compile, wire up, and pass the tests under CVODE.

5. **Chunk 5: Test Fixtures and Parity Checks**
   - Restore modification to test fixtures taking care not to overwrite upstream's new tests.
   - Run the full test gate (`ctest`).
   - **Escalate (to O48/Fable)** if there are *any* deviations in spread magnitude or statistical test failures compared to the `softrain/cl2d-docs` baseline. The port must not move statistics.