/**
 * @file test_soft_rain_grid_2d.cpp
 * @brief Focused SR-2c tests for deterministic 2D gridded rainfall forcing.
 *
 * @details Exercises the deterministic /location path (SR-2c) and, since PR
 *          SP2, the /spread -> soft-forcing path: the SR-3c lognormal-CV
 *          warning (restored from the pre-port file), the two-plane family
 *          split, CV/SD conversion, and that spread reaches a live 2D ROM.
 *          The CL-1c correlated-coherence case ("InitRecordsCoherenceCorrLen")
 *          is NOT here: correlated coherence is PR SP3.
 *
 * @ingroup engine_2d
 */

#include <gtest/gtest.h>

// updateRainfall() is private (it's an internal step of the co-advance
// cycle, not public API); reach it directly here rather than drive the
// whole engine lifecycle just to exercise one function.
#define private public
#include "2d/SurfaceRouter2D.hpp"
#undef private

#include "2d/mesh/MeshBuilder.hpp"
#include "core/SimulationContext.hpp"
#include "core/SWMMEngine.hpp"
#include "uncertainty/UncertaintyConfig.hpp"

#include <openswmm/engine/openswmm_engine.h>
#include <openswmm/engine/openswmm_2d.h>

#include <hdf5.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using openswmm::FlowUnits;
using openswmm::SimulationContext;
using openswmm::twoD::MeshData;
using openswmm::twoD::SurfaceRouter2D;
using openswmm::uncertainty::GridMapping;
using openswmm::uncertainty::GridTarget;
using openswmm::uncertainty::SoftGridSourceSpec;

namespace {

static MeshData makeUnitSquareMesh() {
    MeshData mesh;
    mesh.resize_vertices(4);
    mesh.vx = {0.0, 1.0, 0.0, 1.0};
    mesh.vy = {0.0, 0.0, 1.0, 1.0};
    mesh.vz = {0.0, 0.0, 0.0, 0.0};

    mesh.resize_triangles(2);
    mesh.tri_v0[0] = 0; mesh.tri_v1[0] = 1; mesh.tri_v2[0] = 3; // centroid (2/3,1/3)
    mesh.tri_v0[1] = 0; mesh.tri_v1[1] = 3; mesh.tri_v2[1] = 2; // centroid (1/3,2/3)
    mesh.mannings_n[0] = 0.035;
    mesh.mannings_n[1] = 0.035;

    buildMeshTopology(mesh);
    return mesh;
}

void writeStringAttr(hid_t loc, const char* name, const char* value) {
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, std::strlen(value));
    H5Tset_strpad(atype, H5T_STR_NULLTERM);
    hid_t aspace = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(loc, name, atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, value);
    H5Aclose(attr);
    H5Sclose(aspace);
    H5Tclose(atype);
}

std::string makeGridFile(const char* family = "NORMAL",
                         const char* spread_kind = "SD",
                         float s0 = 0.0f, float s1 = 0.0f,
                         float s2 = 0.0f, float s3 = 0.0f) {
    std::string path = "/tmp/test_soft_rain_grid_2d.h5";
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    writeStringAttr(file, "family", family);
    writeStringAttr(file, "spread_kind", spread_kind);
    writeStringAttr(file, "units", "mm/hr");

    // time = [0]
    {
        double time_data[1] = {0.0};
        hsize_t dims[1] = {1};
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t ds = H5Dcreate2(file, "/time", H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_data);
        H5Dclose(ds);
        H5Sclose(space);
    }

    // x = [0,1], y = [0,1]
    {
        double x_data[2] = {0.0, 1.0};
        hsize_t dims[1] = {2};
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t ds = H5Dcreate2(file, "/x", H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_data);
        H5Dclose(ds);
        H5Sclose(space);
    }
    {
        double y_data[2] = {0.0, 1.0};
        hsize_t dims[1] = {2};
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t ds = H5Dcreate2(file, "/y", H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_data);
        H5Dclose(ds);
        H5Sclose(space);
    }

    // spread dummy [1,2,2]
    {
        float spread[4] = {s0, s1, s2, s3};
        hsize_t dims[3] = {1, 2, 2};
        hid_t space = H5Screate_simple(3, dims, nullptr);
        hid_t ds = H5Dcreate2(file, "/spread", H5T_NATIVE_FLOAT, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, spread);
        H5Dclose(ds);
        H5Sclose(space);
    }

    // location [1,2,2] row-major: [[10,20],[30,40]] mm/hr
    {
        float loc[4] = {10.0f, 20.0f, 30.0f, 40.0f};
        hsize_t dims[3] = {1, 2, 2};
        hid_t space = H5Screate_simple(3, dims, nullptr);
        hid_t ds = H5Dcreate2(file, "/location", H5T_NATIVE_FLOAT, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, loc);
        H5Dclose(ds);
        H5Sclose(space);
    }

    H5Fclose(file);
    return path;
}


// General HDF5 grid writer for the SP2 tests: n x n pixels on [0, extent],
// optional /location and /family_code planes, single time plane.
struct GridSpec {
    const char* family      = "NORMAL";
    const char* spread_kind = "SD";
    int    n                = 2;
    double extent           = 1.0;
    std::vector<float>   loc;      // empty => no /location dataset
    std::vector<float>   spread;   // n*n, row-major (iy, ix)
    std::vector<uint8_t> fcode;    // empty => no /family_code (n*n otherwise)
};

std::string writeGrid(const std::string& path, const GridSpec& g) {
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    writeStringAttr(file, "family", g.family);
    writeStringAttr(file, "spread_kind", g.spread_kind);
    writeStringAttr(file, "units", "mm/hr");

    auto writeVec = [&](const char* name, const std::vector<double>& v) {
        hsize_t dims[1] = {static_cast<hsize_t>(v.size())};
        hid_t sp = H5Screate_simple(1, dims, nullptr);
        hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
        H5Dclose(ds); H5Sclose(sp);
    };
    writeVec("/time", {0.0});
    std::vector<double> axis(static_cast<std::size_t>(g.n));
    for (int k = 0; k < g.n; ++k)
        axis[static_cast<std::size_t>(k)] = (g.n > 1) ? g.extent * k / (g.n - 1) : 0.0;
    writeVec("/x", axis);
    writeVec("/y", axis);

    auto writePlane = [&](const char* name, const std::vector<float>& v) {
        hsize_t dims[3] = {1, static_cast<hsize_t>(g.n), static_cast<hsize_t>(g.n)};
        hid_t sp = H5Screate_simple(3, dims, nullptr);
        hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_FLOAT, sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
        H5Dclose(ds); H5Sclose(sp);
    };
    writePlane("/spread", g.spread);
    if (!g.loc.empty()) writePlane("/location", g.loc);
    if (!g.fcode.empty()) {
        hsize_t dims[2] = {static_cast<hsize_t>(g.n), static_cast<hsize_t>(g.n)};
        hid_t sp = H5Screate_simple(2, dims, nullptr);
        hid_t ds = H5Dcreate2(file, "/family_code", H5T_NATIVE_UINT8, sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, g.fcode.data());
        H5Dclose(ds); H5Sclose(sp);
    }
    H5Fclose(file);
    return path;
}

} // anonymous namespace

TEST(SoftRainGrid2D, NoGridSourceUsesLegacyGageFallback) {
    SimulationContext ctx;
    ctx.options.flow_units = FlowUnits::CMS; // SI => mm/hr -> m/s
    ctx.gage_names.add("RG1");
    ctx.gages.resize(1);
    ctx.gages.rainfall[0] = 12.0; // mm/hr

    SurfaceRouter2D router;
    router.mesh() = makeUnitSquareMesh();
    router.state().resize(router.mesh().n_triangles(), router.mesh().n_vertices());

    router.updateRainfall(ctx);

    const double expected = 12.0 * 0.001 / 3600.0;
    ASSERT_EQ(router.state().rainfall.size(), 2u);
    EXPECT_DOUBLE_EQ(router.state().rainfall[0], expected);
    EXPECT_DOUBLE_EQ(router.state().rainfall[1], expected);
}

TEST(SoftRainGrid2D, ForceLocationGridOverridesGageFallback) {
    SimulationContext ctx;
    ctx.options.flow_units = FlowUnits::CMS; // SI => mm/hr -> m/s
    ctx.gage_names.add("RG1");
    ctx.gages.resize(1);
    ctx.gages.rainfall[0] = 999.0; // should be ignored when grid forcing active

    SurfaceRouter2D router;
    router.mesh() = makeUnitSquareMesh();
    router.state().resize(router.mesh().n_triangles(), router.mesh().n_vertices());

    const std::string grid_path = makeGridFile();
    SoftGridSourceSpec spec;
    spec.file_path = grid_path;
    spec.target = GridTarget::TWO_D;
    spec.mapping = GridMapping::CENTROID;
    spec.force_location = true;

    ASSERT_TRUE(router.initGridRainfall(spec, ""));
    ASSERT_TRUE(router.gridRainfallActive());

    router.updateRainfall(ctx);

    // T0 centroid (2/3,1/3) -> nearest cell center (1,0) -> flat index 1 -> 20 mm/hr
    // T1 centroid (1/3,2/3) -> nearest cell center (0,1) -> flat index 2 -> 30 mm/hr
    const double expected_t0 = 20.0 * 0.001 / 3600.0;
    const double expected_t1 = 30.0 * 0.001 / 3600.0;
    ASSERT_EQ(router.state().rainfall.size(), 2u);
    EXPECT_DOUBLE_EQ(router.state().rainfall[0], expected_t0);
    EXPECT_DOUBLE_EQ(router.state().rainfall[1], expected_t1);

    std::remove(grid_path.c_str());
}

TEST(SoftRainGrid2D, LognormalHighCvWarnsOnce) {
    SimulationContext ctx;
    ctx.options.flow_units = FlowUnits::CMS;
    ctx.gage_names.add("RG1");
    ctx.gages.resize(1);
    ctx.gages.rainfall[0] = 0.0;

    SurfaceRouter2D router;
    router.mesh() = makeUnitSquareMesh();
    router.state().resize(router.mesh().n_triangles(), router.mesh().n_vertices());

    // With loc values 20/30 mm/hr after centroid mapping, spread values of 20/30
    // imply CV = 1.0 on both active cells, exceeding the 0.5 SR-3c guard.
    const std::string grid_path = makeGridFile("LOGNORMAL", "CV", 0.0f, 20.0f, 30.0f, 0.0f);
    SoftGridSourceSpec spec;
    spec.file_path = grid_path;
    spec.target = GridTarget::TWO_D;
    spec.mapping = GridMapping::CENTROID;
    spec.force_location = true;

    ASSERT_TRUE(router.initGridRainfall(spec, ""));
    router.updateRainfall(ctx);
    ASSERT_EQ(ctx.warnings.size(), 1u);
    EXPECT_NE(ctx.warnings[0].find("CV > 0.5"), std::string::npos);

    // Second call should not emit an additional warning.
    router.updateRainfall(ctx);
    EXPECT_EQ(ctx.warnings.size(), 1u);

    std::remove(grid_path.c_str());
}

// ============================================================================
// SP2: /spread -> two coefficient planes (router level, no ROM needed)
// ============================================================================
//
// The unit-square mesh maps T0 -> flat pixel 1 and T1 -> flat pixel 2 under
// CENTROID (see ForceLocationGridOverridesGageFallback above).

namespace {

constexpr double kToMs = 0.001 / 3600.0;   // SI mm/hr -> m/s

struct RouterFixture {
    SimulationContext ctx;
    SurfaceRouter2D   router;
    RouterFixture(double gage_mm_hr = 0.0) {
        ctx.options.flow_units = FlowUnits::CMS;
        ctx.gage_names.add("RG1");
        ctx.gages.resize(1);
        ctx.gages.rainfall[0] = gage_mm_hr;
        router.mesh() = makeUnitSquareMesh();
        router.state().resize(router.mesh().n_triangles(), router.mesh().n_vertices());
    }
    bool init(const std::string& path, bool force_location,
              openswmm::uncertainty::Coherence coh = openswmm::uncertainty::Coherence::FULL) {
        SoftGridSourceSpec spec;
        spec.file_path = path;
        spec.target = GridTarget::TWO_D;
        spec.mapping = GridMapping::CENTROID;
        spec.force_location = force_location;
        spec.coherence = coh;
        spec.corr_len = (coh == openswmm::uncertainty::Coherence::CORR_LEN) ? 100.0 : 0.0;
        return router.initGridRainfall(spec, "");
    }
};

}  // namespace

TEST(SoftRainGrid2D, MixedGridSplitsIntoDisjointPlanesWithNoPrescale) {
    RouterFixture fx;
    GridSpec g; g.family = "MIXED"; g.spread_kind = "SD";
    g.loc = {10, 20, 30, 40};
    g.spread = {0, 4, 6, 0};
    g.fcode  = {0, 0, 2, 0};      // T0 (px1) NORMAL, T1 (px2) UNIFORM
    const std::string path = writeGrid("/tmp/sp2_mixed.h5", g);
    ASSERT_TRUE(fx.init(path, /*force_location=*/true));
    fx.router.updateRainfall(fx.ctx);

    const auto& A = fx.router.gridSoftSpreadPlaneA();
    const auto& B = fx.router.gridSoftSpreadPlaneB();
    ASSERT_EQ(A.size(), 2u);
    ASSERT_EQ(B.size(), 2u);
    EXPECT_DOUBLE_EQ(A[0], 4.0 * kToMs);
    EXPECT_DOUBLE_EQ(B[0], 0.0);
    EXPECT_DOUBLE_EQ(B[1], 6.0 * kToMs)
        << "UNIFORM cell must carry its own spread exactly -- the pre-port 3.0x pre-scale is gone";
    EXPECT_DOUBLE_EQ(A[1], 0.0);
    // Disjoint planes that sum to the original spread field.
    EXPECT_DOUBLE_EQ(A[0] + B[0], 4.0 * kToMs);
    EXPECT_DOUBLE_EQ(A[1] + B[1], 6.0 * kToMs);
    std::remove(path.c_str());
}

TEST(SoftRainGrid2D, SingleFamilyGridUsesExactlyOnePlane) {
    {
        RouterFixture fx;
        GridSpec g; g.family = "NORMAL"; g.spread_kind = "SD";
        g.loc = {10, 20, 30, 40}; g.spread = {0, 4, 6, 0};
        const std::string path = writeGrid("/tmp/sp2_normal.h5", g);
        ASSERT_TRUE(fx.init(path, true));
        fx.router.updateRainfall(fx.ctx);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneA()[0], 4.0 * kToMs);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneA()[1], 6.0 * kToMs);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneB()[0], 0.0);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneB()[1], 0.0);
        std::remove(path.c_str());
    }
    {
        RouterFixture fx;
        GridSpec g; g.family = "UNIFORM"; g.spread_kind = "HALFRANGE";
        g.loc = {10, 20, 30, 40}; g.spread = {0, 4, 6, 0};
        const std::string path = writeGrid("/tmp/sp2_uniform.h5", g);
        ASSERT_TRUE(fx.init(path, true));
        fx.router.updateRainfall(fx.ctx);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneB()[0], 4.0 * kToMs);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneB()[1], 6.0 * kToMs);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneA()[0], 0.0);
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneA()[1], 0.0);
        std::remove(path.c_str());
    }
}

TEST(SoftRainGrid2D, CvSpreadScalesByTheForcedLocationRate) {
    RouterFixture fx;
    GridSpec g; g.family = "NORMAL"; g.spread_kind = "CV";
    g.loc = {10, 20, 30, 40}; g.spread = {0, 0.25f, 0.5f, 0};
    const std::string path = writeGrid("/tmp/sp2_cv_forced.h5", g);
    ASSERT_TRUE(fx.init(path, true));
    fx.router.updateRainfall(fx.ctx);
    EXPECT_NEAR(fx.router.gridSoftSpreadPlaneA()[0], 0.25 * 20.0 * kToMs, 1e-18);
    EXPECT_NEAR(fx.router.gridSoftSpreadPlaneA()[1], 0.50 * 30.0 * kToMs, 1e-18);
    std::remove(path.c_str());
}

TEST(SoftRainGrid2D, SpreadOnlyGridLeavesTheGageRainfallAsTheLocation) {
    // No FORCE_LOCATION: the model's own rainfall stays the deterministic
    // location and the grid supplies only spread. This is what the USER_GUIDE
    // has always documented; before SP2 such a grid was never even opened.
    RouterFixture fx(/*gage_mm_hr=*/12.0);
    GridSpec g; g.family = "NORMAL"; g.spread_kind = "SD";
    g.loc = {10, 20, 30, 40}; g.spread = {0, 4, 6, 0};
    const std::string path = writeGrid("/tmp/sp2_spread_only.h5", g);
    ASSERT_TRUE(fx.init(path, /*force_location=*/false));
    ASSERT_TRUE(fx.router.gridRainfallActive());
    fx.router.updateRainfall(fx.ctx);
    // Deterministic rainfall is the GAGE value, not the grid's /location.
    EXPECT_DOUBLE_EQ(fx.router.state().rainfall[0], 12.0 * kToMs);
    EXPECT_DOUBLE_EQ(fx.router.state().rainfall[1], 12.0 * kToMs);
    // ... while the spread still reaches the plane.
    EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneA()[0], 4.0 * kToMs);
    std::remove(path.c_str());
}

TEST(SoftRainGrid2D, CvSpreadWithNoLocationPlaneUsesTheGageRate) {
    // A grid with no /location dataset. The pre-port block multiplied a CV
    // spread by the (absent) /location value, giving zero spread; the rainfall
    // actually in force (here the gage) is the right reference.
    RouterFixture fx(/*gage_mm_hr=*/12.0);
    GridSpec g; g.family = "NORMAL"; g.spread_kind = "CV";
    g.spread = {0, 0.5f, 0.25f, 0};       // no g.loc
    const std::string path = writeGrid("/tmp/sp2_cv_noloc.h5", g);
    ASSERT_TRUE(fx.init(path, /*force_location=*/false));
    fx.router.updateRainfall(fx.ctx);
    EXPECT_NEAR(fx.router.gridSoftSpreadPlaneA()[0], 0.5  * 12.0 * kToMs, 1e-18);
    EXPECT_NEAR(fx.router.gridSoftSpreadPlaneA()[1], 0.25 * 12.0 * kToMs, 1e-18);
    EXPECT_DOUBLE_EQ(fx.router.state().rainfall[0], 12.0 * kToMs);
    std::remove(path.c_str());
}

TEST(SoftRainGrid2D, LastPlaneIsHeldAfterItsTimestamp) {
    // Regression: GridFileReader::time_next() returns the CURRENT plane's own
    // time when there is no next plane, so once t_now passed the last plane the
    // advance loop ran off the end and the grid silently stopped forcing (a
    // single-plane grid only worked at t = 0). Every earlier test ran at t = 0.
    RouterFixture fx(/*gage_mm_hr=*/999.0);   // must stay ignored
    GridSpec g; g.family = "NORMAL"; g.spread_kind = "SD";
    g.loc = {10, 20, 30, 40}; g.spread = {0, 4, 6, 0};   // single plane at t = 0
    const std::string path = writeGrid("/tmp/sp2_last_plane.h5", g);
    ASSERT_TRUE(fx.init(path, /*force_location=*/true));

    for (double t : {0.0, 30.0, 600.0, 7200.0}) {
        fx.ctx.current_time = t;
        fx.router.updateRainfall(fx.ctx);
        EXPECT_DOUBLE_EQ(fx.router.state().rainfall[0], 20.0 * kToMs)
            << "grid /location lost at t=" << t;
        EXPECT_DOUBLE_EQ(fx.router.gridSoftSpreadPlaneA()[0], 4.0 * kToMs)
            << "grid /spread lost at t=" << t;
    }
    std::remove(path.c_str());
}

// CL-1c (restored from the pre-port file, SP3): initGridRainfall records the
// coherence correlation length; FULL / default leaves it at 0 (comonotone).
TEST(SoftRainGrid2D, InitRecordsCoherenceCorrLen) {
    GridSpec g; g.loc = {10, 20, 30, 40}; g.spread = {0, 4, 6, 0};
    const std::string path = writeGrid("/tmp/sp3_corrlen_init.h5", g);
    {
        RouterFixture fx;
        ASSERT_TRUE(fx.init(path, true));
        EXPECT_DOUBLE_EQ(fx.router.gridSoftCorrLen(), 0.0);
    }
    {
        RouterFixture fx;
        ASSERT_TRUE(fx.init(path, true, openswmm::uncertainty::Coherence::CORR_LEN));
        EXPECT_DOUBLE_EQ(fx.router.gridSoftCorrLen(), 100.0);
        fx.router.updateRainfall(fx.ctx);
        EXPECT_TRUE(fx.ctx.warnings.empty()) << "no 'not wired' warning since SP3: "
                                             << (fx.ctx.warnings.empty() ? "" : fx.ctx.warnings[0]);
    }
    std::remove(path.c_str());
}

TEST(SoftRainGrid2D, CorrLenWithMixedGridIsRefusedAtInit) {
    // H7 scope guard at the first point the file's family is known.
    RouterFixture fx;
    GridSpec g; g.family = "MIXED"; g.loc = {10, 20, 30, 40}; g.spread = {0, 4, 6, 0};
    g.fcode = {0, 0, 2, 0};
    const std::string path = writeGrid("/tmp/sp3_corrlen_mixed.h5", g);
    EXPECT_FALSE(fx.init(path, true, openswmm::uncertainty::Coherence::CORR_LEN));
    EXPECT_FALSE(fx.router.gridRainfallActive());
    EXPECT_NE(fx.router.gridInitError().find("CORR_LEN"), std::string::npos)
        << fx.router.gridInitError();
    // Same file under FULL is fine.
    RouterFixture fy;
    EXPECT_TRUE(fy.init(path, true));
    EXPECT_TRUE(fy.router.gridInitError().empty());
    std::remove(path.c_str());
}

// ============================================================================
// End-to-end: [SOFT_RAINFALL_GRID] parsing -> initHydraulics() trigger ->
// SurfaceRouter2D::initGridRainfall(), through the real engine lifecycle.
// The two tests above call initGridRainfall() directly; this one is the
// guard that the wiring reaching it — registerSoftRainfallGridSection and
// the uncertainty_config_.grid_sources scan in SWMMEngine::initHydraulics()
// — actually connects to a real .inp.
// ============================================================================

namespace {

std::string buildGridModel(const std::string& grid_path) {
    return
        "[OPTIONS]\n"
        "FLOW_UNITS           CMS\n"
        "FLOW_ROUTING         DYNWAVE\n"
        "START_DATE           01/01/2026\n"
        "START_TIME           00:00:00\n"
        "END_DATE             01/01/2026\n"
        "END_TIME             00:05:00\n"
        "REPORT_STEP          00:01:00\n"
        "ROUTING_STEP         6\n"
        "\n"
        "[JUNCTIONS]\n"
        ";;Name  Elev  MaxDepth  InitDepth  SurDepth  Aponded\n"
        "J1      0.0   1.0       0          0         0\n"
        "\n"
        "[OUTFALLS]\n"
        ";;Name  Elev   Type  Gated\n"
        "O1     -0.5    FREE  NO\n"
        "\n"
        "[CONDUITS]\n"
        ";;Name  From  To  Length  Roughness  InOffset  OutOffset  InitFlow\n"
        "C1      J1    O1  30.0    0.013      0         0          0\n"
        "\n"
        "[XSECTIONS]\n"
        ";;Link  Shape     Geom1  Geom2  Geom3  Geom4  Barrels\n"
        "C1      CIRCULAR  0.3    0      0      0      1\n"
        "\n"
        "[2D_VERTICES]\n"
        ";;X    Y    Z\n"
        " 0.0   0.0  0.0\n"
        " 1.0   0.0  0.0\n"
        " 0.0   1.0  0.0\n"
        " 1.0   1.0  0.0\n"
        "\n"
        "[2D_TRIANGLES]\n"
        ";;V1  V2  V3  MANNINGS_N\n"
        "0     1   3   0.035\n"
        "0     3   2   0.035\n"
        "\n"
        "[2D_VERTEX_NODE_MAP]\n;;Vertex  Node  Cd   Area\n0  J1  0.7  1.0\n"
        "\n[SOFT_RAINFALL_GRID]\n"
        ";;Target File Mapping [Options]\n"
        "2D  \"" + grid_path + "\"  CENTROID  FORCE_LOCATION\n";
}

}  // namespace

TEST(SoftRainGrid2D, InpSectionReachesInitGridRainfallThroughTheRealEngine) {
    namespace fs = std::filesystem;
    const fs::path dir = fs::current_path() / "soft_rain_grid_2d_out";
    fs::create_directories(dir);

    const std::string grid_path = makeGridFile();
    const fs::path inp = dir / "grid_end_to_end.inp";
    const fs::path rpt = dir / "grid_end_to_end.rpt";
    const fs::path out = dir / "grid_end_to_end.out";
    { std::ofstream f(inp); f << buildGridModel(grid_path); }

    SWMM_Engine eng = swmm_engine_create();
    ASSERT_EQ(swmm_engine_open(eng, inp.string().c_str(), rpt.string().c_str(),
                               out.string().c_str(), nullptr), SWMM_OK);
    ASSERT_EQ(swmm_engine_initialize(eng), SWMM_OK);

    int active = 0;
    swmm_2d_is_active(eng, &active);
    ASSERT_TRUE(active);

    auto* impl = static_cast<openswmm::SWMMEngine*>(eng);
    EXPECT_TRUE(impl->surfaceRouter2D().gridRainfallActive())
        << "[SOFT_RAINFALL_GRID] FORCE_LOCATION did not reach initGridRainfall()";

    swmm_engine_close(eng);
    swmm_engine_destroy(eng);
    std::remove(grid_path.c_str());
}

// ============================================================================
// SP2: the planes reach a LIVE 2D ROM through the real engine
// ============================================================================
//
// 4-triangle sloped domain (the Phase-8 ROM fixture). Manning and rainfall
// perturbation are set to ZERO, so the only possible source of ensemble spread
// is the grid's /spread plane -- the no-grid baseline below proves it is zero
// otherwise, which makes a nonzero band attributable to SP2.

namespace {

std::string buildRomGridModel(const std::string& grid_path, bool force_location) {
    std::string m =
        "[OPTIONS]\n"
        "FLOW_UNITS           CMS\n"
        "FLOW_ROUTING         DYNWAVE\n"
        "START_DATE           01/01/2025\n"
        "START_TIME           00:00:00\n"
        "REPORT_START_DATE    01/01/2025\n"
        "REPORT_START_TIME    00:00:00\n"
        "END_DATE             01/01/2025\n"
        "END_TIME             00:05:00\n"
        "REPORT_STEP          00:01:00\n"
        "ROUTING_STEP         0:00:30\n"
        "MIN_SURFAREA         1.0\n"
        "MINIMUM_STEP         0.5\n"
        "THREADS              1\n\n"
        "[JUNCTIONS]\nJ1 1.0 10.0 0.0 0.0 0\n\n"
        "[OUTFALLS]\nO1 0.0 FREE NO\n\n"
        "[CONDUITS]\nC1 J1 O1 100.0 0.035 0.0 0.0\n\n"
        "[XSECTIONS]\nC1 CIRCULAR 1.0 0.0 0.0 0.0 1\n\n"
        "[2D_OPTIONS]\nMAX_TIMESTEP 30.0\nMIN_TIMESTEP 0.01\nDRY_DEPTH 0.001\n\n"
        "[2D_VERTICES]\n"
        "0.0 0.0 0.20\n3.0 0.0 0.15\n3.0 3.0 0.10\n0.0 3.0 0.05\n1.5 1.5 0.12\n\n"
        "[2D_TRIANGLES]\n0 1 4 0.035\n1 2 4 0.035\n2 3 4 0.035\n3 0 4 0.035\n\n"
        "[2D_ROM]\nENABLE YES\nMEMBERS 20\nMODES 4\nMANNINGS_PERT 0.0\nRAINFALL_PERT 0.0\n\n";
    if (!grid_path.empty())
        m += "[SOFT_RAINFALL_GRID]\n2D \"" + grid_path + "\" CENTROID"
             + (force_location ? "  FORCE_LOCATION" : "") + "\n";
    return m;
}

struct RomProbe {
    bool   ok = false;
    bool   grid_active = false;
    bool   has_coeff_a = false, has_coeff_b = false;
    std::string planes_error;
    double band = 0.0;           // max over cells of q95 - q05
    int    corr_len_warnings = 0;
    int    mixed_family_warnings = 0;
};

RomProbe runRomModel(const std::string& tag, const std::string& grid_path,
                     bool force_location) {
    namespace fs = std::filesystem;
    RomProbe pr;
    const fs::path dir = fs::current_path() / "soft_rain_grid_2d_out";
    fs::create_directories(dir);
    const fs::path inp = dir / (tag + ".inp");
    const fs::path rpt = dir / (tag + ".rpt");
    const fs::path out = dir / (tag + ".out");
    { std::ofstream f(inp); f << buildRomGridModel(grid_path, force_location); }

    SWMM_Engine eng = swmm_engine_create();
    if (!eng) return pr;
    if (swmm_engine_open(eng, inp.string().c_str(), rpt.string().c_str(),
                         out.string().c_str(), nullptr) != SWMM_OK ||
        swmm_engine_initialize(eng) != SWMM_OK ||
        swmm_engine_start(eng, 0) != SWMM_OK) {
        swmm_engine_destroy(eng);
        return pr;
    }
    double elapsed = 1.0;
    for (int n = 0; elapsed > 0.0 && n < 2000; ++n)
        if (swmm_engine_step(eng, &elapsed) != SWMM_OK) break;

    auto* impl = static_cast<openswmm::SWMMEngine*>(eng);
    pr.grid_active = impl->surfaceRouter2D().gridRainfallActive();
    if (const auto* rom = impl->surfaceRouter2D().rom()) {
        if (rom->is_ready()) {
            pr.ok = true;
            pr.has_coeff_a = !rom->softCoeff().empty();
            pr.has_coeff_b = !rom->softCoeffB().empty();
            pr.planes_error = rom->softPlanesError();
            for (std::size_t i = 0; i < rom->q95.size(); ++i)
                pr.band = std::max(pr.band, rom->q95[i] - rom->q05[i]);
        }
    }
    for (const auto& w : impl->context().warnings) {
        if (w.find("CORR_LEN") != std::string::npos) ++pr.corr_len_warnings;
        if (w.find("first") != std::string::npos &&
            w.find("famil") != std::string::npos) ++pr.mixed_family_warnings;
    }
    swmm_engine_end(eng);
    swmm_engine_close(eng);
    swmm_engine_destroy(eng);
    return pr;
}

// 3x3 grid on [0,3]: the four triangle centroids hit four DISTINCT pixels
// (flat 1, 3, 5, 7), so a non-uniform spread survives the zero-mean modes.
GridSpec romGrid(const char* family, const char* kind) {
    GridSpec g; g.family = family; g.spread_kind = kind; g.n = 3; g.extent = 3.0;
    g.loc    = {0, 30, 0,  20, 0, 10,  0, 40, 0};
    g.spread = {0,  9, 0,   3, 0,  6,  0, 12, 0};
    return g;
}

}  // namespace

TEST(SoftRainGrid2DRom, WithoutAGridTheBandIsExactlyZero) {
    // Baseline: perturbations are zero, so nothing but the grid can widen it.
    const RomProbe pr = runRomModel("rom_baseline", "", false);
    ASSERT_TRUE(pr.ok) << "ROM did not come up";
    // (The ROM fills its coefficient column at seed time whether or not a grid
    // is present, so has_coeff_a says nothing here; the band is the evidence.)
    EXPECT_EQ(pr.band, 0.0);
}

TEST(SoftRainGrid2DRom, SpreadOnlyGridWithoutForceLocationDrivesTheRom) {
    const std::string path = writeGrid("/tmp/sp2_rom_normal.h5", romGrid("NORMAL", "SD"));
    const RomProbe pr = runRomModel("rom_spread_only", path, /*force_location=*/false);
    ASSERT_TRUE(pr.ok);
    EXPECT_TRUE(pr.grid_active) << "a TWO_D grid without FORCE_LOCATION must still open";
    EXPECT_TRUE(pr.has_coeff_a);
    EXPECT_FALSE(pr.has_coeff_b) << "single-family grid must leave the second plane unarmed";
    EXPECT_GT(pr.band, 1.0e-9) << "grid /spread never reached the ROM";
    std::remove(path.c_str());
}

TEST(SoftRainGrid2DRom, ZeroSpreadGridCollapsesTheBandExactly) {
    GridSpec g = romGrid("NORMAL", "SD");
    std::fill(g.spread.begin(), g.spread.end(), 0.0f);
    const std::string path = writeGrid("/tmp/sp2_rom_zero.h5", g);
    const RomProbe pr = runRomModel("rom_zero", path, true);
    ASSERT_TRUE(pr.ok);
    EXPECT_EQ(pr.band, 0.0) << "zero spread must collapse the band exactly (deviation form)";
    std::remove(path.c_str());
}

TEST(SoftRainGrid2DRom, MixedGridArmsBothPlanesWithNoFallback) {
    GridSpec g = romGrid("MIXED", "SD");
    // Pixels 1,3 (T0,T3) NORMAL; pixels 5,7 (T1,T2) UNIFORM.
    g.fcode = {0, 0, 0,  0, 0, 2,  0, 2, 0};
    const std::string path = writeGrid("/tmp/sp2_rom_mixed.h5", g);
    const RomProbe pr = runRomModel("rom_mixed", path, true);
    ASSERT_TRUE(pr.ok);
    EXPECT_TRUE(pr.has_coeff_a);
    EXPECT_TRUE(pr.has_coeff_b) << "UNIFORM cells never reached the second coefficient plane";
    EXPECT_TRUE(pr.planes_error.empty()) << pr.planes_error;
    EXPECT_EQ(pr.mixed_family_warnings, 0);
    EXPECT_GT(pr.band, 1.0e-9);
    std::remove(path.c_str());
}

TEST(SoftRainGrid2DRom, CorrLenGridTakesTheReducedPathAndChangesTheBand) {
    namespace fs = std::filesystem;
    // Same NORMAL grid under FULL and under CORR_LEN <<< centroid spacing
    // (~1.5 m): the correlated path must engage and give a band that differs
    // from comonotone. On 4 points the SPDE basis retains K_s >= M = 20 modes,
    // so this takes the MATERIALIZED branch, not the reduced one -- both are
    // the correlated path. Direction is MEASURED, not asserted: here the band
    // WIDENS (2.02x, 2026-10-03). Comonotone forcing on a tiny flat patch is
    // mostly a uniform shift, which lives in the discarded constant mode;
    // per-cell independent coefficients project more onto the zero-mean modes.
    // The 1D chain narrows downstream (CL-1e); the two are not in conflict.
    const std::string path = writeGrid("/tmp/sp3_rom_cl.h5", romGrid("NORMAL", "SD"));
    const RomProbe full = runRomModel("rom_cl_full", path, true);
    ASSERT_TRUE(full.ok);
    // Rewrite the .inp with COHERENCE CORR_LEN 0.1 appended to the grid line.
    const fs::path dir = fs::current_path() / "soft_rain_grid_2d_out";
    const fs::path inp = dir / "rom_cl_corr.inp";
    { std::ofstream f(inp);
      std::string m = buildRomGridModel(path, true);
      m.replace(m.rfind("\n"), 1, "  COHERENCE CORR_LEN 0.1\n");
      f << m; }
    SWMM_Engine eng = swmm_engine_create();
    ASSERT_EQ(swmm_engine_open(eng, inp.string().c_str(), (dir / "rom_cl_corr.rpt").string().c_str(),
                               (dir / "rom_cl_corr.out").string().c_str(), nullptr), SWMM_OK)
        << swmm_get_last_error_msg(eng);
    ASSERT_EQ(swmm_engine_initialize(eng), SWMM_OK) << swmm_get_last_error_msg(eng);
    ASSERT_EQ(swmm_engine_start(eng, 0), SWMM_OK);
    double elapsed = 1.0;
    for (int n = 0; elapsed > 0.0 && n < 2000; ++n)
        if (swmm_engine_step(eng, &elapsed) != SWMM_OK) break;
    auto* impl = static_cast<openswmm::SWMMEngine*>(eng);
    const auto& router = impl->surfaceRouter2D();
    EXPECT_DOUBLE_EQ(router.gridSoftCorrLen(), 0.1);
    EXPECT_TRUE(router.gridSoftCorrelatedActive())
        << "correlated path never engaged (neither reduced nor materialized)";
    const int ks = router.gridSoftBasisModes();
    const bool reduced = router.gridSoftReduced();
    double band_corr = 0.0;
    const auto* rom = router.rom();
    ASSERT_NE(rom, nullptr);
    for (std::size_t i = 0; i < rom->q95.size(); ++i)
        band_corr = std::max(band_corr, rom->q95[i] - rom->q05[i]);
    int fallback_warnings = 0;
    for (const auto& w : impl->context().warnings)
        if (w.find("CORR_LEN") != std::string::npos) ++fallback_warnings;
    swmm_engine_end(eng); swmm_engine_close(eng); swmm_engine_destroy(eng);

    std::printf("[SP3 2D] band FULL=%.6e  CORR_LEN(0.1 m)=%.6e  ratio=%.3f  K_s=%d reduced=%d\n",
                full.band, band_corr, band_corr / full.band, ks, reduced ? 1 : 0);
    EXPECT_EQ(fallback_warnings, 0);
    EXPECT_GT(band_corr, 1.0e-9);
    EXPECT_NE(band_corr, full.band) << "CORR_LEN had no effect on the band";
    std::remove(path.c_str());
}

TEST(SoftRainGrid2DRom, UniformOnlyGridUsesThePrimaryPlaneAlone) {
    const std::string path = writeGrid("/tmp/sp2_rom_uniform.h5", romGrid("UNIFORM", "HALFRANGE"));
    const RomProbe pr = runRomModel("rom_uniform", path, true);
    ASSERT_TRUE(pr.ok);
    EXPECT_TRUE(pr.has_coeff_a);
    EXPECT_FALSE(pr.has_coeff_b);
    EXPECT_GT(pr.band, 1.0e-9);
    std::remove(path.c_str());
}

