#include "kairoq/ipc2581.hpp"
#include "kairoq/mesh.hpp"
#include "kairoq/model_io.hpp"
#include "kairoq/open.hpp"

#include <doctest/doctest.h>

#include <filesystem>

using namespace kairoq;
namespace fs = std::filesystem;

namespace {

const std::string kData = KAIROQ_TESTDATA;

Board small_board() {
    Board b;
    b.name = "small";
    b.nets = {"", "GND", "VCC"};
    b.net_source = "x2";
    b.stackup = {{"F.Cu", LayerType::Copper, 35000, "copper", {}, {}, 5.8e7, 1},
                 {"core", LayerType::Dielectric, 1.5e6, "FR-4", 4.5, 0.02, {}, 0},
                 {"B.Cu", LayerType::Copper, 35000, "copper", {}, {}, {}, 2}};
    const Path64 sq{Point64(0, 0), Point64(10000000, 0), Point64(10000000, 10000000), Point64(0, 10000000)};
    const Path64 hole{Point64(4000000, 4000000), Point64(4000000, 6000000), Point64(6000000, 6000000), Point64(6000000, 4000000)};
    b.outline = {Polygon{sq, {}}};
    b.copper = {{1, "F.Cu", Side::Top, {{1, {Polygon{sq, {hole}}}}, {2, {Polygon{hole, {}}}}}}, {2, "B.Cu", Side::Bottom, {}}};
    std::reverse(b.copper[0].nets[1].polygons[0].outer.begin(), b.copper[0].nets[1].polygons[0].outer.end());
    b.masks = {{Side::Top, {Polygon{hole, {}}}}};
    std::reverse(b.masks[0].openings[0].outer.begin(), b.masks[0].openings[0].outer.end());
    b.pads = {{"R1", "1", "SMDPad", 2, 1, Point64(5000000, 5000000), {Polygon{sq, {}}}}};
    b.drills = {{Point64(1000000, 1000000), std::nullopt, 300000, true, DrillFunction::Via, 1, 2, 1},
                {Point64(2000000, 1000000), Point64(3000000, 1000000), 800000, false, DrillFunction::Mechanical, 1, 2, 0}};
    b.components = {{"R1", "1k", "R0603", Point64(5000000, 5000000), 90.0, Side::Top, "placement"}};
    b.warnings = {"example"};
    return b;
}

void check_same(const Board& a, const Board& b) {
    CHECK(a.name == b.name);
    CHECK(a.nets == b.nets);
    CHECK(a.stackup.size() == b.stackup.size());
    CHECK(*b.stackup[1].epsilon_r == doctest::Approx(4.5));
    CHECK(*b.stackup[0].conductivity_s_per_m == doctest::Approx(5.8e7));
    REQUIRE(b.copper.size() == 2);
    REQUIRE(b.copper[0].nets.size() == 2);
    CHECK(b.copper[0].nets[0].polygons[0].outer == a.copper[0].nets[0].polygons[0].outer);
    CHECK(b.copper[0].nets[0].polygons[0].holes == a.copper[0].nets[0].polygons[0].holes);
    CHECK(area_nm2(b.copper[0].nets[1].polygons) == doctest::Approx(area_nm2(a.copper[0].nets[1].polygons)));
    CHECK(b.pads.size() == 1);
    CHECK(b.pads[0].shape[0].outer == a.pads[0].shape[0].outer);
    REQUIRE(b.drills.size() == 2);
    CHECK(b.drills[1].slot_end == a.drills[1].slot_end);
    CHECK(b.drills[0].function == DrillFunction::Via);
    CHECK(b.components[0].rotation_deg == 90.0);
    CHECK(b.masks.size() == 1);
    CHECK(b.warnings == a.warnings);
}

}  // namespace

TEST_CASE("WKB round trip is exact to the nanometre") {
    const Board b = small_board();
    const PolygonSet& p = b.copper[0].nets[0].polygons;
    const std::string wkb = encode_wkb(p);
    const PolygonSet q = decode_wkb(wkb);
    REQUIRE(q.size() == 1);
    CHECK(q[0].outer == p[0].outer);
    CHECK(q[0].holes == p[0].holes);
    // Header: little endian, MultiPolygon, one polygon.
    CHECK(wkb[0] == 1);
    CHECK(wkb[1] == 6);
}

TEST_CASE("canonical model round trip through a directory and a zip archive") {
    const Board b = small_board();
    const fs::path dir = fs::temp_directory_path() / "kairoq_test_model";
    fs::remove_all(dir);
    save_board(b, dir.string());
    CHECK(is_board_model(dir.string()));
    check_same(b, load_board(dir.string()));
    const std::string zip = (fs::temp_directory_path() / "kairoq_test_model.zip").string();
    save_board(b, zip);
    CHECK(is_board_model(zip));
    check_same(b, load_board(zip));
    fs::remove_all(dir);
    fs::remove(zip);
}

TEST_CASE("display meshes cover the copper area") {
    const Board b = small_board();
    const Mesh m = board_mesh(b, "copper/1");
    // Top faces: triangles at the top z of the copper layer.
    const float ztop = static_cast<float>(to_mm(stackup_z(b)[0].top_nm));
    double area = 0.0;
    for (std::size_t t = 0; t + 2 < m.indices.size(); t += 3) {
        const float* a = &m.positions[3 * m.indices[t]];
        const float* p = &m.positions[3 * m.indices[t + 1]];
        const float* c = &m.positions[3 * m.indices[t + 2]];
        if (a[2] != ztop || p[2] != ztop || c[2] != ztop) continue;
        area += std::abs((p[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (p[1] - a[1])) / 2.0;
    }
    CHECK(area == doctest::Approx(100.0).epsilon(1e-6));
    CHECK(encode_mesh(m).size() == 8 + m.positions.size() * 4 + m.indices.size() * 4 + m.ids.size() * 4);
    for (const auto& key : mesh_keys(b)) CHECK_NOTHROW(board_mesh(b, key));
}

TEST_CASE("the ecc83 demo board imports from Gerber X2 and from IPC-2581") {
    const Board g = open_board({kData + "/kicad/ecc83/gerber"});
    CHECK(g.net_source == "x2");
    CHECK(g.stackup_source == "gbrjob");
    CHECK(g.copper.size() == 2);
    CHECK(g.drills.size() > 0);
    CHECK(g.pads.size() > 0);
    CHECK(g.components.size() > 0);
    CHECK(g.outline.size() == 1);
    CHECK(area_nm2(g.outline) / 1e12 == doctest::Approx(52.07 * 46.355).epsilon(1e-6));

    const Board x = open_board({kData + "/kicad/ecc83/ecc83.xml"});
    CHECK(x.net_source == "ipc2581");
    CHECK(x.copper.size() == 2);
    CHECK(x.drills.size() == g.drills.size());
    CHECK(area_nm2(x.outline) == doctest::Approx(area_nm2(g.outline)).epsilon(1e-6));
}
