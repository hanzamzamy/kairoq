#include "kairoq/mesh2d.hpp"

#include <doctest/doctest.h>

#include <cmath>

using namespace kairoq;

namespace {

Path64 rect(std::int64_t x0, std::int64_t y0, std::int64_t x1, std::int64_t y1) {
    return {Point64(x0, y0), Point64(x1, y0), Point64(x1, y1), Point64(x0, y1)};
}

double mesh_area(const TriMesh& m) {
    double a = 0;
    for (std::size_t t = 0; t < m.triangles.size(); ++t) a += triangle_area(m, t);
    return a;
}

}  // namespace

TEST_CASE("quality mesh of a plate with a hole and a port on its edge") {
    // 10 x 5 mm plate, 1 mm square hole, 1 x 2 mm port flush with the left edge.
    Polygon plate{rect(0, 0, 10000000, 5000000), {}};
    Path64 hole = rect(4000000, 2000000, 5000000, 3000000);
    std::reverse(hole.begin(), hole.end());
    plate.holes.push_back(hole);
    const PolygonSet port{Polygon{rect(0, 1500000, 1000000, 3500000), {}}};
    MeshSettings s;
    s.boundary_tolerance_nm = 0;
    s.max_edge_nm = 500000;
    const TriMesh m = mesh_region(plate, port, s);
    REQUIRE(m.triangles.size() > 100);
    CHECK(mesh_area(m) == doctest::Approx(50e12 - 1e12).epsilon(1e-12));
    double worst = 180;
    for (std::size_t t = 0; t < m.triangles.size(); ++t) worst = std::min(worst, min_angle_deg(m, t));
    CHECK(worst >= 19.9);
    // The port is an exact union of triangles.
    double in_port = 0;
    for (std::size_t t = 0; t < m.triangles.size(); ++t) {
        const auto& tr = m.triangles[t];
        const double cx = (m.x[tr[0]] + m.x[tr[1]] + m.x[tr[2]]) / 3, cy = (m.y[tr[0]] + m.y[tr[1]] + m.y[tr[2]]) / 3;
        if (cx < 1000000 && cy > 1500000 && cy < 3500000) in_port += triangle_area(m, t);
    }
    CHECK(in_port == doctest::Approx(2e12).epsilon(1e-12));
    // No triangle is degenerate and every vertex is used.
    std::vector<int> used(m.x.size(), 0);
    for (std::size_t t = 0; t < m.triangles.size(); ++t) {
        CHECK(triangle_area(m, t) > 0);
        for (auto v : m.triangles[t]) used[v] = 1;
    }
    for (int u : used) CHECK(u == 1);
}

TEST_CASE("seeds refine the mesh where they are placed") {
    const Polygon plate{rect(0, 0, 10000000, 10000000), {}};
    MeshSettings s;
    s.max_edge_nm = 2000000;
    s.grading = 0;
    const TriMesh coarse = mesh_region(plate, {}, s);
    std::vector<std::array<double, 3>> seeds;
    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 20; ++j) seeds.push_back({1e5 + i * 1e5, 1e5 + j * 1e5, 1e5});
    const TriMesh fine = mesh_region(plate, {}, s, seeds);
    CHECK(fine.triangles.size() > coarse.triangles.size() + 600);
    CHECK(mesh_area(fine) == doctest::Approx(1e14).epsilon(1e-12));
    // Seeds on top of each other or on the boundary are dropped, not duplicated.
    const std::vector<std::array<double, 3>> bad{{5e6, 5e6, 1e5}, {5e6 + 1, 5e6, 1e5}, {0.0, 5e6, 1e5}};
    const TriMesh m = mesh_region(plate, {}, s, bad);
    for (std::size_t t = 0; t < m.triangles.size(); ++t) CHECK(triangle_area(m, t) > 1e6);
}
