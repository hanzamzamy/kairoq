#include "kairoq/feature.hpp"
#include "kairoq/geometry.hpp"

#include <doctest/doctest.h>

#include <cmath>

using namespace kairoq;

namespace {
PolygonSet one(const Path64& p) { return union_polygons({p}); }
}  // namespace

TEST_CASE("circle polygons keep the circle area") {
    for (double r : {0.05e6, 0.3e6, 2.0e6, 25.0e6}) {
        const Path64 c = circle_path(Point64(0, 0), r, 1000.0);
        CHECK(signed_area(c) > 0);
        // Vertices are rounded to 1 nm, which bounds the agreement for small circles.
        CHECK(signed_area(c) == doctest::Approx(kPi * r * r).epsilon(2e-5));
    }
}

TEST_CASE("area-preserving arcs enclose the sector area") {
    const double r = 2e6;
    for (bool ccw : {true, false}) {
        Path64 ring = arc_points(Point64(2000000, 0), Point64(-2000000, 0), Point64(0, 0), ccw, 1000.0, true);
        const Path64 back = arc_points(Point64(-2000000, 0), Point64(2000000, 0), Point64(0, 0), ccw, 1000.0, true);
        ring.insert(ring.end(), back.begin() + 1, back.end() - 1);
        CHECK(std::abs(signed_area(ring)) == doctest::Approx(kPi * r * r).epsilon(1e-6));
    }
}

TEST_CASE("capsules have the exact width and the area of a stadium") {
    const double r = 0.4e6, len = 10e6;
    const Path64 c = capsule(Point64(0, 0), Point64(std::llround(len), 0LL), r, 1000.0);
    CHECK(signed_area(c) > 0);
    CHECK(signed_area(c) == doctest::Approx(2 * r * len + kPi * r * r).epsilon(1e-6));
    const Rect64 b = Clipper2Lib::GetBounds(c);
    // The straight part is exactly 2r wide; the caps stay within the tolerance.
    CHECK(std::abs((b.bottom - b.top) - 2 * r) <= 2000);
    // Diagonal capsule: same area.
    const Path64 d = capsule(Point64(0, 0), Point64(7000000, 7000000), r, 1000.0);
    CHECK(signed_area(d) == doctest::Approx(2 * r * std::hypot(7e6, 7e6) + kPi * r * r).epsilon(1e-6));
}

TEST_CASE("arc points follow the radius and end exactly") {
    const Path64 a = arc_points(Point64(10000000, 0), Point64(0, 10000000), Point64(0, 0), true, 1000.0);
    CHECK(a.front() == Point64(10000000, 0));
    CHECK(a.back() == Point64(0, 10000000));
    for (const auto& p : a) CHECK(std::hypot(double(p.x), double(p.y)) == doctest::Approx(1e7).epsilon(1e-6));
    // Chord deviation within tolerance: sagitta of each step.
    for (std::size_t i = 0; i + 1 < a.size(); ++i) {
        const double mx = (a[i].x + a[i + 1].x) / 2.0, my = (a[i].y + a[i + 1].y) / 2.0;
        CHECK(1e7 - std::hypot(mx, my) <= 1000.0 + 1.0);
    }
    // Clockwise from the same start goes the long way round.
    const Path64 cw = arc_points(Point64(10000000, 0), Point64(0, 10000000), Point64(0, 0), false, 1000.0);
    CHECK(cw.size() > 2 * a.size());
}

TEST_CASE("boolean operations and point containment") {
    const Path64 sq{Point64(0, 0), Point64(10, 0), Point64(10, 10), Point64(0, 10)};
    const Path64 inner{Point64(2, 2), Point64(8, 2), Point64(8, 8), Point64(2, 8)};
    const PolygonSet ring = difference(one(sq), one(inner));
    REQUIRE(ring.size() == 1);
    CHECK(ring[0].holes.size() == 1);
    CHECK(area_nm2(ring) == doctest::Approx(64.0));
    CHECK(contains(ring, Point64(1, 1)));
    CHECK_FALSE(contains(ring, Point64(5, 5)));
    CHECK(contains(ring, Point64(0, 5)));  // boundary
}

TEST_CASE("grid index returns the boxes containing a point") {
    std::vector<Rect64> boxes{Rect64(0, 0, 10, 10), Rect64(5, 5, 20, 20), Rect64(100, 100, 110, 110)};
    GridIndex index(boxes, 3.0);
    auto hits = index.query(Point64(7, 7));
    CHECK(hits.size() == 2);
    CHECK(index.query(Point64(105, 105)) == std::vector<std::size_t>{2});
    CHECK(index.query(Point64(50, 50)).empty());
    CHECK(index.query(Point64(-5, 0)).empty());
}

TEST_CASE("placement transforms mirror, rotate and offset") {
    Placement t;
    t.rotation_deg = 90;
    t.offset = Point64(100, 0);
    CHECK(transform(Point64(10, 0), t) == Point64(100, 10));
    t.mirror_x = true;
    CHECK(transform(Point64(10, 0), t) == Point64(100, -10));
    // Mirroring reverses the ring order so that the orientation stays positive.
    const Path64 sq{Point64(0, 0), Point64(10, 0), Point64(10, 10), Point64(0, 10)};
    CHECK(signed_area(transform(sq, t)) > 0);
}

TEST_CASE("outline chaining joins unordered segments and keeps cut-outs") {
    auto seg = [](Point64 a, Point64 b) {
        Feature f;
        f.kind = FeatureKind::Stroke;
        f.at = a;
        f.end = b;
        return f;
    };
    const std::int64_t L = 1000000;
    auto P = [](std::int64_t x, std::int64_t y) { return Point64(x, y); };
    std::vector<Feature> fs{seg(P(0, 0), P(L, 0)), seg(P(L, L), P(0, L)), seg(P(L, 0), P(L, L)), seg(P(0, 0), P(0, L))};
    Feature hole;
    hole.kind = FeatureKind::Stroke;
    hole.at = hole.end = Point64(600000, 500000);
    hole.arc_centre = Point64(500000, 500000);
    hole.arc_ccw = true;
    fs.push_back(hole);
    const auto out = outline_from_features(fs, 10.0, 1.0);
    CHECK(out.open_chains == 0);
    REQUIRE(out.polygons.size() == 1);
    CHECK(out.polygons[0].holes.size() == 1);
    CHECK(area_nm2(out.polygons) == doctest::Approx(1e12 - kPi * 1e10).epsilon(1e-4));
}
