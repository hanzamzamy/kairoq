// Planar geometry: integer coordinates in nanometres, polygons with holes, and the
// arc and stroke constructions shared by the importers.
#pragma once

#include <clipper2/clipper.h>

#include <cmath>
#include <cstdint>
#include <vector>

namespace kairoq {

using Clipper2Lib::Path64;
using Clipper2Lib::Paths64;
using Clipper2Lib::Point64;
using Clipper2Lib::Rect64;

// Internal length unit: 1 nm. Exchange formats use millimetres.
constexpr double kNmPerMm = 1.0e6;
constexpr double kPi = 3.14159265358979323846;

inline double to_mm(double nm) { return nm / kNmPerMm; }
inline std::int64_t nm_from_mm(double mm) { return std::llround(mm * kNmPerMm); }

// A polygon with holes. The outer ring has positive (counter-clockwise) orientation,
// the holes negative. Rings are not closed explicitly.
struct Polygon {
    Path64 outer;
    Paths64 holes;
};
using PolygonSet = std::vector<Polygon>;

struct GeometryOptions {
    // Largest distance between an arc and the chords that replace it (nm).
    double arc_tolerance_nm = 1000.0;
};

// Number of chords for a full circle of the given radius.
int circle_segments(double radius_nm, double tolerance_nm);

// Regular polygon replacing a circle. The vertex radius is enlarged so that the
// polygon has the area of the circle.
Path64 circle_path(Point64 centre, double radius_nm, double tolerance_nm);

// Points along an arc from start to end about centre, both ends included. When start
// equals end the arc is a full circle. The radius is interpolated linearly between the
// start and end radii, which absorbs rounding in the source data.
// With area_preserving, the interior points move outwards so that the chords enclose
// the area of the arc's sector (for boundaries of filled areas); otherwise they lie on
// the arc (for centre lines).
Path64 arc_points(Point64 start, Point64 end, Point64 centre, bool ccw, double tolerance_nm,
                  bool area_preserving = false);

// Segment stroked with a circle of the given radius (a capsule).
Path64 capsule(Point64 a, Point64 b, double radius_nm, double tolerance_nm);

// Open polyline stroked with a circle of the given radius: round joins and ends.
Paths64 stroke_polyline(const Path64& centreline, double radius_nm, double tolerance_nm);

// Closed ring stroked with a circle of the given radius (an annulus for a circle).
Paths64 stroke_ring(const Path64& ring, double radius_nm, double tolerance_nm);

// Segment stroked with a convex shape: the convex hull of the shape placed at both ends.
Path64 stroke_convex(Point64 a, Point64 b, const Path64& convex_shape);

Path64 convex_hull(Path64 points);

// Orientation helpers.
void make_positive(Path64& path);
double signed_area(const Path64& path);

// Union of paths under the non-zero fill rule, as polygons with holes.
PolygonSet union_polygons(const Paths64& paths);
PolygonSet polytree_to_polygons(const Clipper2Lib::PolyTree64& tree);
// Even-odd fill: nested rings alternate between material and holes.
PolygonSet even_odd_polygons(const Paths64& paths);

Paths64 to_paths(const PolygonSet& polygons);
PolygonSet difference(const PolygonSet& a, const PolygonSet& b);
PolygonSet intersection(const PolygonSet& a, const PolygonSet& b);
PolygonSet unite(const PolygonSet& a, const PolygonSet& b);

double area_nm2(const PolygonSet& polygons);
double area_nm2(const Polygon& polygon);
Rect64 bounds(const PolygonSet& polygons);
Rect64 bounds(const Polygon& polygon);

// True if the point lies inside the material of the polygon (on the boundary counts).
bool contains(const Polygon& polygon, Point64 point);
bool contains(const PolygonSet& polygons, Point64 point);

// Uniform grid over bounding boxes, for point queries against many polygons.
class GridIndex {
public:
    GridIndex(const std::vector<Rect64>& boxes, double cell_nm);
    // Indices of the boxes containing the point.
    std::vector<std::size_t> query(Point64 point) const;

private:
    std::vector<Rect64> boxes_;
    double cell_;
    std::int64_t x0_ = 0, y0_ = 0, nx_ = 0, ny_ = 0;
    std::vector<std::vector<std::size_t>> cells_;
};

// Affine placement used by flashes and IPC-2581 transforms:
// p' = offset + rotate(scale * mirror_x(p)).
struct Placement {
    double scale = 1.0;
    bool mirror_x = false;  // x -> -x before rotation
    bool mirror_y = false;  // y -> -y before rotation
    double rotation_deg = 0.0;  // counter-clockwise
    Point64 offset{0, 0};
};
Path64 transform(const Path64& path, const Placement& placement);
Paths64 transform(const Paths64& paths, const Placement& placement);
Point64 transform(Point64 point, const Placement& placement);

}  // namespace kairoq
