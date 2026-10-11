#include "kairoq/geometry.hpp"

#include <algorithm>
#include <cmath>

namespace kairoq {

using namespace Clipper2Lib;

namespace {

Point64 rounded(double x, double y) { return Point64(std::llround(x), std::llround(y)); }

}  // namespace

int circle_segments(double radius_nm, double tolerance_nm) {
    if (radius_nm <= tolerance_nm) return 8;
    const double step = 2.0 * std::acos(1.0 - tolerance_nm / radius_nm);
    const int n = static_cast<int>(std::ceil(2.0 * kPi / step));
    return std::clamp(n, 8, 4096);
}

Path64 circle_path(Point64 centre, double radius_nm, double tolerance_nm) {
    Path64 path;
    if (radius_nm <= 0.0) return path;
    const int n = circle_segments(radius_nm, tolerance_nm);
    const double step = 2.0 * kPi / n;
    // Equal-area radius: n/2 R^2 sin(2 pi/n) = pi r^2.
    const double r = radius_nm * std::sqrt(step / std::sin(step));
    path.reserve(n);
    for (int i = 0; i < n; ++i) {
        const double a = i * step;
        path.push_back(rounded(centre.x + r * std::cos(a), centre.y + r * std::sin(a)));
    }
    return path;
}

Path64 arc_points(Point64 start, Point64 end, Point64 centre, bool ccw, double tolerance_nm, bool area_preserving) {
    const double sx = static_cast<double>(start.x - centre.x);
    const double sy = static_cast<double>(start.y - centre.y);
    const double ex = static_cast<double>(end.x - centre.x);
    const double ey = static_cast<double>(end.y - centre.y);
    const double r0 = std::hypot(sx, sy);
    const double r1 = std::hypot(ex, ey);
    const double a0 = std::atan2(sy, sx);
    double sweep = std::atan2(ey, ex) - a0;
    if (ccw) {
        if (sweep <= 0.0) sweep += 2.0 * kPi;
    } else {
        if (sweep >= 0.0) sweep -= 2.0 * kPi;
    }
    if (start == end) sweep = ccw ? 2.0 * kPi : -2.0 * kPi;
    const double r = std::max(r0, r1);
    const int full = circle_segments(r, tolerance_nm);
    const int n = std::max(1, static_cast<int>(std::ceil(full * std::abs(sweep) / (2.0 * kPi))));
    // Interior radius factor k: the fan of chords has the sector area when
    // (n - 2) k^2 + 2 k = n theta / sin(theta), theta = |sweep| / n.
    double k = 1.0;
    if (area_preserving && n >= 2) {
        const double theta = std::abs(sweep) / n;
        const double rhs = n * theta / std::sin(theta);
        k = n == 2 ? rhs / 2.0 : (-2.0 + std::sqrt(4.0 + 4.0 * (n - 2) * rhs)) / (2.0 * (n - 2));
    }
    Path64 path;
    path.reserve(n + 1);
    path.push_back(start);
    for (int i = 1; i < n; ++i) {
        const double t = static_cast<double>(i) / n;
        const double a = a0 + sweep * t;
        const double ri = (r0 + (r1 - r0) * t) * k;
        path.push_back(rounded(centre.x + ri * std::cos(a), centre.y + ri * std::sin(a)));
    }
    path.push_back(end);
    return path;
}

Path64 capsule(Point64 a, Point64 b, double radius_nm, double tolerance_nm) {
    if (a == b) return circle_path(a, radius_nm, tolerance_nm);
    const double dx = static_cast<double>(b.x - a.x);
    const double dy = static_cast<double>(b.y - a.y);
    const double heading = std::atan2(dy, dx);
    const int full = circle_segments(radius_nm, tolerance_nm);
    const int half = std::max(2, full / 2);
    const double step = kPi / half;
    // The cap end points stay at the exact radius, so the straight part has the exact
    // width. The interior cap vertices sit at a radius ri chosen so that each cap has
    // the area of a half circle:
    //   1/2 sin(step) ((half - 2) ri^2 + 2 R ri) = pi R^2 / 2.
    const double R = radius_nm;
    const double m = half - 2;
    const double rhs = kPi * R * R / std::sin(step);
    const double ri = m > 0 ? (-2.0 * R + std::sqrt(4.0 * R * R + 4.0 * m * rhs)) / (2.0 * m) : R;
    Path64 path;
    path.reserve(2 * (half + 1));
    // Cap around b from heading - pi/2 to heading + pi/2, then around a.
    for (int i = 0; i <= half; ++i) {
        const double ang = heading - kPi / 2.0 + i * step;
        const double r = (i == 0 || i == half) ? R : ri;
        path.push_back(rounded(b.x + r * std::cos(ang), b.y + r * std::sin(ang)));
    }
    for (int i = 0; i <= half; ++i) {
        const double ang = heading + kPi / 2.0 + i * step;
        const double r = (i == 0 || i == half) ? R : ri;
        path.push_back(rounded(a.x + r * std::cos(ang), a.y + r * std::sin(ang)));
    }
    return path;
}

Paths64 stroke_polyline(const Path64& centreline, double radius_nm, double tolerance_nm) {
    if (centreline.size() == 2) return {capsule(centreline[0], centreline[1], radius_nm, tolerance_nm)};
    ClipperOffset offset(2.0, tolerance_nm);
    offset.AddPath(centreline, JoinType::Round, EndType::Round);
    Paths64 out;
    offset.Execute(radius_nm, out);
    return out;
}

Paths64 stroke_ring(const Path64& ring, double radius_nm, double tolerance_nm) {
    ClipperOffset offset(2.0, tolerance_nm);
    offset.AddPath(ring, JoinType::Round, EndType::Joined);
    Paths64 out;
    offset.Execute(radius_nm, out);
    return out;
}

Path64 convex_hull(Path64 pts) {
    std::sort(pts.begin(), pts.end(), [](const Point64& a, const Point64& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    if (pts.size() < 3) return pts;
    auto cross = [](const Point64& o, const Point64& a, const Point64& b) {
        return static_cast<double>(a.x - o.x) * static_cast<double>(b.y - o.y) -
               static_cast<double>(a.y - o.y) * static_cast<double>(b.x - o.x);
    };
    Path64 hull(2 * pts.size());
    size_t k = 0;
    for (size_t i = 0; i < pts.size(); ++i) {
        while (k >= 2 && cross(hull[k - 2], hull[k - 1], pts[i]) <= 0) --k;
        hull[k++] = pts[i];
    }
    for (size_t i = pts.size() - 1, t = k + 1; i > 0; --i) {
        while (k >= t && cross(hull[k - 2], hull[k - 1], pts[i - 1]) <= 0) --k;
        hull[k++] = pts[i - 1];
    }
    hull.resize(k - 1);
    return hull;
}

Path64 stroke_convex(Point64 a, Point64 b, const Path64& shape) {
    Path64 pts;
    pts.reserve(2 * shape.size());
    for (const auto& p : shape) {
        pts.emplace_back(p.x + a.x, p.y + a.y);
        pts.emplace_back(p.x + b.x, p.y + b.y);
    }
    return convex_hull(std::move(pts));
}

double signed_area(const Path64& path) { return Area(path); }

void make_positive(Path64& path) {
    if (Area(path) < 0.0) std::reverse(path.begin(), path.end());
}

namespace {

void collect(const PolyPath64& node, PolygonSet& out) {
    for (size_t i = 0; i < node.Count(); ++i) {
        const PolyPath64* outer = node.Child(i);
        Polygon poly;
        poly.outer = outer->Polygon();
        for (size_t j = 0; j < outer->Count(); ++j) {
            const PolyPath64* hole = outer->Child(j);
            poly.holes.push_back(hole->Polygon());
            collect(*hole, out);  // islands inside holes
        }
        out.push_back(std::move(poly));
    }
}

PolygonSet boolean(ClipType type, FillRule rule, const Paths64& subject, const Paths64& clip) {
    Clipper64 clipper;
    clipper.AddSubject(subject);
    if (!clip.empty()) clipper.AddClip(clip);
    PolyTree64 tree;
    clipper.Execute(type, rule, tree);
    return polytree_to_polygons(tree);
}

}  // namespace

PolygonSet polytree_to_polygons(const PolyTree64& tree) {
    PolygonSet out;
    collect(tree, out);
    return out;
}

PolygonSet union_polygons(const Paths64& paths) {
    return boolean(ClipType::Union, FillRule::NonZero, paths, {});
}

PolygonSet even_odd_polygons(const Paths64& paths) {
    return boolean(ClipType::Union, FillRule::EvenOdd, paths, {});
}

Paths64 to_paths(const PolygonSet& polygons) {
    Paths64 out;
    for (const auto& p : polygons) {
        out.push_back(p.outer);
        out.insert(out.end(), p.holes.begin(), p.holes.end());
    }
    return out;
}

PolygonSet difference(const PolygonSet& a, const PolygonSet& b) {
    if (b.empty()) return a;
    return boolean(ClipType::Difference, FillRule::NonZero, to_paths(a), to_paths(b));
}

PolygonSet intersection(const PolygonSet& a, const PolygonSet& b) {
    if (a.empty() || b.empty()) return {};
    return boolean(ClipType::Intersection, FillRule::NonZero, to_paths(a), to_paths(b));
}

PolygonSet unite(const PolygonSet& a, const PolygonSet& b) {
    return boolean(ClipType::Union, FillRule::NonZero, to_paths(a), to_paths(b));
}

double area_nm2(const Polygon& polygon) {
    double a = Area(polygon.outer);
    for (const auto& h : polygon.holes) a += Area(h);
    return a;
}

double area_nm2(const PolygonSet& polygons) {
    double a = 0.0;
    for (const auto& p : polygons) a += area_nm2(p);
    return a;
}

Rect64 bounds(const Polygon& polygon) { return GetBounds(polygon.outer); }

Rect64 bounds(const PolygonSet& polygons) {
    Rect64 r(INT64_MAX, INT64_MAX, INT64_MIN, INT64_MIN);
    for (const auto& p : polygons) {
        const Rect64 b = GetBounds(p.outer);
        r.left = std::min(r.left, b.left);
        r.top = std::min(r.top, b.top);
        r.right = std::max(r.right, b.right);
        r.bottom = std::max(r.bottom, b.bottom);
    }
    return r;
}

bool contains(const Polygon& polygon, Point64 point) {
    const Rect64 b = GetBounds(polygon.outer);
    if (point.x < b.left || point.x > b.right || point.y < b.top || point.y > b.bottom) return false;
    if (PointInPolygon(point, polygon.outer) == PointInPolygonResult::IsOutside) return false;
    for (const auto& h : polygon.holes)
        if (PointInPolygon(point, h) == PointInPolygonResult::IsInside) return false;
    return true;
}

bool contains(const PolygonSet& polygons, Point64 point) {
    for (const auto& p : polygons)
        if (contains(p, point)) return true;
    return false;
}

GridIndex::GridIndex(const std::vector<Rect64>& boxes, double cell_nm) : boxes_(boxes), cell_(cell_nm) {
    if (boxes_.empty()) return;
    Rect64 all(INT64_MAX, INT64_MAX, INT64_MIN, INT64_MIN);
    for (const auto& b : boxes_) {
        all.left = std::min(all.left, b.left);
        all.top = std::min(all.top, b.top);
        all.right = std::max(all.right, b.right);
        all.bottom = std::max(all.bottom, b.bottom);
    }
    // Keep the grid below about four million cells.
    const double span = std::max<double>(all.right - all.left, all.bottom - all.top);
    cell_ = std::max(cell_, span / 2000.0);
    x0_ = all.left;
    y0_ = all.top;
    nx_ = static_cast<std::int64_t>((all.right - all.left) / cell_) + 1;
    ny_ = static_cast<std::int64_t>((all.bottom - all.top) / cell_) + 1;
    cells_.resize(static_cast<std::size_t>(nx_ * ny_));
    for (std::size_t i = 0; i < boxes_.size(); ++i) {
        const auto& b = boxes_[i];
        const auto ix0 = static_cast<std::int64_t>((b.left - x0_) / cell_);
        const auto ix1 = static_cast<std::int64_t>((b.right - x0_) / cell_);
        const auto iy0 = static_cast<std::int64_t>((b.top - y0_) / cell_);
        const auto iy1 = static_cast<std::int64_t>((b.bottom - y0_) / cell_);
        for (auto iy = iy0; iy <= iy1; ++iy)
            for (auto ix = ix0; ix <= ix1; ++ix) cells_[static_cast<std::size_t>(iy * nx_ + ix)].push_back(i);
    }
}

std::vector<std::size_t> GridIndex::query(Point64 p) const {
    std::vector<std::size_t> out;
    if (cells_.empty()) return out;
    const auto ix = static_cast<std::int64_t>(std::floor((p.x - x0_) / cell_));
    const auto iy = static_cast<std::int64_t>(std::floor((p.y - y0_) / cell_));
    if (ix < 0 || iy < 0 || ix >= nx_ || iy >= ny_) return out;
    for (std::size_t i : cells_[static_cast<std::size_t>(iy * nx_ + ix)]) {
        const auto& b = boxes_[i];
        if (p.x >= b.left && p.x <= b.right && p.y >= b.top && p.y <= b.bottom) out.push_back(i);
    }
    return out;
}

Point64 transform(Point64 p, const Placement& t) {
    double x = static_cast<double>(p.x) * t.scale;
    double y = static_cast<double>(p.y) * t.scale;
    if (t.mirror_x) x = -x;
    if (t.mirror_y) y = -y;
    if (t.rotation_deg != 0.0) {
        const double a = t.rotation_deg * kPi / 180.0;
        // Exact quarter turns keep integer coordinates exact.
        const double q = t.rotation_deg / 90.0;
        if (q == std::floor(q)) {
            const int k = ((static_cast<int>(q) % 4) + 4) % 4;
            const double xs = x, ys = y;
            if (k == 1) { x = -ys; y = xs; }
            else if (k == 2) { x = -xs; y = -ys; }
            else if (k == 3) { x = ys; y = -xs; }
        } else {
            const double c = std::cos(a), s = std::sin(a);
            const double xr = c * x - s * y;
            y = s * x + c * y;
            x = xr;
        }
    }
    return rounded(x + t.offset.x, y + t.offset.y);
}

Path64 transform(const Path64& path, const Placement& t) {
    Path64 out;
    out.reserve(path.size());
    for (const auto& p : path) out.push_back(transform(p, t));
    if (t.mirror_x != t.mirror_y) std::reverse(out.begin(), out.end());
    return out;
}

Paths64 transform(const Paths64& paths, const Placement& t) {
    Paths64 out;
    out.reserve(paths.size());
    for (const auto& p : paths) out.push_back(transform(p, t));
    return out;
}

}  // namespace kairoq
