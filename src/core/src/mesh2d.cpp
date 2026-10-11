#include "kairoq/mesh2d.hpp"

#include "segment_grid.hpp"

#include <CDT.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <unordered_set>

namespace kairoq {

using namespace Clipper2Lib;

namespace {

using detail::SegmentGrid;

// Element size at a small closed constraint ring, relative to its equivalent radius.
constexpr double kSmallRingFraction = 0.15;
// Element size at sharp re-entrant outline corners, relative to the minimum size, and
// the turning angle above which a corner counts as sharp (20 degrees).
constexpr double kCornerFraction = 0.2;
constexpr double kCornerTurn = 20.0 * kPi / 180.0;


}  // namespace

PolygonSet simplified_region(const Polygon& region, const MeshSettings& s) {
    Paths64 rings{region.outer};
    rings.insert(rings.end(), region.holes.begin(), region.holes.end());
    if (s.boundary_tolerance_nm > 0) rings = SimplifyPaths(rings, s.boundary_tolerance_nm, true);
    // Simplification can make thin parts self-intersect; the union repairs them.
    return union_polygons(rings);
}

TriMesh mesh_region(const Polygon& region, const PolygonSet& constraints, const MeshSettings& s,
                    const std::vector<std::array<double, 3>>& seeds) {
    TriMesh out;
    const PolygonSet domain = simplified_region(region, s);
    if (domain.empty()) return out;
    // Each constraint separately: a union would swallow rings nested in other
    // constraints (a hole wall inside a port pad). Only the outline rings whose boxes
    // overlap a constraint's can cut it; a plane has as many rings as clearances.
    std::vector<std::pair<Rect64, const Path64*>> domain_rings;
    for (const auto& p : domain) {
        domain_rings.emplace_back(GetBounds(p.outer), &p.outer);
        for (const auto& hh : p.holes) domain_rings.emplace_back(GetBounds(hh), &hh);
    }
    PolygonSet inner;
    for (const auto& c : constraints) {
        const Rect64 cb = bounds(c);
        Paths64 near;
        for (const auto& [rb, ring] : domain_rings)
            if (rb.Intersects(cb)) near.push_back(*ring);
        if (near.empty()) continue;
        Clipper64 clipper;
        clipper.AddSubject(near);
        clipper.AddClip(to_paths({c}));
        PolyTree64 tree;
        clipper.Execute(ClipType::Intersection, FillRule::NonZero, tree);
        const PolygonSet part = polytree_to_polygons(tree);
        inner.insert(inner.end(), part.begin(), part.end());
    }

    // Local coordinates in mm about the lower-left corner keep the doubles well scaled.
    using P2 = CDT::V2d<double>;
    const Rect64 box = bounds(domain);
    const double ox = static_cast<double>(box.left), oy = static_cast<double>(box.top);
    const double w_mm = (box.right - box.left) * 1e-6, h_mm = (box.bottom - box.top) * 1e-6;
    auto local = [&](const Point64& p) { return P2((p.x - ox) * 1e-6, (p.y - oy) * 1e-6); };
    const double hmin = std::max(s.min_size_nm, 1.0) * 1e-6, hmax = std::max(s.max_edge_nm * 1e-6, hmin);

    std::vector<std::vector<P2>> rings;  // region outline
    std::vector<char> outer_ring;
    for (const auto& p : domain) {
        for (const Path64* r : [&] {
                 std::vector<const Path64*> v{&p.outer};
                 for (const auto& hh : p.holes) v.push_back(&hh);
                 return v;
             }()) {
            if (r->size() < 3) continue;
            std::vector<P2> ring;
            for (const auto& q : *r) ring.push_back(local(q));
            rings.push_back(std::move(ring));
            outer_ring.push_back(r == &p.outer);
        }
    }
    const double h_outer = s.outer_size_nm > 0 ? s.outer_size_nm * 1e-6 : hmin;
    std::vector<SegmentGrid::Seg> outline_segs;
    std::vector<std::pair<std::size_t, std::size_t>> outline_id;
    for (std::size_t r = 0; r < rings.size(); ++r)
        for (std::size_t k = 0; k < rings[r].size(); ++k) {
            const P2& a = rings[r][k];
            const P2& c = rings[r][(k + 1) % rings[r].size()];
            outline_segs.push_back({a.x, a.y, c.x, c.y, true, hmin});
            outline_id.emplace_back(r, k);
        }
    const SegmentGrid outline_grid(outline_segs, 0.0, 0.0, w_mm, h_mm, std::max(hmin, 0.02));

    // Constraint segments. Those lying on the outline are dropped (the outline is already
    // a mesh edge); ends near it are moved onto it, and the outline segment is split
    // there, so that no two constraint edges overlap or run side by side. "Near" covers
    // the outline simplification, which can move the outline off a coinciding edge.
    const double kSnap = std::max(1e-5, 1.5 * s.boundary_tolerance_nm * 1e-6);  // mm
    struct ConstraintSeg {
        P2 a, b;
        double h;
    };
    std::vector<ConstraintSeg> csegs;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<double, P2>>> splits;
    auto snap = [&](P2& p) {
        const auto near = outline_grid.nearest(p.x, p.y, kSnap);
        if (!near) return;
        const auto& g = outline_segs[near->first];
        p = P2(g.x0 + near->second * (g.x1 - g.x0), g.y0 + near->second * (g.y1 - g.y0));
        if (near->second > 0.0 && near->second < 1.0) splits[outline_id[near->first]].push_back({near->second, p});
    };
    for (const auto& p : inner) {
        std::vector<const Path64*> v{&p.outer};
        for (const auto& hh : p.holes) v.push_back(&hh);
        for (const Path64* r : v) {
            if (r->size() < 3) continue;
            // Finer size at small closed rings (hole walls, small ports), where the
            // current converges: a fixed fraction of the ring's equivalent radius.
            const double radius = std::sqrt(std::abs(Area(*r)) / kPi) * 1e-6;
            const double h = s.constraint_size_nm > 0 ? s.constraint_size_nm * 1e-6
                                                      : std::clamp(kSmallRingFraction * radius, hmin / 8, hmin);
            for (std::size_t k = 0; k < r->size(); ++k) {
                P2 a = local((*r)[k]), b = local((*r)[(k + 1) % r->size()]);
                const P2 m((a.x + b.x) / 2, (a.y + b.y) / 2);
                if (outline_grid.any_within(a.x, a.y, kSnap) && outline_grid.any_within(b.x, b.y, kSnap) &&
                    outline_grid.any_within(m.x, m.y, kSnap))
                    continue;
                snap(a);
                snap(b);
                csegs.push_back({a, b, h});
            }
        }
    }

    std::vector<P2> verts;
    std::vector<CDT::Edge> edges;
    std::vector<SegmentGrid::Seg> segs;
    for (std::size_t r = 0; r < rings.size(); ++r) {
        std::vector<P2> pts;
        for (std::size_t k = 0; k < rings[r].size(); ++k) {
            pts.push_back(rings[r][k]);
            auto it = splits.find({r, k});
            if (it == splits.end()) continue;
            auto sp = it->second;
            std::sort(sp.begin(), sp.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
            for (const auto& q : sp) pts.push_back(q.second);
        }
        const auto base = static_cast<CDT::VertInd>(verts.size());
        verts.insert(verts.end(), pts.begin(), pts.end());
        const auto n = static_cast<CDT::VertInd>(pts.size());
        for (CDT::VertInd i = 0; i < n; ++i) {
            edges.emplace_back(base + i, base + (i + 1) % n);
            const P2& a = pts[i];
            const P2& c = pts[(i + 1) % n];
            double h = outer_ring[r] ? h_outer : hmin;
            if (s.boundary_size) h = s.boundary_size(0.5 * (a.x + c.x) * 1e6 + ox, 0.5 * (a.y + c.y) * 1e6 + oy) * 1e-6;
            segs.push_back({a.x, a.y, c.x, c.y, true, h});
        }
    }
    double h_floor = hmin;
    for (const auto& c : csegs) {
        const auto base = static_cast<CDT::VertInd>(verts.size());
        verts.push_back(c.a);
        verts.push_back(c.b);
        edges.emplace_back(base, base + 1);
        segs.push_back({c.a.x, c.a.y, c.b.x, c.b.y, false, c.h});
        h_floor = std::min(h_floor, c.h);
    }
    // Sharp re-entrant corners of the outline (interior angle above 200 degrees) are
    // field singularities: they get the finest size, as zero-length segments.
    {
        const double h_corner = hmin * kCornerFraction;
        auto corners = [&](const Path64& ring, bool hole) {
            const std::size_t n = ring.size();
            for (std::size_t k = 0; k < n; ++k) {
                const Point64& p0 = ring[(k + n - 1) % n];
                const Point64& p1 = ring[k];
                const Point64& p2 = ring[(k + 1) % n];
                const double ux = double(p1.x - p0.x), uy = double(p1.y - p0.y);
                const double vx = double(p2.x - p1.x), vy = double(p2.y - p1.y);
                // Turning angle; outer rings are counter-clockwise, holes clockwise, so a
                // right turn along an outer ring (left along a hole) is re-entrant.
                const double turn = std::atan2(ux * vy - uy * vx, ux * vx + uy * vy);
                if ((hole ? turn : -turn) > kCornerTurn) {
                    const double x = (p1.x - ox) * 1e-6, y = (p1.y - oy) * 1e-6;
                    const double h = s.boundary_size ? kCornerFraction * s.boundary_size(double(p1.x), double(p1.y)) * 1e-6 : h_corner;
                    segs.push_back({x, y, x, y, false, h});
                    h_floor = std::min(h_floor, h);
                }
            }
        };
        for (const auto& p : domain) {
            corners(p.outer, false);
            for (const auto& hh : p.holes) corners(hh, true);
        }
    }
    const SegmentGrid grid(std::move(segs), 0.0, 0.0, w_mm, h_mm, std::max(hmin, 0.02));

    // Steiner points from the size field: quadtree cells are split while larger than
    // h at their centre; leaf centres inside the region and clear of the boundaries
    // become vertices. The quality refinement afterwards completes the grading.
    // With seeds (a refinement of an earlier mesh), the seeds already carry this grading.
    if (s.grading > 0 && seeds.empty()) {
        struct Cell {
            double x, y, size;
        };
        std::vector<Cell> stack{{0.0, 0.0, std::max(w_mm, h_mm)}};
        const std::size_t budget = s.max_vertices;
        std::size_t seeded = 0;
        while (!stack.empty() && seeded < budget) {
            const Cell c = stack.back();
            stack.pop_back();
            const double cx = c.x + c.size / 2, cy = c.y + c.size / 2;
            if (cx - c.size / 2 > w_mm || cy - c.size / 2 > h_mm) continue;
            const double half_diag = c.size * 0.7072;
            const bool crossed = grid.any_within(cx, cy, half_diag);
            if (!crossed && !grid.inside(cx, cy)) continue;  // entirely outside
            const bool split = c.size > hmax || (c.size > h_floor && grid.finer_than(cx, cy, c.size, s.grading, h_floor));
            if (split) {
                const double h = c.size / 2;
                stack.push_back({c.x, c.y, h});
                stack.push_back({c.x + h, c.y, h});
                stack.push_back({c.x, c.y + h, h});
                stack.push_back({c.x + h, c.y + h, h});
            } else if (!grid.any_within(cx, cy, 0.5 * c.size) && grid.inside(cx, cy)) {
                verts.emplace_back(cx, cy);
                ++seeded;
            }
        }
    }
    // Seeds clear of the boundaries by a small margin, and of each other by 10 nm;
    // closer points would only make slivers.
    {
        std::unordered_set<std::uint64_t> taken;
        auto key = [](double x, double y) {
            return (static_cast<std::uint64_t>(std::llround(x * 1e5)) << 32) ^ static_cast<std::uint64_t>(std::llround(y * 1e5));
        };
        auto near_taken = [&](double x, double y) {
            const auto kx = std::llround(x * 1e5), ky = std::llround(y * 1e5);
            for (long long dx = -1; dx <= 1; ++dx)
                for (long long dy = -1; dy <= 1; ++dy)
                    if (taken.count((static_cast<std::uint64_t>(kx + dx) << 32) ^ static_cast<std::uint64_t>(ky + dy))) return true;
            return false;
        };
        for (const auto& v : verts) taken.insert(key(v.x, v.y));
        for (const auto& p : seeds) {
            const double x = (p[0] - ox) * 1e-6, y = (p[1] - oy) * 1e-6;
            if (x < 0 || y < 0 || x > w_mm || y > h_mm) continue;
            if (near_taken(x, y) || !grid.inside(x, y) || grid.any_within(x, y, std::max(1e-4, 0.3e-6 * p[2]))) continue;
            taken.insert(key(x, y));
            verts.emplace_back(x, y);
        }
    }
    CDT::RemoveDuplicatesAndRemapEdges(verts, edges);

    CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::TryResolve, 1e-7);
    cdt.insertVertices(verts);
    cdt.insertEdges(edges);

    // Triangles outside the region are excluded by their centroid.
    CDT::TriIndUSet erase = cdt.collectSuperTriangle();
    for (CDT::TriInd t = 0; t < cdt.triangles.size(); ++t) {
        if (erase.count(t)) continue;
        const auto& v = cdt.triangles[t].vertices;
        const double cx = (cdt.vertices[v[0]].x + cdt.vertices[v[1]].x + cdt.vertices[v[2]].x) / 3.0;
        const double cy = (cdt.vertices[v[0]].y + cdt.vertices[v[1]].y + cdt.vertices[v[2]].y) / 3.0;
        // Degenerate triangles (collinear outline points) have their centroid on the
        // outline, where the inside test is ambiguous; they carry no area and go.
        const double ux = cdt.vertices[v[1]].x - cdt.vertices[v[0]].x, uy = cdt.vertices[v[1]].y - cdt.vertices[v[0]].y;
        const double wx = cdt.vertices[v[2]].x - cdt.vertices[v[0]].x, wy = cdt.vertices[v[2]].y - cdt.vertices[v[0]].y;
        const double area2 = std::abs(ux * wy - uy * wx);
        const double longest2 = std::max({ux * ux + uy * uy, wx * wx + wy * wy, (wx - ux) * (wx - ux) + (wy - uy) * (wy - uy)});
        if (area2 <= 1e-9 * longest2 || !grid.inside(cx, cy)) erase.insert(t);
    }

    const double min_edge_mm = s.min_edge_nm * 1e-6;
    const auto budget = static_cast<CDT::VertInd>(s.max_vertices);
    CDT::Unrefined u1 = cdt.refineTriangles(budget, CDT::RefinementCriterion::SmallestAngle,
                                            s.min_angle_deg * kPi / 180.0, &erase, min_edge_mm);
    const double max_area_mm2 = std::sqrt(3.0) / 4.0 * std::pow(s.max_edge_nm * 1e-6, 2);
    CDT::Unrefined u2 =
        cdt.refineTriangles(budget, CDT::RefinementCriterion::LargestArea, max_area_mm2, &erase, min_edge_mm);
    out.unrefined = u1.circumcenterOutside + u1.splitVertexInvalid + u2.circumcenterOutside + u2.splitVertexInvalid;
    cdt.finalizeTriangulation(erase);

    // Keep only vertices used by a triangle.
    std::vector<std::uint32_t> index(cdt.vertices.size(), UINT32_MAX);
    out.triangles.reserve(cdt.triangles.size());
    for (const auto& t : cdt.triangles) {
        {
            const auto& a = cdt.vertices[t.vertices[0]];
            const auto& b = cdt.vertices[t.vertices[1]];
            const auto& c = cdt.vertices[t.vertices[2]];
            const double ux = b.x - a.x, uy = b.y - a.y, wx = c.x - a.x, wy = c.y - a.y;
            const double longest2 = std::max({ux * ux + uy * uy, wx * wx + wy * wy, (wx - ux) * (wx - ux) + (wy - uy) * (wy - uy)});
            if (ux * wy - uy * wx <= 1e-9 * longest2) continue;
        }
        std::array<std::uint32_t, 3> tri{};
        for (int k = 0; k < 3; ++k) {
            const auto v = static_cast<std::size_t>(t.vertices[k]);
            if (index[v] == UINT32_MAX) {
                index[v] = static_cast<std::uint32_t>(out.x.size());
                out.x.push_back(cdt.vertices[v].x * 1e6 + ox);
                out.y.push_back(cdt.vertices[v].y * 1e6 + oy);
            }
            tri[k] = index[v];
        }
        out.triangles.push_back(tri);
    }
    return out;
}

double triangle_area(const TriMesh& m, std::size_t t) {
    const auto& v = m.triangles[t];
    return 0.5 * std::abs((m.x[v[1]] - m.x[v[0]]) * (m.y[v[2]] - m.y[v[0]]) -
                          (m.x[v[2]] - m.x[v[0]]) * (m.y[v[1]] - m.y[v[0]]));
}

double min_angle_deg(const TriMesh& m, std::size_t t) {
    const auto& v = m.triangles[t];
    double best = 180.0;
    for (int k = 0; k < 3; ++k) {
        const auto a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3];
        const double ux = m.x[b] - m.x[a], uy = m.y[b] - m.y[a];
        const double wx = m.x[c] - m.x[a], wy = m.y[c] - m.y[a];
        const double ang = std::atan2(std::abs(ux * wy - uy * wx), ux * wx + uy * wy) * 180.0 / kPi;
        best = std::min(best, ang);
    }
    return best;
}

}  // namespace kairoq
