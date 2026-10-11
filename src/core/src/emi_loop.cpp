// Loop inductance by partial inductances (PEEC) in the high-frequency limit: the copper
// around a loop is a set of thin sheets with Rao-Wilton-Glisson (lowest-order
// Raviart-Thomas) currents, vias are tubes, and the loop's components straight
// conductors; the current distribution minimises the magnetic energy under Kirchhoff's
// current law. Method: docs/solvers/emi.md.
#include "kairoq/dc.hpp"
#include "kairoq/emi.hpp"
#include "kairoq/mesh2d.hpp"

#include "fem_util.hpp"
#include "emi_loop.hpp"
#include "peec.hpp"

#include <algorithm>
#include <numeric>
#include <atomic>
#include <cmath>
#include <version>
#include <cstring>
#include <cstdint>
#include <complex>
#include <map>
#include <memory>
#include <tuple>
#include <set>
#include <stdexcept>

namespace kairoq {

namespace {

// Adds v to x in place, safely against other threads: std::atomic_ref where the standard
// library has it, otherwise a compare-and-swap on the value's bits (GCC and Clang builtins;
// Apple's library lacks atomic_ref).
inline void atomic_add(double& x, double v) {
#if defined(__cpp_lib_atomic_ref) && !defined(KAIROQ_TEST_ATOMIC_FALLBACK)
    std::atomic_ref<double>(x).fetch_add(v);
#else
    auto* p = reinterpret_cast<std::uint64_t*>(&x);
    std::uint64_t old = __atomic_load_n(p, __ATOMIC_RELAXED), next = 0;
    do {
        double d;
        std::memcpy(&d, &old, sizeof d);
        d += v;
        std::memcpy(&next, &d, sizeof next);
    } while (!__atomic_compare_exchange_n(p, &old, next, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
#endif
}


using namespace detail;

struct PadGroup {
    std::string ref;                       // REF.PIN
    std::map<int, PolygonSet> by_layer;    // pad shapes
    std::map<int, Point64> centre;         // per layer
    int net = 0;
};

PadGroup find_pad(const Board& b, const std::string& ref) {
    const auto dot = ref.find('.');
    if (dot == std::string::npos) throw std::runtime_error("loop element: pad '" + ref + "' is not REF.PIN");
    PadGroup g;
    g.ref = ref;
    const std::string comp = ref.substr(0, dot), pin = ref.substr(dot + 1);
    for (const auto& p : b.pads)
        if (p.component == comp && p.pin == pin) {
            g.by_layer[p.layer] = unite(g.by_layer[p.layer], p.shape);
            g.centre[p.layer] = p.at;
            g.net = p.net;
        }
    if (g.by_layer.empty()) throw std::runtime_error("loop element: no pad " + ref);
    return g;
}

// Straight current branch outside the sheets.
struct Filament {
    std::array<double, 3> p0, p1;  // m
    double radius = 0;
    double extra = 0;  // H
    // Where the current leaves the sheets into the filament and returns: triangles with
    // the share of the current each takes (by area, as a lead soldered over a pad).
    std::vector<std::pair<std::size_t, double>> end0, end1;
    bool vertical = false;
    int stitch = -1;              // index of the stitch it models
    double capacitance = 0, resistance = 0;  // a stitch's series capacitance (0: none) and resistance
};

struct Tri {
    std::array<Vec2, 3> v;  // m, counter-clockwise
    double z = 0, area = 0;
    int layer = 0;
    long node = -1;
    std::array<long, 3> basis{-1, -1, -1};  // RWG index of the edge opposite each vertex
    std::array<int, 3> sign{0, 0, 0};
    double size = 0;  // longest edge
    Vec2 c;           // centroid
    std::size_t job = 0, local = 0;  // island and index in its mesh
};

double dist_point_segment(double px, double py, double ax, double ay, double bx, double by) {
    const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
    double t = l2 > 0 ? ((px - ax) * ex + (py - ay) * ey) / l2 : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(px - ax - t * ex, py - ay - t * ey);
}

struct Discretisation {
    std::vector<Tri> tris;
    std::vector<Filament> fil;
    std::size_t nrwg = 0;
    std::size_t nnodes = 0;
    std::vector<std::pair<int, std::vector<std::size_t>>> sheet_tris;  // layer -> triangles (for output)
    std::vector<TriMesh> meshes;
    std::vector<int> mesh_layer;
};

class LoopSolver {
public:
    LoopSolver(const Board& b, const LoopSpec& loop, const LoopSettings& s) : b_(b), loop_(loop), s_(s) {}
    // With `sweep`, the frequency sweep runs on the last pass; otherwise it can follow on
    // the last mesh through sweep().
    LoopResult run(bool sweep = true);
    std::vector<LoopImpedance> sweep();
    LoopGeometry geometry();

private:
    const Board& b_;
    const LoopSpec& loop_;
    const LoopSettings s_;
    Discretisation last_;
    double tol_ = 1000.0;
    std::vector<PadGroup> pads_;  // per element: from, to
    std::set<int> loop_nets_, stitch_nets_;
    PolygonSet window_;
    std::vector<double> zmid_;  // m, per copper layer
    std::vector<std::array<double, 4>> path_;  // segments x0 y0 x1 y1 (nm) of the loop
    double h0_ = 0, hmax_ = 0;
    struct Ring {
        std::size_t drill;
        std::vector<int> layers;
        PolygonSet outline;
    };
    std::vector<Ring> rings_;
    std::map<int, std::vector<Polygon>> islands_;  // layer -> islands
    std::map<int, PolygonSet> constraints_;

    void setup();
    bool closes() const;
    // An element as a straight conductor: the layer it lies on, the layers of its two pads,
    // its ends (the pad centres) and its height (m).
    struct ElementPlace {
        int layer, la, lc;
        Point64 a, c;
        double z;
    };
    ElementPlace place_element(std::size_t k) const;
    using Seeds = std::vector<std::vector<std::array<double, 3>>>;  // per island
    Seeds initial_seeds(double scale) const;
    Discretisation discretise(const Seeds& seeds, double scale);
    double solve(const Discretisation& d, LoopResult* out, std::vector<std::array<double, 2>>* sheet, bool sweep = false);
    void filament_blocks(const Discretisation& d, std::vector<double>& F, std::vector<double>& FF) const;
    using Column = std::vector<std::pair<std::size_t, double>>;
    // M^T L M, -M^T L x_p and x_p^T L x_p without the dense L: its rows by blocks of
    // triangles, each multiplied out and dropped.
    void loop_system_streamed(const Discretisation& d, const std::vector<Column>& M, const Column& xp, std::vector<double>& Z,
                              std::vector<double>& rhs, double& xlx) const;
    // Resistance and inductance at the settings' frequencies: the loops' impedance with the
    // sheets' surface impedance and the holes' walls, from the loop inductance matrix ML
    // (M^T L M) and -M^T L x_p.
    void frequency_sweep(const Discretisation& d, const std::vector<std::vector<std::pair<std::size_t, double>>>& M,
                         const std::vector<std::pair<std::size_t, double>>& xp, const std::vector<double>& ML, const std::vector<double>& minus_mlxp,
                         double xlx, std::vector<LoopImpedance>& out) const;
    bool solve_by_loops(const Discretisation& d, const std::vector<double>* L, const std::vector<std::array<long, 2>>& ends, std::size_t src,
                        std::vector<LoopImpedance>* sweep,
                        std::vector<double>& x, double& inductance) const;
    std::vector<std::pair<int, const Polygon*>> jobs_;  // islands to mesh
    double size_at(double x, double y) const {
        double dmin = INFINITY;
        for (const auto& sg : path_) dmin = std::min(dmin, dist_point_segment(x, y, sg[0], sg[1], sg[2], sg[3]));
        return std::min(hmax_, h0_ + 0.35 * dmin);
    }
};

void LoopSolver::setup() {
    tol_ = b_.source.arc_tolerance_nm > 0 ? b_.source.arc_tolerance_nm : 1000.0;
    if (loop_.elements.size() < 1) throw std::runtime_error("loop " + loop_.name + " has no element");
    const int n = static_cast<int>(b_.copper.size());
    zmid_.assign(n + 1, 0.0);
    const auto zs = stackup_z(b_);
    for (std::size_t i = 0; i < b_.stackup.size(); ++i)
        if (b_.stackup[i].type == LayerType::Copper) zmid_[b_.stackup[i].copper_index] = 0.5 * (zs[i].top_nm + zs[i].bottom_nm) * 1e-9;
    for (const auto& e : loop_.elements) {
        pads_.push_back(find_pad(b_, e.from_pad));
        pads_.push_back(find_pad(b_, e.to_pad));
    }
    const std::size_t ne = loop_.elements.size();
    for (std::size_t k = 0; k < ne; ++k) {
        const auto& to = pads_[2 * k + 1];
        const auto& next = pads_[2 * ((k + 1) % ne)];
        if (to.net != next.net || to.net == 0)
            throw std::runtime_error("loop " + loop_.name + " is not closed: " + to.ref + " and " + next.ref + " are not on one net");
        loop_nets_.insert(to.net);
    }
    // Stitches: their pads follow the elements', and their nets' copper is modelled.
    for (const auto& st : loop_.stitches) {
        pads_.push_back(find_pad(b_, st.from_pad));
        pads_.push_back(find_pad(b_, st.to_pad));
        stitch_nets_.insert(pads_[pads_.size() - 2].net);
        stitch_nets_.insert(pads_.back().net);
    }
    // Window: the pads' bounding box with the margin, within the board.
    double x0 = INFINITY, y0 = INFINITY, x1 = -INFINITY, y1 = -INFINITY;
    double minpad = INFINITY;
    for (const auto& g : pads_)
        for (const auto& [layer, poly] : g.by_layer) {
            const Rect64 r = bounds(poly);
            x0 = std::min(x0, double(r.left));
            y0 = std::min(y0, double(r.top));
            x1 = std::max(x1, double(r.right));
            y1 = std::max(y1, double(r.bottom));
            minpad = std::min({minpad, double(r.right - r.left), double(r.bottom - r.top)});
        }
    h0_ = s_.min_size_nm > 0 ? s_.min_size_nm : std::clamp(0.5 * minpad, 1.0e5, 5.0e5);
    // Corridor of the loop current: the DC current of 1 A through each copper connection
    // of the loop (from one element to the next). Where it flows, the high-frequency current
    // flows too, if less spread out; the corridor sets the window and the element size.
    std::vector<std::array<double, 2>> corridor;
    {
        DcSettings dcs;
        dcs.target_error = 0.3;
        dcs.max_passes = 2;
        dcs.max_triangles = 400000;
        for (std::size_t k = 0; k < ne; ++k) {
            const auto& to = pads_[2 * k + 1];
            const auto& next = pads_[2 * ((k + 1) % ne)];
            if (to.ref == next.ref) continue;
            PortSpec a, c;
            a.name = "A" + std::to_string(k);
            a.pads = {to.ref};
            a.current_a = 1.0;
            c.name = "B" + std::to_string(k);
            c.pads = {next.ref};
            c.current_a = -1.0;
            dcs.ports.push_back(a);
            dcs.ports.push_back(c);
        }
        if (!dcs.ports.empty()) {
            const DcResult dc = solve_dc(b_, dcs);
            double jmax = 0;
            for (const auto& f : dc.layers)
                for (std::size_t t = 0; t < f.jx.size(); ++t) jmax = std::max(jmax, std::hypot(f.jx[t], f.jy[t]));
            std::set<std::pair<long long, long long>> taken;  // one point per cell of h0
            for (const auto& f : dc.layers)
                for (std::size_t t = 0; t < f.jx.size(); ++t) {
                    if (std::hypot(f.jx[t], f.jy[t]) < 0.05 * jmax) continue;
                    const auto& tr = f.triangles[t];
                    const double cx = (f.x[tr[0]] + f.x[tr[1]] + f.x[tr[2]]) / 3, cy = (f.y[tr[0]] + f.y[tr[1]] + f.y[tr[2]]) / 3;
                    if (taken.insert({std::llround(cx / h0_), std::llround(cy / h0_)}).second) corridor.push_back({cx, cy});
                }
        }
    }
    for (const auto& c : corridor) {
        x0 = std::min(x0, c[0]);
        y0 = std::min(y0, c[1]);
        x1 = std::max(x1, c[0]);
        y1 = std::max(y1, c[1]);
    }
    // The window: the margin around the pads and the corridor; grown if the copper inside it
    // does not join each element to the next.
    auto make_window = [&](double m) {
        const Path64 rect{Point64(std::llround(x0 - m), std::llround(y0 - m)), Point64(std::llround(x1 + m), std::llround(y0 - m)),
                          Point64(std::llround(x1 + m), std::llround(y1 + m)), Point64(std::llround(x0 - m), std::llround(y1 + m))};
        return intersection(union_polygons({rect}), b_.outline);
    };
    const Rect64 bb = bounds(b_.outline);
    const double board_size = std::max(double(bb.right - bb.left), double(bb.bottom - bb.top));
    double m = s_.window_margin_nm;
    for (int attempt = 0;; ++attempt) {
        window_ = make_window(m);
        if (closes() || m > 2 * board_size) break;
        m = m + std::max(std::max(x1 - x0, y1 - y0), s_.window_margin_nm) * (1 << attempt);
    }
    if (!closes()) throw std::runtime_error("loop " + loop_.name + " is open: no copper path joins its elements");
    const Rect64 wr = bounds(window_);
    const double wsize = std::max(double(wr.right - wr.left), double(wr.bottom - wr.top));
    hmax_ = wsize / 6.0;
    // Path of the loop: through each element, from each element to the next, and the
    // corridor's points.
    auto centre = [&](const PadGroup& g) {
        const auto& c = g.centre.begin()->second;
        return std::array<double, 2>{double(c.x), double(c.y)};
    };
    for (std::size_t k = 0; k < ne; ++k) {
        const auto a = centre(pads_[2 * k]), bb2 = centre(pads_[2 * k + 1]), c = centre(pads_[2 * ((k + 1) % ne)]);
        path_.push_back({a[0], a[1], bb2[0], bb2[1]});
        path_.push_back({bb2[0], bb2[1], c[0], c[1]});
    }
    for (const auto& c : corridor) path_.push_back({c[0], c[1], c[0], c[1]});
    for (std::size_t j = 0; j < loop_.stitches.size(); ++j) {
        const auto a = centre(pads_[2 * (ne + j)]), c = centre(pads_[2 * (ne + j) + 1]);
        path_.push_back({a[0], a[1], c[0], c[1]});
    }
    // Pad shapes are constraints and tie their triangles.
    for (const auto& g : pads_)
        for (const auto& [layer, poly] : g.by_layer) {
            auto& c = constraints_[layer];
            c.insert(c.end(), poly.begin(), poly.end());
        }
    // Copper islands in the window, on the loop's layers and the next depth_layers ones.
    std::set<int> element_layers;
    for (const auto& g : pads_)
        for (const auto& [layer, poly] : g.by_layer) element_layers.insert(layer);
    auto near_layer = [&](int l) {
        for (int e : element_layers)
            if (std::abs(l - e) <= s_.depth_layers) return true;
        return false;
    };
    for (const auto& cl : b_.copper)
        if (near_layer(cl.index))
            for (const auto& nc : cl.nets) {
                const PolygonSet part = intersection(nc.polygons, window_);
                for (const auto& island : part) {
                    const double area_mm2 = area_nm2(island) * 1e-12;
                    if (loop_nets_.count(nc.net) || stitch_nets_.count(nc.net) || area_mm2 >= s_.passive_min_area_mm2)
                        islands_[cl.index].push_back(island);
                }
            }
    // Plated holes in the window that join copper on several layers, near the loop.
    for (std::size_t d = 0; d < b_.drills.size(); ++d) {
        const auto& dr = b_.drills[d];
        if (!dr.plated || !contains(window_, dr.at)) continue;
        double dmin = INFINITY;
        for (const auto& sg : path_) dmin = std::min(dmin, dist_point_segment(double(dr.at.x), double(dr.at.y), sg[0], sg[1], sg[2], sg[3]));
        if (dmin > s_.via_radius_nm) continue;
        Ring r;
        r.drill = d;
        for (int l = std::clamp(dr.from_layer, 1, n); l <= std::clamp(dr.to_layer, 1, n); ++l) {
            auto it = islands_.find(l);
            if (it == islands_.end()) continue;
            for (const auto& isl : it->second)
                if (contains(isl, dr.at)) {
                    r.layers.push_back(l);
                    break;
                }
        }
        if (r.layers.size() < 2) continue;
        // The wall as an octagon of equal area: the current enters around it.
        r.outline = union_polygons({dr.slot_end ? hole_outline(dr, tol_)
                                                : circle_path(dr.at, 0.5 * dr.diameter_nm, 0.08 * dr.diameter_nm)});
        for (int l : r.layers) {
            auto& c = constraints_[l];
            c.insert(c.end(), r.outline.begin(), r.outline.end());
        }
        rings_.push_back(std::move(r));
    }
    // Away from the loop the copper is simplified to the local element size: holes smaller
    // than half of it (antipads in planes) are closed, other outlines simplified to 0.15 of
    // it. Near the loop the size is small and the outlines stay as drawn.
    for (auto& [layer, list] : islands_) {
        std::vector<Polygon> simplified;
        for (const auto& isl : list) {
            Paths64 paths;
            auto local = [&](const Path64& r) {
                const Rect64 bb = bounds(Polygon{r, {}});
                double dmin = INFINITY;
                const double cx = 0.5 * (bb.left + bb.right), cy = 0.5 * (bb.top + bb.bottom);
                for (const auto& sg : path_) dmin = std::min(dmin, dist_point_segment(cx, cy, sg[0], sg[1], sg[2], sg[3]));
                const double extent = 0.5 * std::hypot(double(bb.right - bb.left), double(bb.bottom - bb.top));
                return std::make_pair(std::min(hmax_, h0_ + 0.35 * std::max(0.0, dmin - extent)), extent);
            };
            const auto [ho, eo] = local(isl.outer);
            paths.push_back(Clipper2Lib::SimplifyPath(isl.outer, std::max(1000.0, 0.15 * ho), true));
            for (const auto& hole : isl.holes) {
                const auto [hh, eh] = local(hole);
                if (2 * eh < 0.5 * hh) continue;
                paths.push_back(Clipper2Lib::SimplifyPath(hole, std::max(1000.0, 0.15 * hh), true));
            }
            for (auto& p : union_polygons(paths)) simplified.push_back(std::move(p));
        }
        list = std::move(simplified);
    }
    for (const auto& [layer, list] : islands_)
        for (const auto& isl : list) jobs_.emplace_back(layer, &isl);
}

bool LoopSolver::closes() const {
    // Islands of the loop's nets in the window, joined by plated holes; every connection
    // from an element's end to the next element's start must lie within one group.
    struct Island {
        int layer;
        PolygonSet shape;
    };
    std::vector<Island> isl;
    for (const auto& cl : b_.copper)
        for (const auto& nc : cl.nets) {
            if (!loop_nets_.count(nc.net)) continue;
            for (const auto& p : intersection(nc.polygons, window_)) isl.push_back({cl.index, PolygonSet{p}});
        }
    UnionFind uf(isl.size());
    for (const auto& dr : b_.drills) {
        if (!dr.plated || !contains(window_, dr.at)) continue;
        long first = -1;
        for (std::size_t i = 0; i < isl.size(); ++i)
            if (isl[i].layer >= dr.from_layer && isl[i].layer <= dr.to_layer && contains(isl[i].shape, dr.at)) {
                if (first < 0) first = static_cast<long>(i);
                else uf.unite(i, static_cast<std::size_t>(first));
            }
    }
    auto group = [&](const PadGroup& g) -> long {
        for (const auto& [layer, c] : g.centre)
            for (std::size_t i = 0; i < isl.size(); ++i)
                if (isl[i].layer == layer && contains(isl[i].shape, c)) return static_cast<long>(uf.find(i));
        return -1;
    };
    const std::size_t ne = loop_.elements.size();
    for (std::size_t k = 0; k < ne; ++k) {
        const long a = group(pads_[2 * k + 1]), c = group(pads_[2 * ((k + 1) % ne)]);
        if (a < 0 || a != c) return false;
    }
    return true;
}

LoopSolver::Seeds LoopSolver::initial_seeds(double scale) const {
    // Quadtree leaves over the window, with the element size growing away from the loop.
    const double h0 = h0_ * scale, hmax = hmax_ * scale;
    const Rect64 wb = bounds(window_);
    std::vector<std::array<double, 3>> seeds;
    struct Cell {
        double x, y, s;
    };
    const double side = std::max(double(wb.right - wb.left), double(wb.bottom - wb.top));
    std::vector<Cell> stack{{double(wb.left), double(wb.top), side}};
    while (!stack.empty()) {
        const Cell c = stack.back();
        stack.pop_back();
        const double cx = c.x + c.s / 2, cy = c.y + c.s / 2;
        const double h = std::min(hmax, std::max(h0, (size_at(cx, cy) - h0_) * scale + h0));
        if (c.s > h && c.s > 0.5 * h0) {
            const double hs = c.s / 2;
            stack.push_back({c.x, c.y, hs});
            stack.push_back({c.x + hs, c.y, hs});
            stack.push_back({c.x, c.y + hs, hs});
            stack.push_back({c.x + hs, c.y + hs, hs});
        } else {
            seeds.push_back({cx, cy, c.s});
        }
    }
    Seeds out(jobs_.size());
    parallel_for(jobs_.size(), [&](std::size_t j) {
        const Polygon& isl = *jobs_[j].second;
        const PolygonSet one{isl};
        const Rect64 r = bounds(one);
        for (const auto& p : seeds)
            if (p[0] >= r.left && p[0] <= r.right && p[1] >= r.top && p[1] <= r.bottom &&
                contains(one, Point64(std::llround(p[0]), std::llround(p[1]))))
                out[j].push_back(p);
        // Rows of vertices along the copper edges near the loop, at offsets growing by
        // three from 20 um: in the high-frequency limit the sheet current grows as
        // d^-1/2 towards an edge, which elements long along the edge and thin across it
        // resolve at a small cost.
        auto ring = [&](const Path64& rg) {
            // Points every local element size along the boundary (by arc length, whatever
            // the segments), each row offset along the local inward normal.
            double carry = 0;  // arc length since the last point
            for (std::size_t k = 0; k < rg.size(); ++k) {
                const auto& a = rg[k];
                const auto& c = rg[(k + 1) % rg.size()];
                const double ex = double(c.x - a.x), ey = double(c.y - a.y), len = std::hypot(ex, ey);
                if (len <= 0) continue;
                const double mx = 0.5 * (a.x + c.x), my = 0.5 * (a.y + c.y);
                const double h = std::min(hmax, std::max(h0, (size_at(mx, my) - h0_) * scale + h0));
                if (h > 0.5 * hmax) {
                    carry = 0;
                    continue;
                }
                double nx = -ey / len, ny = ex / len;
                const double probe = 1000.0;
                if (!contains(isl, Point64(std::llround(mx + probe * nx), std::llround(my + probe * ny)))) {
                    nx = -nx;
                    ny = -ny;
                }
                for (double at = h - carry; at < len; at += h) {
                    const double t = at / len;
                    for (double off = 20000.0; off < 0.5 * h; off *= 3.0) {
                        const double px = a.x + t * ex + off * nx, py = a.y + t * ey + off * ny;
                        if (contains(isl, Point64(std::llround(px), std::llround(py)))) out[j].push_back({px, py, off});
                    }
                }
                carry = std::fmod(carry + len, h);
            }
        };
        ring(isl.outer);
        for (const auto& hole : isl.holes) ring(hole);
    });
    return out;
}

LoopSolver::ElementPlace LoopSolver::place_element(std::size_t k) const {
    const int n = static_cast<int>(b_.copper.size());
    const PadGroup &a = pads_[2 * k], &c = pads_[2 * k + 1];
    ElementPlace p;
    // The element lies on the outer layer of its pads (the side of the component).
    p.layer = a.by_layer.begin()->first;
    for (const auto& [l, poly] : a.by_layer)
        if (c.by_layer.count(l) && (l == 1 || l == n)) p.layer = l;
    p.la = a.by_layer.count(p.layer) ? p.layer : a.by_layer.begin()->first;
    p.lc = c.by_layer.count(p.layer) ? p.layer : c.by_layer.begin()->first;
    p.a = a.centre.at(p.la);
    p.c = c.centre.at(p.lc);
    const std::size_t ne = loop_.elements.size();
    const double r = (k < ne ? loop_.elements[k].radius_nm : loop_.stitches[k - ne].radius_nm) * 1e-9;
    p.z = zmid_[p.layer] + (p.layer == n && n > 1 ? -r : r);
    return p;
}

LoopGeometry LoopSolver::geometry() {
    setup();
    LoopGeometry g;
    g.window = window_;
    g.copper = islands_;
    g.z_m = zmid_;
    for (const auto& rg : rings_) {
        const auto& dr = b_.drills[rg.drill];
        g.vias.push_back({dr.at, 0.5 * dr.diameter_nm, rg.layers.front(), rg.layers.back()});
    }
    const std::size_t ne = loop_.elements.size();
    for (std::size_t k = 0; k < ne + loop_.stitches.size(); ++k) {
        const ElementPlace ep = place_element(k);
        const bool stitch = k >= ne;
        g.elements.push_back({pads_[2 * k].ref, pads_[2 * k + 1].ref, ep.a, ep.c, ep.layer, ep.z,
                              (stitch ? loop_.stitches[k - ne].radius_nm : loop_.elements[k].radius_nm) * 1e-9, stitch,
                              stitch ? loop_.stitches[k - ne].capacitance_f : 0.0});
    }
    return g;
}

Discretisation LoopSolver::discretise(const Seeds& seeds, double scale) {
    Discretisation D;
    const double h0 = h0_ * scale, hmax = hmax_ * scale;
    MeshSettings ms;
    ms.boundary_tolerance_nm = std::min(5000.0, 0.05 * h0);
    ms.max_edge_nm = hmax;
    ms.grading = 0.0;
    // No angle criterion: the rows along the edges are thin on purpose, and the seeds keep
    // the interior well shaped.
    ms.min_angle_deg = 0.0;
    ms.min_edge_nm = 0.02 * h0;
    const auto& jobs = jobs_;
    std::vector<TriMesh> meshes(jobs.size());
    parallel_for(jobs.size(), [&](std::size_t j) {
        const auto it = constraints_.find(jobs[j].first);
        meshes[j] = mesh_region(*jobs[j].second, it == constraints_.end() ? PolygonSet{} : it->second, ms, seeds[j]);
    });

    // Triangles and RWG edges; every triangle is a node of the current law.
    for (std::size_t j = 0; j < jobs.size(); ++j) {
        const int layer = jobs[j].first;
        const TriMesh& tm = meshes[j];
        const std::size_t base = D.tris.size();
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::pair<std::size_t, int>> edge_of;
        for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
            Tri T;
            const auto& tr = tm.triangles[t];
            for (int k = 0; k < 3; ++k) T.v[k] = {tm.x[tr[k]] * 1e-9, tm.y[tr[k]] * 1e-9};
            T.z = zmid_[layer];
            T.layer = layer;
            T.area = triangle_area(tm, t) * 1e-18;
            T.c = {(T.v[0].x + T.v[1].x + T.v[2].x) / 3, (T.v[0].y + T.v[1].y + T.v[2].y) / 3};
            for (int k = 0; k < 3; ++k)
                T.size = std::max(T.size, std::hypot(T.v[k].x - T.v[(k + 1) % 3].x, T.v[k].y - T.v[(k + 1) % 3].y));
            T.node = static_cast<long>(base + t);
            T.job = j;
            T.local = t;
            D.tris.push_back(T);
            for (int k = 0; k < 3; ++k) {
                const auto a = tr[(k + 1) % 3], c = tr[(k + 2) % 3];
                const auto key = std::make_pair(std::min(a, c), std::max(a, c));
                auto it = edge_of.find(key);
                if (it == edge_of.end()) {
                    edge_of[key] = {base + t, k};
                } else {
                    // Interior edge: the first triangle is T+, this one T-.
                    const long idx = static_cast<long>(D.nrwg++);
                    D.tris[it->second.first].basis[it->second.second] = idx;
                    D.tris[it->second.first].sign[it->second.second] = +1;
                    D.tris[base + t].basis[k] = idx;
                    D.tris[base + t].sign[k] = -1;
                }
            }
        }
        D.meshes.push_back(tm);
        D.mesh_layer.push_back(layer);
    }
    // Triangles of a layer inside a region, with their shares of the area.
    auto region_tris = [&](int layer, const PolygonSet& region) {
        std::vector<std::pair<std::size_t, double>> out;
        const Rect64 r = bounds(region);
        double total = 0;
        for (std::size_t t = 0; t < D.tris.size(); ++t) {
            const Tri& T = D.tris[t];
            if (T.layer != layer) continue;
            const Point64 cp(std::llround(T.c.x * 1e9), std::llround(T.c.y * 1e9));
            if (cp.x < r.left || cp.x > r.right || cp.y < r.top || cp.y > r.bottom || !contains(region, cp)) continue;
            out.push_back({t, T.area});
            total += T.area;
        }
        for (auto& e : out) e.second /= total;
        return out;
    };
    // Elements.
    for (std::size_t k = 0; k < loop_.elements.size(); ++k) {
        const auto& e = loop_.elements[k];
        const PadGroup &a = pads_[2 * k], &c = pads_[2 * k + 1];
        const ElementPlace ep = place_element(k);
        Filament f;
        f.end0 = region_tris(ep.la, a.by_layer.at(ep.la));
        f.end1 = region_tris(ep.lc, c.by_layer.at(ep.lc));
        if (f.end0.empty() || f.end1.empty())
            throw std::runtime_error("loop " + loop_.name + ": pad " + (f.end0.empty() ? a.ref : c.ref) + " has no copper");
        f.p0 = {ep.a.x * 1e-9, ep.a.y * 1e-9, ep.z};
        f.p1 = {ep.c.x * 1e-9, ep.c.y * 1e-9, ep.z};
        f.radius = e.radius_nm * 1e-9;
        f.extra = e.extra_inductance_nh * 1e-9;
        D.fil.push_back(f);
    }
    // Stitches.
    for (std::size_t j = 0; j < loop_.stitches.size(); ++j) {
        const auto& st = loop_.stitches[j];
        const std::size_t k = loop_.elements.size() + j;
        const PadGroup &a = pads_[2 * k], &c = pads_[2 * k + 1];
        const ElementPlace ep = place_element(k);
        Filament f;
        f.end0 = region_tris(ep.la, a.by_layer.at(ep.la));
        f.end1 = region_tris(ep.lc, c.by_layer.at(ep.lc));
        if (f.end0.empty() || f.end1.empty())
            throw std::runtime_error("loop " + loop_.name + ": pad " + (f.end0.empty() ? a.ref : c.ref) + " of a stitch has no copper in the window");
        f.p0 = {ep.a.x * 1e-9, ep.a.y * 1e-9, ep.z};
        f.p1 = {ep.c.x * 1e-9, ep.c.y * 1e-9, ep.z};
        f.radius = st.radius_nm * 1e-9;
        f.extra = st.extra_inductance_nh * 1e-9;
        f.stitch = static_cast<int>(j);
        f.capacitance = st.capacitance_f;
        f.resistance = st.resistance_ohm;
        D.fil.push_back(f);
    }
    // Barrels.
    for (std::size_t r = 0; r < rings_.size(); ++r) {
        const auto& rg = rings_[r];
        const auto& dr = b_.drills[rg.drill];
        for (std::size_t k = 1; k < rg.layers.size(); ++k) {
            Filament f;
            f.end0 = region_tris(rg.layers[k - 1], rg.outline);
            f.end1 = region_tris(rg.layers[k], rg.outline);
            if (f.end0.empty() || f.end1.empty()) continue;
            f.p0 = {dr.at.x * 1e-9, dr.at.y * 1e-9, zmid_[rg.layers[k - 1]]};
            f.p1 = {dr.at.x * 1e-9, dr.at.y * 1e-9, zmid_[rg.layers[k]]};
            f.radius = 0.5 * dr.diameter_nm * 1e-9;
            f.vertical = true;
            D.fil.push_back(f);
        }
    }
    D.nnodes = D.tris.size();

    return D;
}

// Kernel of a triangle pair: K[i][j] = integral over T and T' of (rho - p_i).(rho' - p'_j) / R.
std::array<std::array<double, 3>, 3> pair_kernel(const Tri& T, const Tri& U) {
    const double dx = T.c.x - U.c.x, dy = T.c.y - U.c.y, dz = T.z - U.z;
    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    const TriangleRule& rule = dist < 2.5 * std::max(T.size, U.size) ? triangle_rule_28() : triangle_rule_7();
    std::array<std::array<double, 3>, 3> K{};
    for (std::size_t q = 0; q < rule.bary.size(); ++q) {
        const auto& bq = rule.bary[q];
        const Vec2 r{bq[0] * T.v[0].x + bq[1] * T.v[1].x + bq[2] * T.v[2].x, bq[0] * T.v[0].y + bq[1] * T.v[1].y + bq[2] * T.v[2].y};
        const auto pot = triangle_potential(U.v, r, T.z - U.z);
        const double w = rule.weight[q] * T.area;
        for (int j = 0; j < 3; ++j) {
            // integral over U of (rho' - p'_j)/R = I1 + (r - p'_j) I0
            const double gx = pot.vector.x + (r.x - U.v[j].x) * pot.scalar;
            const double gy = pot.vector.y + (r.y - U.v[j].y) * pot.scalar;
            for (int i = 0; i < 3; ++i) K[i][j] += w * ((r.x - T.v[i].x) * gx + (r.y - T.v[i].y) * gy);
        }
    }
    return K;
}

// The loop currents of a basis that satisfies the current law by construction: a spanning
// tree of the RWG branches within each connected sheet, and one of the filaments between
// the sheets. The source's ampere returns through the trees (x_p); every other branch with
// the tree path that closes it is a loop (a column of M). The branch currents
// x = x_p + M y minimise the energy x^T L x over y: (M^T L M) y = -M^T L x_p, a dense
// system of the loops alone, about a third of the branches, instead of the branches and
// the current law together. Returns false where a filament's end spans several sheets.
bool LoopSolver::solve_by_loops(const Discretisation& D, const std::vector<double>* Lp, const std::vector<std::array<long, 2>>& ends,
                                std::size_t src, std::vector<LoopImpedance>* sweep, std::vector<double>& x, double& inductance) const {
    const std::size_t nr = D.nrwg, nf = D.fil.size(), nb = nr + nf, nn = D.nnodes;
    // RWG forest, breadth first: parents before children in `order`.
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> adj(nn);
    for (std::size_t k = 0; k < nr; ++k) {
        const long a = ends[k][0], c = ends[k][1];
        if (a < 0 || c < 0 || a == c) continue;
        adj[static_cast<std::size_t>(a)].push_back({static_cast<std::size_t>(c), k});
        adj[static_cast<std::size_t>(c)].push_back({static_cast<std::size_t>(a), k});
    }
    std::vector<long> parent(nn, -1), pedge(nn, -1);
    std::vector<int> comp(nn, -1);
    std::vector<std::size_t> order;
    order.reserve(nn);
    std::vector<char> tree(nb, 0);
    int ncomp = 0;
    for (std::size_t s0 = 0; s0 < nn; ++s0) {
        if (comp[s0] >= 0) continue;
        comp[s0] = ncomp;
        std::size_t head = order.size();
        order.push_back(s0);
        while (head < order.size()) {
            const std::size_t v = order[head++];
            for (const auto& [u, k] : adj[v])
                if (comp[u] < 0) {
                    comp[u] = ncomp;
                    parent[u] = static_cast<long>(v);
                    pedge[u] = static_cast<long>(k);
                    tree[k] = 1;
                    order.push_back(u);
                }
        }
        ++ncomp;
    }
    // Filaments between the sheets: each end within one sheet.
    std::vector<std::array<int, 2>> fc(nf);
    std::vector<std::array<double, 2>> fw(nf);
    for (std::size_t q = 0; q < nf; ++q) {
        const Filament& f = D.fil[q];
        for (int e = 0; e < 2; ++e) {
            const auto& list = e == 0 ? f.end0 : f.end1;
            int c = -1;
            double w = 0;
            for (const auto& [t, wt] : list) {
                if (c < 0) c = comp[t];
                else if (comp[t] != c) return false;
                w += wt;
            }
            if (c < 0 || !(w > 0)) return false;
            fc[q][static_cast<std::size_t>(e)] = c;
            fw[q][static_cast<std::size_t>(e)] = w;
        }
    }
    std::vector<int> cparent(static_cast<std::size_t>(ncomp), -1), corder;
    std::vector<long> cedge(static_cast<std::size_t>(ncomp), -1);
    {
        std::vector<std::vector<std::pair<int, std::size_t>>> cadj(static_cast<std::size_t>(ncomp));
        for (std::size_t q = 0; q < nf; ++q) {
            if (nr + q == src || fc[q][0] == fc[q][1]) continue;
            cadj[static_cast<std::size_t>(fc[q][0])].push_back({fc[q][1], q});
            cadj[static_cast<std::size_t>(fc[q][1])].push_back({fc[q][0], q});
        }
        std::vector<char> seen(static_cast<std::size_t>(ncomp), 0);
        for (int c0 = 0; c0 < ncomp; ++c0) {
            if (seen[static_cast<std::size_t>(c0)]) continue;
            seen[static_cast<std::size_t>(c0)] = 1;
            std::size_t head = corder.size();
            corder.push_back(c0);
            while (head < corder.size()) {
                const int c = corder[head++];
                for (const auto& [u, q] : cadj[static_cast<std::size_t>(c)])
                    if (!seen[static_cast<std::size_t>(u)]) {
                        seen[static_cast<std::size_t>(u)] = 1;
                        cparent[static_cast<std::size_t>(u)] = c;
                        cedge[static_cast<std::size_t>(u)] = static_cast<long>(q);
                        tree[nr + q] = 1;
                        corder.push_back(u);
                    }
            }
        }
    }
    // Branch currents on the trees whose current-law contributions meet r (consumed):
    // filaments between the sheets first, then the RWG branches, leaves first.
    using Column = std::vector<std::pair<std::size_t, double>>;
    auto close = [&](std::vector<double>& r, Column& col) {
        std::vector<double> rc(static_cast<std::size_t>(ncomp), 0.0);
        for (std::size_t v = 0; v < nn; ++v) rc[static_cast<std::size_t>(comp[v])] += r[v];
        for (auto it = corder.rbegin(); it != corder.rend(); ++it) {
            const int c = *it;
            if (cedge[static_cast<std::size_t>(c)] < 0) continue;
            const std::size_t q = static_cast<std::size_t>(cedge[static_cast<std::size_t>(c)]);
            const int pc = cparent[static_cast<std::size_t>(c)];
            // The filament adds +w0 I to the sheet of its first end, -w1 I to the other.
            const double I = fc[q][0] == c ? rc[static_cast<std::size_t>(c)] / fw[q][0] : -rc[static_cast<std::size_t>(c)] / fw[q][1];
            rc[static_cast<std::size_t>(c)] = 0;
            rc[static_cast<std::size_t>(pc)] -= fc[q][0] == pc ? fw[q][0] * I : -fw[q][1] * I;
            if (I == 0) continue;
            col.push_back({nr + q, I});
            for (const auto& [t, w] : D.fil[q].end0) r[t] -= w * I;
            for (const auto& [t, w] : D.fil[q].end1) r[t] += w * I;
        }
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            const std::size_t v = *it;
            if (pedge[v] < 0) continue;
            const auto k = static_cast<std::size_t>(pedge[v]);
            const double I = ends[k][0] == static_cast<long>(v) ? r[v] : -r[v];
            const auto u = static_cast<std::size_t>(parent[v]);
            r[u] -= ends[k][0] == static_cast<long>(u) ? I : -I;
            r[v] = 0;
            if (I != 0) col.push_back({k, I});
        }
        double left = 0;
        for (std::size_t v = 0; v < nn; ++v) left = std::max(left, std::abs(r[v]));
        return left;
    };
    auto required = [&](std::size_t k, std::vector<double>& r) {  // r = -(the branch's column)
        std::fill(r.begin(), r.end(), 0.0);
        if (k < nr) {
            if (ends[k][0] >= 0) r[static_cast<std::size_t>(ends[k][0])] -= 1.0;
            if (ends[k][1] >= 0) r[static_cast<std::size_t>(ends[k][1])] += 1.0;
            return;
        }
        for (const auto& [t, w] : D.fil[k - nr].end0) r[t] -= w;
        for (const auto& [t, w] : D.fil[k - nr].end1) r[t] += w;
    };
    // The source with its return path, and the loops.
    std::vector<double> r(nn);
    Column xp{{src, 1.0}};
    required(src, r);
    if (close(r, xp) > 1e-9) return false;
    std::vector<std::size_t> loop_branch;
    for (std::size_t k = 0; k < nb; ++k)
        if (k != src && !tree[k]) loop_branch.push_back(k);
    const std::size_t nl = loop_branch.size();
    std::vector<Column> M(nl);
    std::atomic<bool> bad{false};
    parallel_for(nl, [&](std::size_t j) {
        std::vector<double> rr(nn);
        M[j].push_back({loop_branch[j], 1.0});
        required(loop_branch[j], rr);
        if (close(rr, M[j]) > 1e-9) bad = true;
    });
    if (bad) return false;
    std::vector<double> Z(nl * nl), rhs(nl);
    double xLx = 0;
    if (!Lp) {
        loop_system_streamed(D, M, xp, Z, rhs, xLx);
    } else {
    const std::vector<double>& L = *Lp;
    // L M (branches x loops), M^T L M and M^T L x_p.
    std::vector<double> LM(nb * nl, 0.0), Lxp(nb, 0.0);
    parallel_for(nl, [&](std::size_t j) {
        for (const auto& [k, v] : M[j]) {
            const double* row = &L[k * nb];
            for (std::size_t i = 0; i < nb; ++i) LM[i * nl + j] += v * row[i];
        }
    });
    for (const auto& [k, v] : xp) {
        const double* row = &L[k * nb];
        for (std::size_t i = 0; i < nb; ++i) Lxp[i] += v * row[i];
    }
    parallel_for(nl, [&](std::size_t i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double t = 0;
            for (const auto& [k, v] : M[i]) t += v * LM[k * nl + j];
            Z[i * nl + j] = Z[j * nl + i] = t;
        }
        double t = 0;
        for (const auto& [k, v] : M[i]) t += v * Lxp[k];
        rhs[i] = -t;
    });
    std::vector<double>().swap(LM);
    for (const auto& [k, v] : xp) xLx += v * Lxp[k];
    }
    const std::vector<double> minus_c = rhs;
    if (sweep) frequency_sweep(D, M, xp, Z, rhs, xLx, *sweep);
    if (nl > 0) {
        if (!cholesky(Z, nl)) throw std::runtime_error("loop " + loop_.name + ": the loop inductance matrix is not positive definite");
        forward_solve(Z, nl, rhs, 1);
        backward_solve(Z, nl, rhs, 1);
    }
    std::vector<double> xf(nb, 0.0);
    for (const auto& [k, v] : xp) xf[k] += v;
    for (std::size_t j = 0; j < nl; ++j)
        for (const auto& [k, v] : M[j]) xf[k] += rhs[j] * v;
    // Inductance: the energy x^T L x of one ampere; without L, x_p^T L x_p + y^T M^T L x_p,
    // the same for the solution y of the loops.
    if (Lp) {
        const std::vector<double>& L = *Lp;
        std::vector<double> part(nb, 0.0);
        parallel_for(nb, [&](std::size_t i) {
            double t = 0;
            const double* row = &L[i * nb];
            for (std::size_t k = 0; k < nb; ++k) t += row[k] * xf[k];
            part[i] = xf[i] * t;
        });
        inductance = std::accumulate(part.begin(), part.end(), 0.0);
    } else {
        inductance = xLx;
        for (std::size_t j = 0; j < nl; ++j) inductance -= rhs[j] * minus_c[j];
    }
    for (std::size_t k = 0; k < nb; ++k)
        if (k != src) x[k < src ? k : k - 1] = xf[k];
    return true;
}

namespace {

using cplx = std::complex<double>;

// Solve A x = b for a complex symmetric matrix (A = A^T, not Hermitian) by LDL^T without
// pivoting; A of the form R + j omega L with R, L positive (semi)definite.
bool complex_symmetric_solve(std::vector<cplx>& a, std::size_t n, std::vector<cplx>& b) {
    std::vector<cplx> w(n);
    for (std::size_t j = 0; j < n; ++j) {
        // w_k = L_jk D_k: row j of L D, used by every row below.
        cplx d = a[j * n + j];
        for (std::size_t k = 0; k < j; ++k) {
            w[k] = a[j * n + k] * a[k * n + k];
            d -= a[j * n + k] * w[k];
        }
        if (std::abs(d) == 0.0) return false;
        a[j * n + j] = d;
        parallel_for(n - j - 1, [&](std::size_t q) {
            const std::size_t i = j + 1 + q;
            const cplx* ai = &a[i * n];
            cplx v = ai[j];
            for (std::size_t k = 0; k < j; ++k) v -= ai[k] * w[k];
            a[i * n + j] = v / d;
        });
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < i; ++k) b[i] -= a[i * n + k] * b[k];
    for (std::size_t i = 0; i < n; ++i) b[i] /= a[i * n + i];
    for (std::size_t i = n; i-- > 0;)
        for (std::size_t k = i + 1; k < n; ++k) b[i] -= a[k * n + i] * b[k];
    return true;
}

// Surface impedance of a copper sheet carrying current on both faces, per square:
// (1 + j) / (2 sigma delta) coth((1 + j) t / (2 delta)); 1 / (sigma t) at low frequency.
cplx sheet_impedance(double f, double sigma, double t) {
    const double delta = 1.0 / std::sqrt(kPi * f * 4e-7 * kPi * sigma);
    const cplx x = cplx(1, 1) * t / (2 * delta);
    if (std::abs(x) < 1e-3) return 1.0 / (sigma * t);
    const cplx e = std::exp(-2.0 * x);
    return cplx(1, 1) / (2 * sigma * delta) * (1.0 + e) / (1.0 - e);
}

// The same for one face of a wall of thickness t (a plated hole's barrel).
cplx wall_impedance(double f, double sigma, double t) {
    const double delta = 1.0 / std::sqrt(kPi * f * 4e-7 * kPi * sigma);
    const cplx x = cplx(1, 1) * t / delta;
    if (std::abs(x) < 1e-3) return 1.0 / (sigma * t);
    const cplx e = std::exp(-2.0 * x);
    return cplx(1, 1) / (sigma * delta) * (1.0 + e) / (1.0 - e);
}

}  // namespace

void LoopSolver::frequency_sweep(const Discretisation& D, const std::vector<std::vector<std::pair<std::size_t, double>>>& M,
                                 const std::vector<std::pair<std::size_t, double>>& xp, const std::vector<double>& ML,
                                 const std::vector<double>& minus_mlxp, double xlx, std::vector<LoopImpedance>& out) const {
    const std::size_t nl = M.size(), nr = D.nrwg;
    // Copper per layer: conductivity and thickness.
    std::map<int, std::pair<double, double>> metal;
    for (const auto& l : b_.stackup)
        if (l.type == LayerType::Copper)
            metal[l.copper_index] = {l.conductivity_s_per_m ? *l.conductivity_s_per_m : 5.8e7, l.thickness_nm > 0 ? l.thickness_nm * 1e-9 : 35e-6};
    // Branch-level resistance pattern: per layer, the RWG functions' overlap on each triangle,
    // int f_i . f_j dA (to be multiplied by the sheet impedance per square).
    // A zero-thickness sheet carries a current growing as d^-1/2 towards its edges, and the
    // resistance of that current grows without bound as the mesh refines; copper of
    // thickness t spreads it over about t / 2. A triangle nearer an edge than t / 2 counts
    // with the share d / (t / 2) of its overlap: the current density held at its value there.
    std::vector<std::vector<std::array<double, 4>>> edges(D.meshes.size());  // per mesh: boundary segments (m)
    for (std::size_t j = 0; j < D.meshes.size(); ++j) {
        const TriMesh& m = D.meshes[j];
        std::map<std::pair<std::uint32_t, std::uint32_t>, int> count;
        for (const auto& t : m.triangles)
            for (int k = 0; k < 3; ++k) ++count[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}];
        for (const auto& [e, c] : count)
            if (c == 1) edges[j].push_back({m.x[e.first] * 1e-9, m.y[e.first] * 1e-9, m.x[e.second] * 1e-9, m.y[e.second] * 1e-9});
    }
    // Per layer, the overlaps as computed (gram[0], uniform current at DC) and with the edge
    // weighting (gram[1], the skin effect's limit); the sweep blends them.
    std::map<int, std::map<std::pair<std::size_t, std::size_t>, double>> gram[2];
    for (const auto& T : D.tris) {
        const double A = T.area;
        double weight = 1.0;
        {
            const double half_t = 0.5 * (metal.count(T.layer) ? metal[T.layer].second : 35e-6);
            double d = INFINITY;
            for (const auto& e : edges[T.job]) d = std::min(d, dist_point_segment(T.c.x, T.c.y, e[0], e[1], e[2], e[3]));
            if (d < half_t) weight = d / half_t;
        }
        double r2 = 0;  // int |r|^2 dA / A
        for (int a = 0; a < 3; ++a)
            for (int c = a; c < 3; ++c) r2 += (T.v[a].x * T.v[c].x + T.v[a].y * T.v[c].y) / 6.0;
        for (int i = 0; i < 3; ++i) {
            if (T.basis[i] < 0) continue;
            for (int j = 0; j < 3; ++j) {
                if (T.basis[j] < 0) continue;
                const Vec2 vi = T.v[i], vj = T.v[j];
                const double integral = A * (r2 - ((vi.x + vj.x) * T.c.x + (vi.y + vj.y) * T.c.y) + vi.x * vj.x + vi.y * vj.y);
                const std::pair<std::size_t, std::size_t> ij{static_cast<std::size_t>(T.basis[i]), static_cast<std::size_t>(T.basis[j])};
                const double g = T.sign[i] * T.sign[j] / (4 * A * A) * integral;
                gram[0][T.layer][ij] += g;
                gram[1][T.layer][ij] += weight * g;
            }
        }
    }
    // Rows of M by branch: the loops through each branch.
    std::vector<std::vector<std::pair<std::size_t, double>>> rows(nr + D.fil.size());
    for (std::size_t j = 0; j < nl; ++j)
        for (const auto& [k, v] : M[j]) rows[k].push_back({j, v});
    std::vector<double> xpv(nr + D.fil.size(), 0.0);
    for (const auto& [k, v] : xp) xpv[k] += v;
    // Reduced per layer: G = M^T Gram M, g = M^T Gram x_p, s = x_p^T Gram x_p.
    struct Reduced {
        std::vector<double> G[2], g[2];
        double s[2] = {0, 0};
        double sigma = 5.8e7, t = 35e-6;
    };
    std::map<int, Reduced> red;
    const std::size_t nb = nr + D.fil.size();
    for (int w = 0; w < 2; ++w)
        for (const auto& [layer, entries] : gram[w]) {
            Reduced& r = red[layer];
            r.G[w].assign(nl * nl, 0.0);
            r.g[w].assign(nl, 0.0);
            if (metal.count(layer)) std::tie(r.sigma, r.t) = metal[layer];
            // G M (branches x loops), then M^T (G M): the loops through a branch near the
            // trees' roots are many, so the products go through the sparse G first.
            std::vector<std::vector<std::pair<std::size_t, double>>> grow(nb);
            for (const auto& [ij, v] : entries) grow[ij.first].push_back({ij.second, v});
            std::vector<double> GM(nb * nl, 0.0), Gx(nb, 0.0);
            parallel_for(nb, [&](std::size_t bi) {
                double* out = &GM[bi * nl];
                for (const auto& [bj, v] : grow[bi]) {
                    for (const auto& [lj, vj] : rows[bj]) out[lj] += v * vj;
                    Gx[bi] += v * xpv[bj];
                }
            });
            parallel_for(nl, [&](std::size_t i) {
                double* out = &r.G[w][i * nl];
                double gi = 0;
                for (const auto& [k, v] : M[i]) {
                    const double* row = &GM[k * nl];
                    for (std::size_t j = 0; j < nl; ++j) out[j] += v * row[j];
                    gi += v * Gx[k];
                }
                r.g[w][i] = gi;
            });
            for (std::size_t k = 0; k < nb; ++k) r.s[w] += xpv[k] * Gx[k];
        }
    // Plated holes (vertical filaments): wall resistance and internal impedance per length.
    for (double f : s_.frequencies_hz) {
        const double w = 2 * kPi * f;
        std::vector<cplx> Z(nl * nl), c(nl);
        for (std::size_t q = 0; q < nl * nl; ++q) Z[q] = cplx(0, w * ML[q]);
        for (std::size_t q = 0; q < nl; ++q) c[q] = cplx(0, -w * minus_mlxp[q]);
        cplx sxx(0, w * xlx);
        for (const auto& [layer, r] : red) {
            const cplx zs = sheet_impedance(f, r.sigma, r.t);
            // The share of the skin effect: 0 at DC, 1 when the impedance has left its DC value.
            const double rdc = 1.0 / (r.sigma * r.t);
            const double chi = std::min(1.0, std::abs(zs - rdc) / std::abs(zs));
            for (std::size_t q = 0; q < nl * nl; ++q) Z[q] += zs * ((1 - chi) * r.G[0][q] + chi * r.G[1][q]);
            for (std::size_t q = 0; q < nl; ++q) c[q] += zs * ((1 - chi) * r.g[0][q] + chi * r.g[1][q]);
            sxx += zs * ((1 - chi) * r.s[0] + chi * r.s[1]);
        }
        for (std::size_t k = 0; k < D.fil.size(); ++k) {
            const Filament& fl = D.fil[k];
            if (!fl.vertical) continue;
            const double len = std::abs(fl.p1[2] - fl.p0[2]);
            const cplx z = len / (2 * kPi * fl.radius) * wall_impedance(f, 5.8e7, s_.plating_nm * 1e-9);
            const auto& rw = rows[nr + k];
            for (const auto& [li, vi] : rw)
                for (const auto& [lj, vj] : rw) Z[li * nl + lj] += z * vi * vj;
            for (const auto& [li, vi] : rw) c[li] += z * vi * xpv[nr + k];
            sxx += z * xpv[nr + k] * xpv[nr + k];
        }
        // Stitches: their series resistance and capacitance.
        for (std::size_t k = 0; k < D.fil.size(); ++k) {
            const Filament& fl = D.fil[k];
            if (fl.stitch < 0) continue;
            const cplx z = fl.resistance + (fl.capacitance > 0 ? cplx(0, -1.0 / (w * fl.capacitance)) : cplx(0));
            if (z == cplx(0)) continue;
            const auto& rw = rows[nr + k];
            for (const auto& [li, vi] : rw)
                for (const auto& [lj, vj] : rw) Z[li * nl + lj] += z * vi * vj;
            for (const auto& [li, vi] : rw) c[li] += z * vi * xpv[nr + k];
            sxx += z * xpv[nr + k] * xpv[nr + k];
        }
        // Loop currents y: Z y = -c; the input impedance x^T Z x = s + y^T c.
        std::vector<cplx> y(nl);
        for (std::size_t q = 0; q < nl; ++q) y[q] = -c[q];
        if (nl > 0 && !complex_symmetric_solve(Z, nl, y)) throw std::runtime_error("loop " + loop_.name + ": singular loop impedance matrix");
        cplx zin = sxx;
        for (std::size_t q = 0; q < nl; ++q) zin += y[q] * c[q];
        LoopImpedance li{f, zin.real(), zin.imag() / w * 1e9, {}};
        li.stitch_current.assign(loop_.stitches.size(), cplx(0));
        for (std::size_t k = 0; k < D.fil.size(); ++k) {
            const Filament& fl = D.fil[k];
            if (fl.stitch < 0) continue;
            cplx I = xpv[nr + k];
            for (const auto& [lj, vj] : rows[nr + k]) I += vj * y[lj];
            li.stitch_current[fl.stitch] = I;
        }
        out.push_back(std::move(li));
    }
}

// Partial inductances of the filaments with the RWG branches (F, nf x nr, a row per filament)
// and with each other (FF, nf x nf), H.
void LoopSolver::filament_blocks(const Discretisation& D, std::vector<double>& F, std::vector<double>& FF) const {
    // Filaments with the sheets (horizontal filaments only: vertical ones are orthogonal).
    const std::size_t nr = D.nrwg, nf = D.fil.size(), nt = D.tris.size();
    F.assign(nf * nr, 0.0);
    FF.assign(nf * nf, 0.0);
    static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
    static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    for (std::size_t k = 0; k < nf; ++k) {
        const Filament& f = D.fil[k];
        if (f.vertical) continue;
        const double ux = f.p1[0] - f.p0[0], uy = f.p1[1] - f.p0[1];
        const double len = std::hypot(ux, uy);
        if (len <= 0) continue;
        const double tx = ux / len, ty = uy / len;
        std::vector<double> col(nr, 0.0);
        std::vector<std::vector<double>> part(nt);
        parallel_for(nt, [&](std::size_t t) {
            const Tri& T = D.tris[t];
            const double dist = std::max(1e-9, dist_point_segment(T.c.x, T.c.y, f.p0[0], f.p0[1], f.p1[0], f.p1[1]));
            const double dd = std::hypot(dist, f.p0[2] - T.z);
            const int pieces = std::clamp(static_cast<int>(std::ceil(2.0 * len / std::max(dd, T.size))), 1, 64);
            std::array<double, 3> acc{0, 0, 0};
            for (int p = 0; p < pieces; ++p)
                for (int g = 0; g < 8; ++g) {
                    const double s = (p + 0.5 * (gx[g] + 1)) / pieces;
                    const double w = gw[g] * 0.5 / pieces * len;
                    const Vec2 r{f.p0[0] + s * ux, f.p0[1] + s * uy};
                    const auto pot = triangle_potential(T.v, r, f.p0[2] - T.z);
                    for (int j = 0; j < 3; ++j) {
                        const double ax = pot.vector.x + (r.x - T.v[j].x) * pot.scalar;
                        const double ay = pot.vector.y + (r.y - T.v[j].y) * pot.scalar;
                        acc[j] += w * (tx * ax + ty * ay);
                    }
                }
            part[t].assign(acc.begin(), acc.end());
        });
        for (std::size_t t = 0; t < nt; ++t) {
            const Tri& T = D.tris[t];
            for (int j = 0; j < 3; ++j)
                if (T.basis[j] >= 0) col[static_cast<std::size_t>(T.basis[j])] += kMu0Over4Pi * T.sign[j] / (2.0 * T.area) * part[t][j];
        }
        std::copy(col.begin(), col.end(), F.begin() + static_cast<std::ptrdiff_t>(k * nr));
    }
    // Filaments with filaments.
    for (std::size_t a = 0; a < nf; ++a)
        for (std::size_t c = a; c < nf; ++c) {
            const Filament &f = D.fil[a], &g = D.fil[c];
            double m = 0;
            if (a == c) {
                const double len = std::sqrt(std::pow(f.p1[0] - f.p0[0], 2) + std::pow(f.p1[1] - f.p0[1], 2) + std::pow(f.p1[2] - f.p0[2], 2));
                m = parallel_filaments(0, len, 0, len, f.radius) + f.extra;
            } else if (f.vertical && g.vertical) {
                const double rho = std::max(std::hypot(f.p0[0] - g.p0[0], f.p0[1] - g.p0[1]), std::max(f.radius, g.radius));
                const double s1 = f.p1[2] > f.p0[2] ? 1.0 : -1.0, s2 = g.p1[2] > g.p0[2] ? 1.0 : -1.0;
                m = s1 * s2 * parallel_filaments(std::min(f.p0[2], f.p1[2]), std::max(f.p0[2], f.p1[2]), std::min(g.p0[2], g.p1[2]),
                                                 std::max(g.p0[2], g.p1[2]), rho);
            } else if (!f.vertical && !g.vertical) {
                m = segment_mutual(f.p0, f.p1, g.p0, g.p1);
            }
            FF[a * nf + c] = FF[c * nf + a] = m;
        }

}

void LoopSolver::loop_system_streamed(const Discretisation& D, const std::vector<Column>& M, const Column& xp, std::vector<double>& Z,
                                      std::vector<double>& rhs, double& xlx) const {
    const std::size_t nr = D.nrwg, nf = D.fil.size(), nb = nr + nf, nt = D.tris.size(), nl = M.size();
    std::vector<double> F, FF;
    filament_blocks(D, F, FF);
    std::vector<double> xpv(nb, 0.0);
    for (const auto& [k, v] : xp) xpv[k] += v;
    std::fill(Z.begin(), Z.end(), 0.0);
    std::fill(rhs.begin(), rhs.end(), 0.0);
    xlx = 0;
    // Rows of a block: a branch and its part of L's row (W, rows x nb).
    struct Row {
        std::size_t branch;
    };
    std::vector<Row> rows;
    std::vector<double> W, Y, yx;
    // Loops through each branch, for the rows of a block: rows_of[branch] lists its rows.
    std::vector<std::array<long, 2>> rows_of(nb, {-1, -1});
    auto multiply_out = [&]() {
        const std::size_t nrow = rows.size();
        Y.assign(nrow * nl, 0.0);
        yx.assign(nrow, 0.0);
        parallel_for(nrow, [&](std::size_t r) {
            const double* w = &W[r * nb];
            double* y = &Y[r * nl];
            for (std::size_t j = 0; j < nl; ++j) {
                double t = 0;
                for (const auto& [k, v] : M[j]) t += v * w[k];
                y[j] = t;
            }
            double t = 0;
            for (const auto& [k, v] : xp) t += v * w[k];
            yx[r] = t;
        });
        for (std::size_t r = 0; r < nrow; ++r) {
            auto& ro = rows_of[rows[r].branch];
            (ro[0] < 0 ? ro[0] : ro[1]) = static_cast<long>(r);
        }
        parallel_for(nl, [&](std::size_t p) {
            double* z = &Z[p * nl];
            double c = 0;
            for (const auto& [k, v] : M[p])
                for (const long r : rows_of[k]) {
                    if (r < 0) continue;
                    const double* y = &Y[static_cast<std::size_t>(r) * nl];
                    for (std::size_t q = 0; q < nl; ++q) z[q] += v * y[q];
                    c += v * yx[static_cast<std::size_t>(r)];
                }
            rhs[p] -= c;
        });
        for (std::size_t r = 0; r < nrow; ++r) {
            xlx += xpv[rows[r].branch] * yx[r];
            rows_of[rows[r].branch] = {-1, -1};
        }
    };
    // Triangle blocks: each triangle's part of the rows of its RWG branches, against all
    // triangles (for u < t the kernel of (u, t) transposed, as the dense assembly has it),
    // and the filaments' columns once per branch (on its positive triangle).
    const std::size_t per_block = std::max<std::size_t>(1, static_cast<std::size_t>(256e6 / (8.0 * 3.0 * double(nb))));
    for (std::size_t t0 = 0; t0 < nt; t0 += per_block) {
        const std::size_t t1 = std::min(nt, t0 + per_block);
        rows.clear();
        std::vector<std::array<long, 3>> row_of(t1 - t0, {-1, -1, -1});
        for (std::size_t t = t0; t < t1; ++t)
            for (int i = 0; i < 3; ++i)
                if (D.tris[t].basis[i] >= 0) {
                    row_of[t - t0][i] = static_cast<long>(rows.size());
                    rows.push_back({static_cast<std::size_t>(D.tris[t].basis[i])});
                }
        W.assign(rows.size() * nb, 0.0);
        parallel_for(t1 - t0, [&](std::size_t lt) {
            const std::size_t t = t0 + lt;
            const Tri& T = D.tris[t];
            for (std::size_t u = 0; u < nt; ++u) {
                const Tri& U = D.tris[u];
                std::array<std::array<double, 3>, 3> K;
                if (u >= t) {
                    K = pair_kernel(T, U);
                    if (u == t)
                        for (int i = 0; i < 3; ++i)
                            for (int j = 0; j < i; ++j) K[i][j] = K[j][i] = 0.5 * (K[i][j] + K[j][i]);
                } else {
                    const auto KT = pair_kernel(U, T);
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j) K[i][j] = KT[j][i];
                }
                for (int i = 0; i < 3; ++i) {
                    if (T.basis[i] < 0) continue;
                    const double fi = T.sign[i] / (2.0 * T.area);
                    double* w = &W[static_cast<std::size_t>(row_of[lt][i]) * nb];
                    for (int j = 0; j < 3; ++j) {
                        if (U.basis[j] < 0) continue;
                        w[static_cast<std::size_t>(U.basis[j])] += kMu0Over4Pi * fi * U.sign[j] / (2.0 * U.area) * K[i][j];
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                if (T.basis[i] < 0 || T.sign[i] < 0) continue;
                double* w = &W[static_cast<std::size_t>(row_of[lt][i]) * nb];
                const std::size_t a = static_cast<std::size_t>(T.basis[i]);
                for (std::size_t k = 0; k < nf; ++k) w[nr + k] += F[k * nr + a];
            }
        });
        multiply_out();
    }
    // The filaments' rows.
    rows.clear();
    W.assign(nf * nb, 0.0);
    for (std::size_t k = 0; k < nf; ++k) {
        rows.push_back({nr + k});
        std::copy(F.begin() + static_cast<std::ptrdiff_t>(k * nr), F.begin() + static_cast<std::ptrdiff_t>((k + 1) * nr), W.begin() + static_cast<std::ptrdiff_t>(k * nb));
        for (std::size_t c = 0; c < nf; ++c) W[k * nb + nr + c] = FF[k * nf + c];
    }
    if (nf > 0) multiply_out();
    for (std::size_t p = 0; p < nl; ++p)
        for (std::size_t q = 0; q < p; ++q) Z[p * nl + q] = Z[q * nl + p] = 0.5 * (Z[p * nl + q] + Z[q * nl + p]);
}

double LoopSolver::solve(const Discretisation& D, LoopResult* out, std::vector<std::array<double, 2>>* sheet_out, bool sweep) {
    if (out) out->impedance.clear();
    const std::size_t nr = D.nrwg, nf = D.fil.size(), nb = nr + nf;
    const std::size_t nt = D.tris.size();
    // Partial inductance matrix (H), dense, all branches; in the basis of loops, when it
    // would exceed dense_matrix_mb, it is formed by blocks of rows instead.
    std::vector<double> L;
    const bool streamed = s_.loop_basis && double(nb) * double(nb) * 8.0 > s_.dense_matrix_mb * 1048576.0;
    auto assemble = [&]() {
    L.assign(nb * nb, 0.0);
    parallel_for(nt, [&](std::size_t t) {
        const Tri& T = D.tris[t];
        for (std::size_t u = t; u < nt; ++u) {
            const Tri& U = D.tris[u];
            auto K = pair_kernel(T, U);
            if (u == t)
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < i; ++j) K[i][j] = K[j][i] = 0.5 * (K[i][j] + K[j][i]);
            for (int i = 0; i < 3; ++i) {
                if (T.basis[i] < 0) continue;
                const double fi = T.sign[i] / (2.0 * T.area);
                for (int j = 0; j < 3; ++j) {
                    if (U.basis[j] < 0) continue;
                    const double v = kMu0Over4Pi * fi * U.sign[j] / (2.0 * U.area) * K[i][j];
                    const auto a = static_cast<std::size_t>(T.basis[i]), c = static_cast<std::size_t>(U.basis[j]);
                    atomic_add(L[a * nb + c], v);
                    if (u != t) atomic_add(L[c * nb + a], v);
                }
            }
        }
    });
    // Filaments with the sheets and with each other.
    std::vector<double> F, FF;
    filament_blocks(D, F, FF);
    for (std::size_t k = 0; k < nf; ++k)
        for (std::size_t r = 0; r < nr; ++r) L[r * nb + nr + k] = L[(nr + k) * nb + r] = F[k * nr + r];
    for (std::size_t a = 0; a < nf; ++a)
        for (std::size_t c = 0; c < nf; ++c) L[(nr + a) * nb + nr + c] = FF[a * nf + c];
    };
    if (!streamed) assemble();

    // Kirchhoff's current law: an RWG current leaves one triangle and enters another; a
    // filament takes its current from its first set of triangles and returns it to the second.
    std::vector<std::array<long, 2>> ends(nr, {-1, -1});
    for (const auto& T : D.tris)
        for (int k = 0; k < 3; ++k)
            if (T.basis[k] >= 0) ends[static_cast<std::size_t>(T.basis[k])][T.sign[k] > 0 ? 0 : 1] = T.node;
    auto column = [&](std::size_t k, auto&& put) {
        if (k < nr) {
            put(ends[k][0], 1.0);
            put(ends[k][1], -1.0);
            return;
        }
        const Filament& f = D.fil[k - nr];
        for (const auto& [t, w] : f.end0) put(static_cast<long>(t), w);
        for (const auto& [t, w] : f.end1) put(static_cast<long>(t), -w);
    };
    const std::size_t src = nr;  // the first element carries the loop current
    // Components of the node graph without the source; one equation dropped in each.
    UnionFind uf(D.nnodes);
    for (std::size_t k = 0; k < nb; ++k) {
        if (k == src) continue;
        long first = -1;
        column(k, [&](long node, double) {
            if (node < 0) return;
            if (first < 0) first = node;
            else uf.unite(static_cast<std::size_t>(node), static_cast<std::size_t>(first));
        });
    }
    {
        long root = -1;
        bool joined = true;
        column(src, [&](long node, double) {
            const auto r = static_cast<long>(uf.find(static_cast<std::size_t>(node)));
            if (root < 0) root = r;
            joined = joined && r == root;
        });
        if (!joined) throw std::runtime_error("loop " + loop_.name + " is open: no copper path joins its elements");
    }
    // Free branches: all but the source, re-indexed.
    const std::size_t nF = nb - 1;
    auto fidx = [&](std::size_t k) { return k < src ? k : k - 1; };
    std::vector<double> x(nF);
    double ind = 0;
    const bool by_loops = s_.loop_basis && solve_by_loops(D, streamed ? nullptr : &L, ends, src,
                                                          sweep && out && !s_.frequencies_hz.empty() ? &out->impedance : nullptr, x, ind);
    if (!by_loops && L.empty()) assemble();
    if (!by_loops) {
        std::vector<long> row(D.nnodes, -1);
        std::set<std::size_t> dropped;
        std::size_t nrows = 0;
        for (std::size_t v = 0; v < D.nnodes; ++v) {
            if (dropped.insert(uf.find(v)).second) continue;  // the first node of each component
            row[v] = static_cast<long>(nrows++);
        }
        std::vector<double> LF(nF * nF), l(nF), ds(nrows, 0.0);
        for (std::size_t a = 0; a < nb; ++a) {
            if (a == src) continue;
            for (std::size_t c = 0; c < nb; ++c)
                if (c != src) LF[fidx(a) * nF + fidx(c)] = L[a * nb + c];
            l[fidx(a)] = L[a * nb + src];
        }
        const double lss = L[src * nb + src];
        column(src, [&](long node, double w) {
            const long r = row[static_cast<std::size_t>(node)];
            if (r >= 0) ds[static_cast<std::size_t>(r)] += w;
        });
        std::vector<double>().swap(L);
        if (!cholesky(LF, nF)) throw std::runtime_error("loop " + loop_.name + ": the inductance matrix is not positive definite");
        // W = C^-1 D_F^T (nF x nrows), w0 = C^-1 l.
        std::vector<double> W(nF * nrows, 0.0);
        for (std::size_t k = 0; k < nb; ++k) {
            if (k == src) continue;
            column(k, [&](long node, double w) {
                if (node < 0) return;
                const long r = row[static_cast<std::size_t>(node)];
                if (r >= 0) W[fidx(k) * nrows + static_cast<std::size_t>(r)] += w;
            });
        }
        forward_solve(LF, nF, W, nrows);
        std::vector<double> w0(l);
        forward_solve(LF, nF, w0, 1);
        // S = W^T W, rhs = d_s - W^T w0.
        std::vector<double> S(nrows * nrows, 0.0), rhs(nrows);
        parallel_for(nrows, [&](std::size_t i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double t = 0;
                for (std::size_t k = 0; k < nF; ++k) t += W[k * nrows + i] * W[k * nrows + j];
                S[i * nrows + j] = S[j * nrows + i] = t;
            }
            double t = 0;
            for (std::size_t k = 0; k < nF; ++k) t += W[k * nrows + i] * w0[k];
            rhs[i] = ds[i] - t;
        });
        if (!cholesky(S, nrows)) throw std::runtime_error("loop " + loop_.name + ": the current law has dependent equations");
        std::vector<double> lam(rhs);
        forward_solve(S, nrows, lam, 1);
        backward_solve(S, nrows, lam, 1);
        // x = -C^-T (w0 + W lambda)
        parallel_for(nF, [&](std::size_t k) {
            double t = w0[k];
            for (std::size_t i = 0; i < nrows; ++i) t += W[k * nrows + i] * lam[i];
            x[k] = -t;
        });
        backward_solve(LF, nF, x, 1);
        ind = lss;
        for (std::size_t k = 0; k < nF; ++k) ind += l[k] * x[k];
        for (std::size_t i = 0; i < nrows; ++i) ind += ds[i] * lam[i];

    }
    {
        // Kirchhoff's current law at every node, with the source.
        std::vector<double> kcl(D.nnodes, 0.0);
        for (std::size_t k = 0; k < nb; ++k) {
            const double I = k == src ? 1.0 : x[fidx(k)];
            column(k, [&](long node, double w) {
                if (node >= 0) kcl[static_cast<std::size_t>(node)] += w * I;
            });
        }
        double worst = 0;
        for (double v : kcl) worst = std::max(worst, std::abs(v));
        if (out && worst > 1e-6) out->warnings.push_back("the current law holds to " + std::to_string(worst) + " A only");
    }
    auto current = [&](std::size_t k) { return k == src ? 1.0 : x[fidx(k)]; };
    if (sheet_out) {
        sheet_out->assign(nt, {0.0, 0.0});
        for (std::size_t t = 0; t < nt; ++t) {
            const Tri& T = D.tris[t];
            for (int j = 0; j < 3; ++j) {
                if (T.basis[j] < 0) continue;
                const double I = current(static_cast<std::size_t>(T.basis[j])) * T.sign[j] / (2.0 * T.area);
                (*sheet_out)[t][0] += I * (T.c.x - T.v[j].x);
                (*sheet_out)[t][1] += I * (T.c.y - T.v[j].y);
            }
        }
    }
    if (out) {
        out->triangles = nt;
        out->unknowns = nb;
        out->filaments = nf;
        out->stitch_current.assign(loop_.stitches.size(), 0.0);
        for (std::size_t k = 0; k < nf; ++k)
            if (D.fil[k].stitch >= 0) out->stitch_current[D.fil[k].stitch] = current(nr + k);
        // Sheet currents at the centroids.
        std::vector<std::array<double, 2>> sheet(nt, {0.0, 0.0});
        for (std::size_t t = 0; t < nt; ++t) {
            const Tri& T = D.tris[t];
            for (int j = 0; j < 3; ++j) {
                if (T.basis[j] < 0) continue;
                const double I = current(static_cast<std::size_t>(T.basis[j])) * T.sign[j] / (2.0 * T.area);
                sheet[t][0] += I * (T.c.x - T.v[j].x);
                sheet[t][1] += I * (T.c.y - T.v[j].y);
            }
        }
        // Current elements for the radiated field and the moment: each triangle's sheet current times its
        // area at its centroid, each conductor in pieces of at most 1 mm.
        out->currents.clear();
        for (std::size_t t = 0; t < nt; ++t) {
            const Tri& T = D.tris[t];
            if (sheet[t][0] == 0 && sheet[t][1] == 0) continue;
            out->currents.push_back({{T.c.x, T.c.y, T.z}, {sheet[t][0] * T.area, sheet[t][1] * T.area, 0.0}});
        }
        for (std::size_t k = 0; k < nf; ++k) {
            const Filament& f = D.fil[k];
            const double I = current(nr + k);
            const std::array<double, 3> dl{f.p1[0] - f.p0[0], f.p1[1] - f.p0[1], f.p1[2] - f.p0[2]};
            const double len = std::sqrt(dl[0] * dl[0] + dl[1] * dl[1] + dl[2] * dl[2]);
            const int pieces = std::max(1, static_cast<int>(std::ceil(len / 1e-3)));
            for (int q = 0; q < pieces; ++q) {
                const double u = (q + 0.5) / pieces;
                out->currents.push_back({{f.p0[0] + u * dl[0], f.p0[1] + u * dl[1], f.p0[2] + u * dl[2]},
                                         {I * dl[0] / pieces, I * dl[1] / pieces, I * dl[2] / pieces}});
            }
            // The current enters the conductor from the centroid of its share of the sheet
            // and leaves it to the other end's: short pieces join them, so that the current
            // is closed exactly. An offset of micrometres would otherwise radiate as an
            // electric dipole.
            auto join = [&](const std::vector<std::pair<std::size_t, double>>& end, const std::array<double, 3>& p, double sign) {
                std::array<double, 3> c{0, 0, 0};
                double w = 0;
                for (const auto& [t, share] : end) {
                    c[0] += share * D.tris[t].c.x;
                    c[1] += share * D.tris[t].c.y;
                    c[2] += share * D.tris[t].z;
                    w += share;
                }
                if (!(w > 0)) return;
                for (auto& v : c) v /= w;
                out->currents.push_back({{0.5 * (c[0] + p[0]), 0.5 * (c[1] + p[1]), 0.5 * (c[2] + p[2])},
                                         {sign * I * (p[0] - c[0]), sign * I * (p[1] - c[1]), sign * I * (p[2] - c[2])}});
            };
            join(f.end0, f.p0, 1.0);
            join(f.end1, f.p1, -1.0);
        }
        // Magnetic moment m = 1/2 sum of r x I dl of the closed current: with the joining
        // pieces it does not depend on the origin. A triangle's linear current adds nothing
        // to its centroid's term.
        out->moment_m2 = {0, 0, 0};
        for (const auto& e : out->currents) {
            const auto& x = e.position_m;
            const auto& d = e.moment_a_m;
            out->moment_m2[0] += 0.5 * (x[1] * d[2] - x[2] * d[1]);
            out->moment_m2[1] += 0.5 * (x[2] * d[0] - x[0] * d[2]);
            out->moment_m2[2] += 0.5 * (x[0] * d[1] - x[1] * d[0]);
        }
        // Fields per layer for the viewer.
        std::size_t t0 = 0;
        std::map<int, std::size_t> slot;
        for (std::size_t j = 0; j < D.meshes.size(); ++j) {
            const TriMesh& tm = D.meshes[j];
            const int layer = D.mesh_layer[j];
            if (!slot.count(layer)) {
                slot[layer] = out->sheets.size();
                out->sheets.push_back({});
                out->sheets.back().layer = layer;
            }
            auto& f = out->sheets[slot[layer]];
            const auto base = static_cast<std::uint32_t>(f.x.size());
            f.x.insert(f.x.end(), tm.x.begin(), tm.x.end());
            f.y.insert(f.y.end(), tm.y.begin(), tm.y.end());
            for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
                const auto& tr = tm.triangles[t];
                f.triangles.push_back({base + tr[0], base + tr[1], base + tr[2]});
                f.current_a_per_m.push_back(sheet[t0 + t]);
            }
            t0 += tm.triangles.size();
        }
        std::sort(out->sheets.begin(), out->sheets.end(), [](const auto& a, const auto& c) { return a.layer < c.layer; });
    }
    return ind;
}

LoopResult LoopSolver::run(bool sweep) {
    const auto t0 = Clock::now();
    LoopResult res;
    res.name = loop_.name;
    setup();
    for (const auto& e : loop_.elements) res.element_inductance_nh += e.extra_inductance_nh;
    // The first mesh from the a-priori size field, with about a quarter of the triangle
    // budget; then adaptive passes, each refining where the recovered sheet current
    // differs most from the computed one, until the budget is used.
    double scale = 1.0;
    Seeds seeds;
    Discretisation D;
    const std::size_t first = std::max<std::size_t>(200, s_.max_triangles / 4);
    for (int attempt = 0; attempt < 8; ++attempt) {
        seeds = initial_seeds(scale);
        D = discretise(seeds, scale);
        const std::size_t nt = D.tris.size();
        if (nt <= first * 1.2 && nt >= first * 0.6) break;
        // The outlines set a floor to the count; the elements at the loop are not made more
        // than four times coarser than the a-priori size to meet the budget.
        const double next = scale * std::clamp(std::sqrt(double(nt) / first), 0.3, 3.0);
        if (next > 4.0 && scale >= 4.0) break;
        scale = std::min(next, 4.0);
    }
    // At the limit of four the budget, not the loop, sets the elements: the inductance is
    // then several per cent high (docs/solvers/emi.md, "Accuracy and cost").
    char coarse[240];
    std::snprintf(coarse, sizeof coarse,
                  "the triangle budget allows only the coarsest elements at the loop (%.2f mm): the inductance is likely several per cent "
                  "high; raise loop_settings.max_triangles",
                  h0_ * scale * 1e-6);
    const bool at_coarsest = scale >= 4.0;
    if (D.tris.size() > s_.max_triangles)
        res.warnings.push_back("the copper outlines alone need " + std::to_string(D.tris.size()) + " triangles, above the budget of " +
                               std::to_string(s_.max_triangles));
    if (D.tris.size() > 3 * s_.max_triangles)
        throw std::runtime_error("loop " + loop_.name + ": the copper around the loop needs " + std::to_string(D.tris.size()) +
                                 " triangles at the coarsest size; reduce window_margin_mm or raise max_triangles");
    std::vector<double> history;
    for (int pass = 0;; ++pass) {
        LoopResult r;
        std::vector<std::array<double, 2>> sheet;
        // The frequency sweep on the last pass only: the budget is reached or the passes end.
        const bool last = sweep && (pass >= 7 || D.tris.size() >= 0.7 * s_.max_triangles);
        const double ind = solve(D, &r, &sheet, last);
        history.push_back(ind);
        r.name = res.name;
        r.element_inductance_nh = res.element_inductance_nh;
        r.warnings.insert(r.warnings.begin(), res.warnings.begin(), res.warnings.end());
        res = std::move(r);
        const std::size_t nt = D.tris.size();
        if (pass >= 7 || nt >= 0.7 * s_.max_triangles) break;
        // Recovery estimate: sheet current averaged at the vertices against each triangle's.
        std::vector<std::vector<double>> eta2(jobs_.size());
        std::vector<std::vector<std::array<double, 3>>> vsum(jobs_.size());  // per vertex: area, Kx A, Ky A
        for (std::size_t j = 0; j < jobs_.size(); ++j) {
            eta2[j].assign(D.meshes[j].triangles.size(), 0.0);
            vsum[j].assign(D.meshes[j].x.size(), {0, 0, 0});
        }
        for (std::size_t t = 0; t < nt; ++t) {
            const Tri& T = D.tris[t];
            for (auto v : D.meshes[T.job].triangles[T.local]) {
                auto& a = vsum[T.job][v];
                a[0] += T.area;
                a[1] += T.area * sheet[t][0];
                a[2] += T.area * sheet[t][1];
            }
        }
        double total = 0, energy = 0;
        for (std::size_t t = 0; t < nt; ++t) {
            const Tri& T = D.tris[t];
            double e = 0;
            for (auto v : D.meshes[T.job].triangles[T.local]) {
                const auto& a = vsum[T.job][v];
                const double dx = a[1] / a[0] - sheet[t][0], dy = a[2] / a[0] - sheet[t][1];
                e += (dx * dx + dy * dy) / 3.0;
            }
            eta2[T.job][T.local] = e * T.area;
            total += e * T.area;
            energy += (sheet[t][0] * sheet[t][0] + sheet[t][1] * sheet[t][1]) * T.area;
        }
        std::vector<RefinePart> parts;
        for (std::size_t j = 0; j < jobs_.size(); ++j) parts.push_back({&D.meshes[j], &eta2[j], j});
        bool budget = false;
        const std::size_t cap = std::min(s_.max_triangles, static_cast<std::size_t>(2.2 * nt));
        const std::size_t marked = mark_and_seed(parts, total, energy, 1e-3, nt, cap, 0.02 * h0_ * scale, seeds, budget);
        if (marked == 0) break;
        D = discretise(seeds, scale);
    }
    last_ = std::move(D);
    res.inductance_nh = history.back() * 1e9;
    res.element_size_mm = h0_ * scale * 1e-6;
    if (at_coarsest) res.warnings.push_back(coarse);
    if (s_.estimate_error && history.size() >= 2) res.inductance_coarse_nh = history[history.size() - 2] * 1e9;
    res.window_area_mm2 = area_nm2(window_) * 1e-12;
    res.window_margin_mm = s_.window_margin_nm * 1e-6;
    res.seconds = since(t0);
    return res;
}

// Radiated field of the harmonics: the field per ampere on a logarithmic grid of 24 points
// per decade, interpolated to each harmonic.
void radiate(LoopResult& res, const LoopSpec& loop, const LoopSettings& s) {
    if (!(loop.waveform.amplitude_a > 0 && loop.waveform.frequency_hz > 0)) return;
    const auto harmonics = trapezoid_harmonics(loop.waveform, s.max_frequency_hz);
    if (harmonics.empty()) return;
    res.field_site = s.field_site;
    res.field_distance_m = s.distance_m;
    const double f0 = harmonics.front().first, f1 = harmonics.back().first;
    const int n = std::max(2, static_cast<int>(std::ceil(24 * std::log10(f1 / f0))) + 1);
    std::vector<double> lf(n), le(n);
    for (int i = 0; i < n; ++i) {
        const double f = f0 * std::pow(f1 / f0, n > 1 ? double(i) / (n - 1) : 0.0);
        const double e = loop_field_per_ampere(res.currents, res.moment_m2, f, s);
        res.field_per_ampere.push_back({f, e});
        lf[i] = std::log(f);
        le[i] = std::log(std::max(e, 1e-300));
    }
    for (const auto& [f, a] : harmonics) {
        const double x = std::log(f);
        std::size_t i = std::upper_bound(lf.begin(), lf.end(), x) - lf.begin();
        i = std::clamp<std::size_t>(i, 1, n - 1);
        const double u = std::clamp((x - lf[i - 1]) / (lf[i] - lf[i - 1]), 0.0, 1.0);
        LoopHarmonic h{f, a, a * std::exp(le[i - 1] + u * (le[i] - le[i - 1]))};
        if (h.field_v_per_m > res.max_field_v_per_m) {
            res.max_field_v_per_m = h.field_v_per_m;
            res.max_field_frequency_hz = f;
        }
        res.harmonics.push_back(h);
    }
}

std::vector<LoopImpedance> LoopSolver::sweep() {
    LoopResult r;
    std::vector<std::array<double, 2>> sheet;
    solve(last_, &r, &sheet, true);
    return r.impedance;
}

}  // namespace

LoopGeometry detail::loop_geometry(const Board& board, const LoopSpec& loop, const LoopSettings& settings) {
    return LoopSolver(board, loop, settings).geometry();
}

LoopResult solve_loop(const Board& board, const LoopSpec& loop, const LoopSettings& settings) {
    const auto t0 = Clock::now();
    if (!(settings.moment_tolerance > 0)) {
        LoopResult r = LoopSolver(board, loop, settings).run();
        radiate(r, loop, settings);
        return r;
    }
    // The moment over a plane depends on how much of the plane the window holds: grow it,
    // without the frequency sweep, which then runs once on the last mesh of the window that
    // settled.
    auto norm = [](const std::array<double, 3>& m) { return std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]); };
    LoopSettings s = settings;
    auto solver = std::make_unique<LoopSolver>(board, loop, s);
    LoopResult r = solver->run(false);
    const double board_area = area_nm2(board.outline) * 1e-12;
    auto whole_board = [&](const LoopResult& x) { return board_area > 0 && x.window_area_mm2 >= board_area * (1 - 1e-6); };
    bool converged = whole_board(r);
    for (int i = 0; !converged && i < settings.max_window_doublings; ++i) {
        s.window_margin_nm *= 2;
        solver = std::make_unique<LoopSolver>(board, loop, s);
        LoopResult g = solver->run(false);
        const bool same_window = std::abs(g.window_area_mm2 - r.window_area_mm2) <= 1e-6 * r.window_area_mm2;
        const double d = std::sqrt(std::pow(g.moment_m2[0] - r.moment_m2[0], 2) + std::pow(g.moment_m2[1] - r.moment_m2[1], 2) +
                                   std::pow(g.moment_m2[2] - r.moment_m2[2], 2));
        g.moment_change = norm(g.moment_m2) > 0 ? d / norm(g.moment_m2) : 0.0;
        r = std::move(g);
        if (same_window || whole_board(r) || r.moment_change < settings.moment_tolerance) {
            converged = true;
            break;
        }
    }
    if (!settings.frequencies_hz.empty()) r.impedance = solver->sweep();
    radiate(r, loop, settings);
    if (!converged) {
        char text[200];
        std::snprintf(text, sizeof text, "the moment changed by %.1f %% when the window's margin last doubled (to %.1f mm); it is uncertain by about as much",
                      100 * r.moment_change, r.window_margin_mm);
        r.warnings.push_back(text);
    }
    r.seconds = since(t0);
    return r;
}

std::vector<std::pair<double, double>> trapezoid_harmonics(const LoopWaveform& w, double f_max_hz) {
    // One period from t = 0: rise, high, fall, low. The duty cycle is the width at half
    // amplitude. Coefficients of a continuous piecewise-linear periodic function from the
    // jumps of its slope: c_n = -1/(T (n w0)^2) sum dk exp(-j n w0 tk).
    const double T = 1.0 / w.frequency_hz, tr = w.rise_s, tf = w.fall_s < 0 ? w.rise_s : w.fall_s;
    const double width = w.duty * T;
    if (tr < 0 || tf < 0 || width - 0.5 * (tr + tf) < -1e-15 || width + 0.5 * (tr + tf) > T + 1e-15)
        throw std::runtime_error("waveform: rise and fall times do not fit the pulse and the period");
    const double a = w.amplitude_a;
    // Corners: start of rise, top, start of fall, bottom; slope jumps.
    const double t1 = 0, t2 = tr, t3 = width + 0.5 * (tr - tf), t4 = t3 + tf;
    std::vector<std::pair<double, double>> corner;  // time, jump of slope
    const double up = tr > 0 ? a / tr : 0, down = tf > 0 ? -a / tf : 0;
    corner.push_back({t1, up});
    corner.push_back({t2, -up});
    corner.push_back({t3, down});
    corner.push_back({t4, -down});
    std::vector<std::pair<double, double>> out;
    const double w0 = 2 * kPi / T;
    for (std::size_t n = 1; n * w.frequency_hz <= f_max_hz * (1 + 1e-12); ++n) {
        double re = 0, im = 0;
        if (tr > 0 && tf > 0) {
            for (const auto& [t, d] : corner) {
                re += d * std::cos(n * w0 * t);
                im -= d * std::sin(n * w0 * t);
            }
            re *= -1.0 / (T * n * n * w0 * w0);
            im *= -1.0 / (T * n * n * w0 * w0);
        } else {
            // Ideal edges: a rectangular pulse of the given width.
            const double tau = width, phase = n * w0 * tau;
            re = a / (T * n * w0) * std::sin(phase);
            im = a / (T * n * w0) * (std::cos(phase) - 1.0);
        }
        out.push_back({n * w.frequency_hz, 2.0 * std::hypot(re, im)});
    }
    return out;
}

double small_loop_field(double moment_a_m2, double frequency_hz, double distance_m) {
    // E = eta0 k^2 m / (4 pi r) in the plane of the loop.
    constexpr double c0 = 299792458.0, eta0 = 376.730313668;
    const double k = 2 * kPi * frequency_hz / c0;
    return eta0 * k * k * moment_a_m2 / (4 * kPi * distance_m);
}

double loop_field_per_ampere(const std::vector<CurrentElement>& currents, const std::array<double, 3>& moment_m2, double frequency_hz,
                             const LoopSettings& s, bool far_only) {
    constexpr double c0 = 299792458.0, mu0 = 4e-7 * kPi;
    const double r = s.distance_m;
    if (s.field_site == "small_loop")
        return small_loop_field(std::sqrt(moment_m2[0] * moment_m2[0] + moment_m2[1] * moment_m2[1] + moment_m2[2] * moment_m2[2]),
                                frequency_hz, r);
    if (s.field_site != "free" && s.field_site != "ground") throw std::runtime_error("field_site must be free, ground or small_loop");
    if (currents.empty()) return 0.0;
    const double w = 2 * kPi * frequency_hz, k = w / c0;
    // The currents are divergence-free (the loop is closed), so the scalar potential vanishes
    // and E = -j w A with A = mu0 / (4 pi) sum I dl exp(-j k R) / R, near field included.
    const cplx scale = cplx(0, -w * mu0 / (4 * kPi));
    // Centre: the centroid of the elements weighted by their moments.
    std::array<double, 3> c{0, 0, 0};
    double wsum = 0;
    for (const auto& e : currents) {
        const double m = std::sqrt(e.moment_a_m[0] * e.moment_a_m[0] + e.moment_a_m[1] * e.moment_a_m[1] + e.moment_a_m[2] * e.moment_a_m[2]);
        for (int i = 0; i < 3; ++i) c[i] += m * e.position_m[i];
        wsum += m;
    }
    if (wsum > 0)
        for (auto& v : c) v /= wsum;
    // A closed current has no symmetric first moment: int (r_i J_j + r_j J_i) dV
    // = -int r_i r_j div J dV = 0. Pads and triangles taken as points leave one of the order
    // of their size squared, which matters where a plane has cancelled the moment; its field
    // is removed: A_sym = mu0 / (4 pi) (j k + 1/R) exp(-j k R) / R Q R^, Q the symmetric part
    // of sum (x - c) I dl^T about the centre.
    std::array<std::array<double, 3>, 3> Q{};
    for (const auto& e : currents)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) Q[i][j] += 0.5 * ((e.position_m[i] - c[i]) * e.moment_a_m[j] + (e.position_m[j] - c[j]) * e.moment_a_m[i]);
    using V3 = std::array<cplx, 3>;
    auto remove_q = [&](V3& a, const std::array<double, 3>& p, const std::array<double, 3>& centre, bool image) {
        const double dx = p[0] - centre[0], dy = p[1] - centre[1], dz = p[2] - centre[2];
        const double R = std::sqrt(dx * dx + dy * dy + dz * dz);
        const std::array<double, 3> u{dx / R, dy / R, dz / R};
        // The image's Q: positions mirrored in z, moments with x and y reversed.
        const std::array<double, 3> sm{1, 1, -1};
        const cplx g = std::polar(1.0 / R, -k * R) * cplx(1.0 / R, k);
        for (int i = 0; i < 3; ++i) {
            double qu = 0;
            for (int j = 0; j < 3; ++j) qu += (image ? -sm[i] * sm[j] : 1.0) * Q[i][j] * u[j];
            a[i] -= qu * g;
        }
    };
    auto field_at = [&](const std::array<double, 3>& p, bool image) {
        V3 a{0.0, 0.0, 0.0};
        for (const auto& e : currents) {
            std::array<double, 3> x{e.position_m[0] - c[0], e.position_m[1] - c[1], e.position_m[2] - c[2]};
            std::array<double, 3> m = e.moment_a_m;
            if (s.field_site == "ground") x[2] += s.board_height_m;
            for (int pass = 0; pass < (image ? 2 : 1); ++pass) {
                if (pass == 1) {  // the image in the ground: horizontal currents reversed
                    x[2] = -x[2];
                    m[0] = -m[0];
                    m[1] = -m[1];
                }
                const double dx = p[0] - x[0], dy = p[1] - x[1], dz = p[2] - x[2];
                const double R = std::sqrt(dx * dx + dy * dy + dz * dz);
                const cplx g = std::polar(1.0 / R, -k * R);
                for (int i = 0; i < 3; ++i) a[i] += m[i] * g;
            }
        }
        const double h = s.field_site == "ground" ? s.board_height_m : 0.0;
        remove_q(a, p, {0, 0, h}, false);
        if (image) remove_q(a, p, {0, 0, -h}, true);
        for (auto& v : a) v *= scale;
        return a;
    };
    std::vector<std::array<double, 3>> points;
    std::vector<int> kind;  // per point: 0 total, 1 horizontal and vertical polarisation
    if (s.field_site == "free") {
        for (int it = 0; it <= 36; ++it)
            for (int ip = 0; ip < 36; ++ip) {
                if ((it == 0 || it == 36) && ip > 0) continue;
                const double th = it * kPi / 36, ph = ip * kPi / 18;
                points.push_back({std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th)});
            }
    } else {
        const double h0 = s.antenna_min_height_m, h1 = std::max(s.antenna_max_height_m, h0);
        const int nh = std::max(1, static_cast<int>(std::ceil((h1 - h0) / 0.25)) + 1);
        for (int ih = 0; ih < nh; ++ih)
            for (int ip = 0; ip < 36; ++ip) {
                const double ph = ip * kPi / 18, h = nh > 1 ? h0 + (h1 - h0) * ih / (nh - 1) : h0;
                points.push_back({r * std::cos(ph), r * std::sin(ph), h});
            }
    }
    std::vector<double> best(points.size(), 0.0);
    parallel_for(points.size(), [&](std::size_t q) {
        const auto& u = points[q];
        double e2 = 0;
        if (s.field_site == "ground") {
            const V3 E = field_at(u, true);
            const double ph = std::atan2(u[1], u[0]);
            const cplx eh = -std::sin(ph) * E[0] + std::cos(ph) * E[1];
            e2 = std::max(std::norm(eh), std::norm(E[2]));
        } else if (far_only) {
            // Radiation vector N = sum I dl exp(j k u . x); E = -j w mu0 / (4 pi r) N_perp.
            V3 n{0.0, 0.0, 0.0};
            for (const auto& e : currents) {
                const double ux = u[0] * (e.position_m[0] - c[0]) + u[1] * (e.position_m[1] - c[1]) + u[2] * (e.position_m[2] - c[2]);
                const cplx g = std::polar(1.0, k * ux);
                for (int i = 0; i < 3; ++i) n[i] += e.moment_a_m[i] * g;
            }
            for (int i = 0; i < 3; ++i) n[i] -= cplx(0, k) * (Q[i][0] * u[0] + Q[i][1] * u[1] + Q[i][2] * u[2]);
            const cplx un = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
            for (int i = 0; i < 3; ++i) e2 += std::norm((n[i] - un * u[i]) * scale / r);
        } else {
            const V3 E = field_at({r * u[0], r * u[1], r * u[2]}, false);
            e2 = std::norm(E[0]) + std::norm(E[1]) + std::norm(E[2]);
        }
        best[q] = e2;
    });
    return std::sqrt(*std::max_element(best.begin(), best.end()));
}

}  // namespace kairoq
