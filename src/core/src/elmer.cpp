// Export of a DC problem to Elmer (docs/solvers/elmer.md).
//
// One triangulation of the footprint of the ports' nets, with the copper outlines of every
// layer, the ports and the plated holes as constraints, is stacked at the faces of the
// copper layers: a triangle inside a layer's copper gives a wedge through that layer, a
// triangle inside a plated hole a column of wedges over the hole's span, of the plating's
// conductance spread over the hole. All layers share the nodes, so the volumes join
// conformally.
#include "kairoq/elmer.hpp"

#include "fem_util.hpp"
#include "segment_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace fs = std::filesystem;
using namespace detail;

namespace {

// Point-in-region tests for a polygon set.
struct Region {
    std::optional<SegmentGrid> grid;
    explicit Region(const PolygonSet& ps) {
        if (ps.empty()) return;
        std::vector<SegmentGrid::Seg> segs;
        auto add = [&](const Path64& p) {
            for (std::size_t i = 0; i < p.size(); ++i) {
                const auto& a = p[i];
                const auto& b = p[(i + 1) % p.size()];
                segs.push_back({double(a.x), double(a.y), double(b.x), double(b.y), true, 0.0});
            }
        };
        for (const auto& poly : ps) {
            add(poly.outer);
            for (const auto& h : poly.holes) add(h);
        }
        const Rect64 r = bounds(ps);
        grid.emplace(std::move(segs), double(r.left) - 1, double(r.top) - 1, double(r.right) + 1, double(r.bottom) + 1, 2.0e5);
    }
    bool inside(double x, double y) const { return grid && grid->inside(x, y); }
    bool near_boundary(double x, double y, double r) const { return grid && grid->any_within(x, y, r); }
};

// Makes nearly coincident outlines identical: vertices within eps of each other become one,
// and every edge takes the vertices within eps of it. Outlines of different layers that
// meet within rounding (a track's edge where it joins a pad on one layer, the pad's arc on
// the next) would otherwise give the shared triangulation slivers.
void clean_paths(const std::vector<Path64*>& paths, double eps) {
    const double cell = eps;
    auto key = [&](double x, double y) { return std::make_pair(static_cast<long long>(std::floor(x / cell)), static_cast<long long>(std::floor(y / cell))); };
    std::map<std::pair<long long, long long>, std::vector<Point64>> reps;
    auto rep_of = [&](Point64 p) {
        const auto [kx, ky] = key(double(p.x), double(p.y));
        for (long long dx = -1; dx <= 1; ++dx)
            for (long long dy = -1; dy <= 1; ++dy) {
                const auto it = reps.find({kx + dx, ky + dy});
                if (it == reps.end()) continue;
                for (const auto& r : it->second)
                    if (std::hypot(double(r.x - p.x), double(r.y - p.y)) <= eps) return r;
            }
        reps[{kx, ky}].push_back(p);
        return p;
    };
    for (auto* path : paths) {
        Path64 out;
        for (const auto& v : *path) {
            const Point64 r = rep_of(v);
            if (out.empty() || out.back() != r) out.push_back(r);
        }
        while (out.size() > 1 && out.front() == out.back()) out.pop_back();
        *path = std::move(out);
    }
    // Vertices on a coarser grid for the edge queries.
    const double vcell = std::max(eps, 50000.0);
    std::map<std::pair<long long, long long>, std::vector<Point64>> verts;
    for (const auto& [k, list] : reps)
        for (const auto& r : list) verts[{static_cast<long long>(std::floor(r.x / vcell)), static_cast<long long>(std::floor(r.y / vcell))}].push_back(r);
    for (auto* path : paths) {
        if (path->size() < 3) continue;
        Path64 out;
        for (std::size_t i = 0; i < path->size(); ++i) {
            const Point64 a = (*path)[i], b = (*path)[(i + 1) % path->size()];
            out.push_back(a);
            const double dx = double(b.x - a.x), dy = double(b.y - a.y), len = std::hypot(dx, dy);
            if (len <= 2 * eps) continue;
            std::set<std::pair<long long, long long>> cells;
            const int steps = static_cast<int>(std::ceil(len / vcell)) + 1;
            for (int k = 0; k <= steps; ++k) {
                const double x = a.x + dx * k / steps, y = a.y + dy * k / steps;
                const long long cx = static_cast<long long>(std::floor(x / vcell)), cy = static_cast<long long>(std::floor(y / vcell));
                for (long long ox = -1; ox <= 1; ++ox)
                    for (long long oy = -1; oy <= 1; ++oy) cells.insert({cx + ox, cy + oy});
            }
            std::vector<std::pair<double, Point64>> on;
            for (const auto& c : cells) {
                const auto it = verts.find(c);
                if (it == verts.end()) continue;
                for (const auto& v : it->second) {
                    if (v == a || v == b) continue;
                    const double t = ((v.x - a.x) * dx + (v.y - a.y) * dy) / (len * len);
                    if (t <= 0 || t >= 1) continue;
                    const double ex = a.x + t * dx - v.x, ey = a.y + t * dy - v.y;
                    if (std::hypot(ex, ey) <= eps) on.push_back({t, v});
                }
            }
            std::sort(on.begin(), on.end(), [](const auto& u, const auto& w) { return u.first < w.first; });
            for (const auto& [t, v] : on)
                if (out.back() != v) out.push_back(v);
        }
        while (out.size() > 1 && out.front() == out.back()) out.pop_back();
        *path = std::move(out);
    }
}

// A 6-node wedge: triangle `tri` of mesh `mesh` between two levels.
struct Wedge {
    int body;
    std::size_t mesh, tri, lo, hi;
};

// A boundary triangle: the top or bottom face of a wedge.
struct Face {
    int boundary;
    std::size_t wedge;
    bool top;
};

// Writes Elmer's mesh files into dir/mesh. Nodes are the (mesh, vertex, level) triples in
// use, numbered from 1 in that order. Returns the nodes of each wedge.
std::vector<std::array<std::size_t, 6>> write_mesh(const std::string& dir, const std::vector<TriMesh>& meshes, const std::vector<double>& levels,
                                                   const std::vector<Wedge>& wedges, const std::vector<Face>& faces, std::size_t& nodes) {
    const std::size_t nl = levels.size();
    std::vector<std::vector<std::size_t>> id(meshes.size());
    for (std::size_t m = 0; m < meshes.size(); ++m) id[m].assign(meshes[m].x.size() * nl, 0);
    for (const auto& w : wedges)
        for (auto v : meshes[w.mesh].triangles[w.tri]) id[w.mesh][v * nl + w.lo] = id[w.mesh][v * nl + w.hi] = 1;
    nodes = 0;
    fs::create_directories(fs::path(dir) / "mesh");
    {
        std::ofstream f(fs::path(dir) / "mesh" / "mesh.nodes");
        char line[160];
        for (std::size_t m = 0; m < meshes.size(); ++m)
            for (std::size_t k = 0; k < id[m].size(); ++k) {
                if (!id[m][k]) continue;
                id[m][k] = ++nodes;
                const std::size_t v = k / nl, l = k % nl;
                std::snprintf(line, sizeof line, "%zu -1 %.12g %.12g %.12g\n", nodes, meshes[m].x[v] * 1e-9, meshes[m].y[v] * 1e-9, levels[l]);
                f << line;
            }
    }
    std::vector<std::array<std::size_t, 6>> elem(wedges.size());
    for (std::size_t w = 0; w < wedges.size(); ++w) {
        const auto& W = wedges[w];
        const auto& tr = meshes[W.mesh].triangles[W.tri];
        for (std::size_t k = 0; k < 3; ++k) {
            elem[w][k] = id[W.mesh][tr[k] * nl + W.lo];
            elem[w][k + 3] = id[W.mesh][tr[k] * nl + W.hi];
        }
    }
    {
        std::ofstream f(fs::path(dir) / "mesh" / "mesh.elements");
        for (std::size_t w = 0; w < wedges.size(); ++w) {
            f << w + 1 << ' ' << wedges[w].body << " 706";
            for (auto v : elem[w]) f << ' ' << v;
            f << '\n';
        }
    }
    {
        std::ofstream f(fs::path(dir) / "mesh" / "mesh.boundary");
        for (std::size_t i = 0; i < faces.size(); ++i) {
            const auto& e = elem[faces[i].wedge];
            f << i + 1 << ' ' << faces[i].boundary << ' ' << faces[i].wedge + 1 << " 0 303";
            if (faces[i].top)
                f << ' ' << e[3] << ' ' << e[4] << ' ' << e[5] << '\n';
            else
                f << ' ' << e[0] << ' ' << e[2] << ' ' << e[1] << '\n';
        }
    }
    {
        std::ofstream f(fs::path(dir) / "mesh" / "mesh.header");
        f << nodes << ' ' << wedges.size() << ' ' << faces.size() << "\n2\n706 " << wedges.size() << "\n303 " << faces.size() << "\n";
    }
    return elem;
}

}  // namespace

ElmerExport export_elmer_dc(const Board& b, const DcSettings& dc, const ElmerSettings& es, const std::string& dir) {
    ElmerExport out;
    const int n = static_cast<int>(b.copper.size());
    if (dc.ports.empty()) throw std::runtime_error("elmer: the DC problem has no ports");
    if (es.copper_divisions < 1) throw std::runtime_error("elmer: copper_divisions at least 1");
    const double tol = 1000.0;

    // Ports and their nets: the nets of the copper the port's area touches.
    std::vector<PortShape> shapes;
    std::set<int> nets;
    for (const auto& p : dc.ports) {
        shapes.push_back(resolve_port(b, p, tol));
        out.ports.push_back(p.name);
        for (const auto& [layer, poly] : shapes.back().by_layer)
            for (const auto& nc : b.copper[static_cast<std::size_t>(layer - 1)].nets)
                if (nc.net != 0 && !intersection(nc.polygons, poly).empty()) nets.insert(nc.net);
    }
    if (nets.empty()) throw std::runtime_error("elmer: the ports touch no copper with a net");

    // Copper of those nets per layer, and their plated holes.
    std::vector<PolygonSet> copper(static_cast<std::size_t>(n) + 1);
    for (const auto& cl : b.copper)
        for (const auto& nc : cl.nets)
            if (nets.count(nc.net)) copper[static_cast<std::size_t>(cl.index)] = unite(copper[static_cast<std::size_t>(cl.index)], nc.polygons);
    struct Hole {
        const Drill* d;
        PolygonSet disk;
    };
    std::vector<Hole> holes;
    for (const auto& d : b.drills)
        if (d.plated && nets.count(d.net) && !d.slot_end && d.from_layer < d.to_layer)
            holes.push_back({&d, union_polygons({circle_path(d.at, 0.5 * d.diameter_nm, 0.25 * es.min_size_nm)})});
    std::size_t slots = 0;
    for (const auto& d : b.drills) slots += d.plated && nets.count(d.net) && d.slot_end;
    if (slots) out.warnings.push_back(std::to_string(slots) + " plated slots are left out");

    // The footprint and its constraints: every layer's copper, the ports, the holes.
    PolygonSet footprint, constraints;
    for (int l = 1; l <= n; ++l) footprint = unite(footprint, copper[static_cast<std::size_t>(l)]);
    for (const auto& h : holes) footprint = unite(footprint, h.disk);
    for (int l = 1; l <= n; ++l) constraints.insert(constraints.end(), copper[static_cast<std::size_t>(l)].begin(), copper[static_cast<std::size_t>(l)].end());
    // A port's outline that follows the copper's own (a contact disc the size of its pad)
    // is left to the copper's: two nearly equal constraints make slivers.
    std::vector<Region> outline_of;
    for (int l = 0; l <= n; ++l) outline_of.emplace_back(copper[static_cast<std::size_t>(l)]);
    const double merge = 0.5 * es.min_size_nm;
    for (const auto& s : shapes)
        for (const auto& [layer, poly] : s.by_layer)
            for (const auto& c : intersection(poly, copper[static_cast<std::size_t>(layer)])) {
                bool follows = true;
                for (const auto& v : c.outer)
                    if (!outline_of[static_cast<std::size_t>(layer)].near_boundary(double(v.x), double(v.y), merge)) {
                        follows = false;
                        break;
                    }
                if (!follows) constraints.push_back(c);
            }
    for (const auto& h : holes) constraints.insert(constraints.end(), h.disk.begin(), h.disk.end());

    {
        std::vector<Path64*> paths;
        for (auto* set : {&footprint, &constraints})
            for (auto& poly : *set) {
                paths.push_back(&poly.outer);
                for (auto& h : poly.holes) paths.push_back(&h);
            }
        clean_paths(paths, std::max(1000.0, 0.02 * es.min_size_nm));
    }

    MeshSettings ms;
    ms.max_edge_nm = es.max_edge_nm;
    ms.min_size_nm = es.min_size_nm;
    ms.grading = es.grading;
    ms.boundary_tolerance_nm = 10.0;  // the outlines are cleaned; simplifying them again would part them
    ms.min_edge_nm = 0.02 * es.min_size_nm;
    std::vector<TriMesh> meshes(footprint.size());
    parallel_for(footprint.size(), [&](std::size_t i) {
        const Rect64 r = bounds(footprint[i]);
        PolygonSet local;
        for (const auto& c : constraints) {
            const Rect64 cr = bounds(c);
            if (cr.right >= r.left && cr.left <= r.right && cr.bottom >= r.top && cr.top <= r.bottom) local.push_back(c);
        }
        meshes[i] = mesh_region(footprint[i], local, ms);
    });

    // Classification of each triangle: copper per layer, plated hole.
    std::vector<Region> in_copper;
    in_copper.reserve(static_cast<std::size_t>(n) + 1);
    for (int l = 0; l <= n; ++l) in_copper.emplace_back(copper[static_cast<std::size_t>(l)]);
    std::vector<Region> in_hole;
    for (const auto& h : holes) in_hole.emplace_back(h.disk);

    // Heights: the faces of the copper layers, with copper_divisions within each.
    const auto zs = stackup_z(b);
    std::vector<double> zb(static_cast<std::size_t>(n) + 1), zt(static_cast<std::size_t>(n) + 1), sigma(static_cast<std::size_t>(n) + 1);
    const double thermal = 1 + dc.temperature_coefficient * (dc.temperature_c - 20.0);
    for (std::size_t i = 0; i < b.stackup.size(); ++i) {
        const auto& l = b.stackup[i];
        if (l.type != LayerType::Copper || l.copper_index < 1 || l.copper_index > n) continue;
        zb[static_cast<std::size_t>(l.copper_index)] = zs[i].bottom_nm * 1e-9;
        zt[static_cast<std::size_t>(l.copper_index)] = zs[i].top_nm * 1e-9;
        const double s20 = dc.copper_conductivity_s_per_m ? *dc.copper_conductivity_s_per_m : l.conductivity_s_per_m.value_or(dc.default_conductivity_s_per_m);
        sigma[static_cast<std::size_t>(l.copper_index)] = s20 / thermal;
    }
    std::vector<double> levels;
    for (int l = 1; l <= n; ++l)
        for (int k = 0; k <= es.copper_divisions; ++k)
            levels.push_back(zb[static_cast<std::size_t>(l)] + (zt[static_cast<std::size_t>(l)] - zb[static_cast<std::size_t>(l)]) * k / es.copper_divisions);
    // The dielectric gaps, which only the holes' columns cross, in steps of at most four
    // minimum sizes: tall thin columns make the linear system ill-conditioned.
    {
        std::vector<double> faces = levels;
        std::sort(faces.begin(), faces.end());
        const double step = 4 * es.min_size_nm * 1e-9;
        for (int l = 1; l < n; ++l) {
            const double lo = zt[static_cast<std::size_t>(l + 1)], hi = zb[static_cast<std::size_t>(l)];
            const int k = static_cast<int>(std::ceil((hi - lo) / step));
            for (int i = 1; i < k; ++i) levels.push_back(lo + (hi - lo) * i / k);
        }
    }
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end(), [](double a, double c) { return std::abs(a - c) < 1e-12; }), levels.end());
    auto level_of = [&](double z) {
        return static_cast<std::size_t>(std::lower_bound(levels.begin(), levels.end(), z - 1e-12) - levels.begin());
    };

    // Elements. Bodies 1..n: the copper of each layer; then one body per hole conductance.
    std::vector<Wedge> wedges;
    std::map<long long, int> hole_body;  // conductivity (rounded) -> body
    std::vector<double> body_sigma(static_cast<std::size_t>(n) + 1);
    for (int l = 1; l <= n; ++l) body_sigma[static_cast<std::size_t>(l)] = sigma[static_cast<std::size_t>(l)];
    const double sigma_plating = dc.plating_conductivity_s_per_m / thermal, tp = dc.plating_thickness_nm * 1e-9;
    std::vector<Face> faces;
    for (std::size_t m = 0; m < meshes.size(); ++m) {
        const auto& tm = meshes[m];
        out.triangles += tm.triangles.size();
        for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
            const auto& tr = tm.triangles[t];
            const double cx = (tm.x[tr[0]] + tm.x[tr[1]] + tm.x[tr[2]]) / 3, cy = (tm.y[tr[0]] + tm.y[tr[1]] + tm.y[tr[2]]) / 3;
            // A plated hole: a column over its span, its conductance spread over the hole.
            int hole = -1;
            for (std::size_t h = 0; h < holes.size() && hole < 0; ++h) {
                const auto& d = *holes[h].d;
                if (std::hypot(cx - double(d.at.x), cy - double(d.at.y)) < 0.5 * d.diameter_nm && in_hole[h].inside(cx, cy)) hole = int(h);
            }
            if (hole >= 0) {
                const auto& d = *holes[static_cast<std::size_t>(hole)].d;
                const double dm = d.diameter_nm * 1e-9;
                const double s_eff = sigma_plating * kPi * tp * (dm + tp) / (kPi * dm * dm / 4);
                const long long key = std::llround(s_eff);
                auto it = hole_body.find(key);
                if (it == hole_body.end()) {
                    body_sigma.push_back(s_eff);
                    it = hole_body.emplace(key, int(body_sigma.size()) - 1).first;
                }
                const std::size_t lo = level_of(zb[static_cast<std::size_t>(d.to_layer)]), hi = level_of(zt[static_cast<std::size_t>(d.from_layer)]);
                for (std::size_t k = lo; k < hi; ++k) wedges.push_back({it->second, m, t, k, k + 1});
                continue;
            }
            for (int l = 1; l <= n; ++l) {
                if (!in_copper[static_cast<std::size_t>(l)].inside(cx, cy)) continue;
                const std::size_t lo = level_of(zb[static_cast<std::size_t>(l)]), hi = level_of(zt[static_cast<std::size_t>(l)]);
                for (std::size_t k = lo; k < hi; ++k) {
                    wedges.push_back({l, m, t, k, k + 1});
                    // Port faces: the outer face of the copper inside a port's area.
                    const bool top_face = k + 1 == hi, bottom_face = k == lo;
                    for (std::size_t p = 0; p < shapes.size(); ++p) {
                        const auto it = shapes[p].by_layer.find(l);
                        if (it == shapes[p].by_layer.end() || !contains(it->second, Point64(std::llround(cx), std::llround(cy)))) continue;
                        const bool use_top = l != n;
                        if ((use_top && top_face) || (!use_top && bottom_face)) faces.push_back({int(p) + 1, wedges.size() - 1, use_top});
                    }
                }
            }
        }
    }
    if (wedges.empty()) throw std::runtime_error("elmer: no elements");
    for (std::size_t p = 0; p < shapes.size(); ++p)
        if (std::none_of(faces.begin(), faces.end(), [&](const Face& f) { return f.boundary == int(p) + 1; }))
            throw std::runtime_error("elmer: port " + dc.ports[p].name + " covers no copper");

    write_mesh(dir, meshes, levels, wedges, faces, out.nodes);
    out.elements = wedges.size();
    out.boundary_elements = faces.size();

    // One case per port: that port at 1 V, the others at 0 V; the ports' currents saved.
    for (std::size_t k = 0; k < dc.ports.size(); ++k) {
        std::ofstream f(fs::path(dir) / ("port_" + std::to_string(k + 1) + ".sif"));
        f << "! kairoq: DC conductances of " << b.name << ", port " << dc.ports[k].name << " at 1 V\n";
        f << "Header\n  Mesh DB \".\" \"mesh\"\nEnd\n\n";
        f << "Simulation\n  Max Output Level = 3\n  Coordinate System = Cartesian\n  Simulation Type = Steady state\n"
             "  Steady State Max Iterations = 1\n  Output Intervals = 0\nEnd\n\n";
        char num[64];
        for (std::size_t body = 1; body < body_sigma.size(); ++body) {
            f << "Body " << body << "\n  Equation = 1\n  Material = " << body << "\nEnd\n";
            if (body <= static_cast<std::size_t>(n)) {
                std::snprintf(num, sizeof num, "%.10g", body_sigma[body]);
                f << "Material " << body << "\n  Electric Conductivity = " << num << "\nEnd\n";
            } else {
                // A hole: its wall conducts around the hole as copper, along it as the
                // plating spread over the hole.
                char lat[64];
                std::snprintf(lat, sizeof lat, "%.10g", sigma_plating);
                std::snprintf(num, sizeof num, "%.10g", body_sigma[body]);
                f << "Material " << body << "\n  Electric Conductivity(3,3) = Real " << lat << " 0 0  0 " << lat << " 0  0 0 " << num << "\nEnd\n";
            }
        }
        f << "\nEquation 1\n  Active Solvers(1) = 1\nEnd\n\n";
        f << "Solver 1\n  Equation = Stat Current Solver\n  Procedure = \"StatCurrentSolve\" \"StatCurrentSolver\"\n  Variable = Potential\n"
             "  Variable DOFs = 1\n  Calculate Volume Current = False\n  Calculate Joule Heating = False\n  Calculate Loads = True\n"
             "  Linear System Solver = Iterative\n  Linear System Iterative Method = CG\n  Linear System Preconditioning = ILU1\n"
             "  Linear System Max Iterations = 10000\n  Linear System Convergence Tolerance = 1.0e-12\n"
             "  Linear System Abort Not Converged = True\nEnd\n\n";
        f << "Solver 2\n  Exec Solver = After Simulation\n  Equation = SaveScalars\n  Procedure = \"SaveData\" \"SaveScalars\"\n"
             "  Filename = \"port_" << k + 1 << ".dat\"\n";
        // One sum of the nodal loads per port, each over the boundaries of its mask.
        for (std::size_t p = 0; p < dc.ports.size(); ++p)
            f << "  Variable " << p + 1 << " = Potential Loads\n  Operator " << p + 1 << " = boundary sum\n  Mask Name " << p + 1
              << " = String \"port " << p + 1 << "\"\n";
        f << "End\n\n";
        for (std::size_t p = 0; p < dc.ports.size(); ++p)
            f << "Boundary Condition " << p + 1 << "\n  Name = \"" << dc.ports[p].name << "\"\n  Target Boundaries(1) = " << p + 1
              << "\n  Potential = " << (p == k ? "1.0" : "0.0") << "\n  port " << p + 1 << " = Logical True\nEnd\n";
    }
    return out;
}

ElmerDcResult run_elmer_dc(const DcSettings& dc, const std::string& dir, const std::string& solver) {
    const auto t0 = Clock::now();
    ElmerDcResult r;
    const std::size_t np = dc.ports.size();
    r.conductance_s.assign(np, std::vector<double>(np, 0.0));
    for (std::size_t k = 0; k < np; ++k) {
        const std::string sif = "port_" + std::to_string(k + 1) + ".sif";
        const std::string cmd = "cd '" + dir + "' && '" + solver + "' " + sif + " > port_" + std::to_string(k + 1) + ".log 2>&1";
        if (std::system(cmd.c_str()) != 0)
            throw std::runtime_error("elmer: " + solver + " failed on " + sif + " (see " + dir + "/port_" + std::to_string(k + 1) + ".log)");
        // The columns named "over bc port <i>" are the sums of the nodal loads on port i: the
        // discrete currents into the copper, exact for the discrete solution.
        const fs::path dat = fs::path(dir) / ("port_" + std::to_string(k + 1) + ".dat");
        std::ifstream f(dat), names(dat.string() + ".names");
        std::vector<double> v;
        double x;
        while (f >> x) v.push_back(x);
        std::string line;
        std::size_t found = 0;
        while (std::getline(names, line)) {
            const auto colon = line.find(':'), over = line.find("over bc port ");
            if (colon == std::string::npos || over == std::string::npos) continue;
            const std::size_t col = std::stoul(line.substr(0, colon)), bc = std::stoul(line.substr(over + 13));
            if (col < 1 || col > v.size() || bc < 1 || bc > np) continue;
            r.conductance_s[bc - 1][k] = v[col - 1];
            ++found;
        }
        if (found != np) throw std::runtime_error("elmer: " + sif + " saved the currents of " + std::to_string(found) + " of " + std::to_string(np) + " ports");
    }
    // The ports' conditions: set voltages, set currents; a group of connected ports without
    // a set voltage takes its first port as 0 V.
    std::vector<int> group(np, -1);
    int ng = 0;
    double gmax = 0;
    for (const auto& row : r.conductance_s)
        for (double g : row) gmax = std::max(gmax, std::abs(g));
    for (std::size_t i = 0; i < np; ++i) {
        if (group[i] >= 0) continue;
        std::vector<std::size_t> stack{i};
        group[i] = ng;
        while (!stack.empty()) {
            const std::size_t a = stack.back();
            stack.pop_back();
            for (std::size_t c = 0; c < np; ++c)
                if (group[c] < 0 && std::abs(r.conductance_s[a][c]) > 1e-9 * gmax) {
                    group[c] = ng;
                    stack.push_back(c);
                }
        }
        ++ng;
    }
    std::vector<bool> fixed(np, false);
    std::vector<double> V(np, 0.0);
    for (std::size_t i = 0; i < np; ++i)
        if (dc.ports[i].voltage_v) {
            fixed[i] = true;
            V[i] = *dc.ports[i].voltage_v;
        }
    for (int g = 0; g < ng; ++g) {
        bool any = false;
        for (std::size_t i = 0; i < np; ++i) any |= group[i] == g && fixed[i];
        if (!any)
            for (std::size_t i = 0; i < np; ++i)
                if (group[i] == g) {
                    fixed[i] = true;
                    V[i] = 0;
                    break;
                }
    }
    std::vector<std::size_t> free;
    for (std::size_t i = 0; i < np; ++i)
        if (!fixed[i]) free.push_back(i);
    const std::size_t nf = free.size();
    std::vector<double> A(nf * nf), rhs(nf);
    for (std::size_t a = 0; a < nf; ++a) {
        const std::size_t i = free[a];
        rhs[a] = dc.ports[i].current_a.value_or(0.0);
        for (std::size_t j = 0; j < np; ++j)
            if (fixed[j]) rhs[a] -= r.conductance_s[i][j] * V[j];
        for (std::size_t c = 0; c < nf; ++c) A[a * nf + c] = r.conductance_s[i][free[c]];
    }
    // Gaussian elimination with partial pivoting.
    for (std::size_t c = 0; c < nf; ++c) {
        std::size_t p = c;
        for (std::size_t q = c + 1; q < nf; ++q)
            if (std::abs(A[q * nf + c]) > std::abs(A[p * nf + c])) p = q;
        for (std::size_t k = 0; k < nf; ++k) std::swap(A[c * nf + k], A[p * nf + k]);
        std::swap(rhs[c], rhs[p]);
        if (A[c * nf + c] == 0) throw std::runtime_error("elmer: a port is not connected to copper");
        for (std::size_t q = c + 1; q < nf; ++q) {
            const double f = A[q * nf + c] / A[c * nf + c];
            for (std::size_t k = c; k < nf; ++k) A[q * nf + k] -= f * A[c * nf + k];
            rhs[q] -= f * rhs[c];
        }
    }
    for (std::size_t c = nf; c-- > 0;) {
        for (std::size_t k = c + 1; k < nf; ++k) rhs[c] -= A[c * nf + k] * rhs[k];
        rhs[c] /= A[c * nf + c];
    }
    for (std::size_t a = 0; a < nf; ++a) V[free[a]] = rhs[a];
    for (std::size_t i = 0; i < np; ++i) {
        DcPortResult p;
        p.name = dc.ports[i].name;
        p.voltage_v = V[i];
        p.voltage_set = dc.ports[i].voltage_v.has_value();
        for (std::size_t j = 0; j < np; ++j) p.current_a += r.conductance_s[i][j] * V[j];
        r.ports.push_back(p);
    }
    r.seconds = since(t0);
    return r;
}

std::string elmer_dc_summary(const ElmerExport& e, const ElmerDcResult* r, const DcResult* native) {
    std::ostringstream o;
    char line[256];
    std::snprintf(line, sizeof line, "elmer model: %zu triangles in the footprint, %zu nodes, %zu wedge elements, %zu port faces\n", e.triangles,
                  e.nodes, e.elements, e.boundary_elements);
    o << line;
    if (r) {
        o << "port                     voltage V (Elmer)   current A (Elmer)";
        if (native) o << "   voltage V (kairoq)   difference mV";
        o << "\n";
        for (std::size_t i = 0; i < r->ports.size(); ++i) {
            const auto& p = r->ports[i];
            std::snprintf(line, sizeof line, "%-24s %17.6f   %17.6f", p.name.c_str(), p.voltage_v, p.current_a);
            o << line;
            if (native)
                for (const auto& q : native->ports)
                    if (q.name == p.name) {
                        std::snprintf(line, sizeof line, "   %20.6f   %13.4f", q.voltage_v, 1e3 * (q.voltage_v - p.voltage_v));
                        o << line;
                    }
            o << "\n";
        }
        std::snprintf(line, sizeof line, "elmer solved %zu cases in %.1f s\n", r->ports.size(), r->seconds);
        o << line;
    }
    for (const auto& w : e.warnings) o << "warning: " << w << "\n";
    return o.str();
}

ElmerExport export_elmer_thermal(const Board& b, const ThermalSettings& ts, const ElmerSettings& es, const std::string& dir) {
    ElmerExport out;
    const int n = static_cast<int>(b.copper.size());
    if (ts.dc) throw std::runtime_error("elmer: Joule heat (thermal.coupling) is not exported");
    if (!ts.heatsinks.empty()) throw std::runtime_error("elmer: heat sinks are not exported");
    for (const auto& h : ts.sources)
        if (h.temperature_c || h.junction_k_per_w || !h.power_w)
            throw std::runtime_error("elmer: source " + h.area.name + ": only a set power (power_w) is exported");
    if (ts.sources.empty()) throw std::runtime_error("elmer: the thermal problem has no sources");
    if (es.copper_divisions < 1) throw std::runtime_error("elmer: copper_divisions at least 1");
    if (ts.top.emissivity > 0 || ts.bottom.emissivity > 0) out.warnings.push_back("radiation is left out");
    if (ts.transient) out.warnings.push_back("the steady state is exported");
    const double tol = b.source.arc_tolerance_nm > 0 ? b.source.arc_tolerance_nm : 1000.0;

    // The board without its non-plated holes, the copper of every layer and the plated holes.
    Paths64 npth;
    for (const auto& d : b.drills)
        if (!d.plated) npth.push_back(hole_outline(d, tol));
    PolygonSet domain = npth.empty() ? b.outline : difference(b.outline, union_polygons(npth));
    std::vector<PolygonSet> copper(static_cast<std::size_t>(n) + 1);
    for (const auto& cl : b.copper) {
        Paths64 all;
        for (const auto& nc : cl.nets) {
            const auto p = to_paths(nc.polygons);
            all.insert(all.end(), p.begin(), p.end());
        }
        copper[static_cast<std::size_t>(cl.index)] = intersection(union_polygons(all), domain);
    }
    struct Hole {
        const Drill* d;
        PolygonSet disk;
    };
    std::vector<Hole> holes;
    std::size_t slots = 0;
    for (const auto& d : b.drills) {
        if (!d.plated || d.from_layer >= d.to_layer) continue;
        if (d.slot_end) {
            ++slots;
            continue;
        }
        holes.push_back({&d, intersection(union_polygons({circle_path(d.at, 0.5 * d.diameter_nm, 0.25 * es.min_size_nm)}), domain)});
    }
    if (slots) out.warnings.push_back(std::to_string(slots) + " plated slots are left out");
    std::vector<std::map<int, PolygonSet>> source_area;
    for (const auto& h : ts.sources) {
        std::map<int, PolygonSet> by_layer;
        for (const auto& [layer, poly] : resolve_port(b, h.area, tol).by_layer) {
            auto c = intersection(intersection(poly, domain), copper[static_cast<std::size_t>(layer)]);
            if (!c.empty()) by_layer[layer] = std::move(c);
        }
        if (by_layer.empty()) throw std::runtime_error("elmer: source " + h.area.name + " covers no copper");
        source_area.push_back(std::move(by_layer));
        out.ports.push_back(h.area.name);
    }

    // Constraints: every layer's copper, the holes, the sources' areas where they do not
    // follow the copper's outline.
    PolygonSet constraints;
    for (int l = 1; l <= n; ++l) constraints.insert(constraints.end(), copper[static_cast<std::size_t>(l)].begin(), copper[static_cast<std::size_t>(l)].end());
    std::vector<Region> outline_of;
    for (int l = 0; l <= n; ++l) outline_of.emplace_back(copper[static_cast<std::size_t>(l)]);
    const double merge = 0.5 * es.min_size_nm;
    for (const auto& area : source_area)
        for (const auto& [layer, poly] : area)
            for (const auto& c : poly) {
                bool follows = true;
                for (const auto& v : c.outer)
                    if (!outline_of[static_cast<std::size_t>(layer)].near_boundary(double(v.x), double(v.y), merge)) {
                        follows = false;
                        break;
                    }
                if (!follows) constraints.push_back(c);
            }
    for (const auto& h : holes) constraints.insert(constraints.end(), h.disk.begin(), h.disk.end());
    {
        std::vector<Path64*> paths;
        for (auto* set : {&domain, &constraints})
            for (auto& poly : *set) {
                paths.push_back(&poly.outer);
                for (auto& h : poly.holes) paths.push_back(&h);
            }
        clean_paths(paths, std::max(1000.0, 0.02 * es.min_size_nm));
    }
    MeshSettings ms;
    ms.max_edge_nm = es.max_edge_nm;
    ms.min_size_nm = es.min_size_nm;
    ms.grading = es.grading;
    ms.boundary_tolerance_nm = 10.0;
    ms.min_edge_nm = 0.02 * es.min_size_nm;
    std::vector<TriMesh> meshes(domain.size());
    parallel_for(domain.size(), [&](std::size_t i) {
        const Rect64 r = bounds(domain[i]);
        PolygonSet local;
        for (const auto& c : constraints)
            if (bounds(c).Intersects(r)) local.push_back(c);
        meshes[i] = mesh_region(domain[i], local, ms);
    });

    // Heights: every layer between the outer faces of the outer copper, copper layers in
    // copper_divisions, dielectric layers in steps of at most two minimum sizes. Solder mask
    // is a series resistance to the surface, as in the layered model.
    const auto zs = stackup_z(b);
    struct Slab {
        double lo, hi;
        int copper;     // copper layer index, or 0
        double kxy, kz; // dielectric
    };
    std::vector<Slab> slabs;
    double r_top = 0, r_bottom = 0;
    bool seen = false;
    std::vector<double> levels;
    for (std::size_t i = 0; i < b.stackup.size(); ++i) {
        const auto& l = b.stackup[i];
        const double lo = zs[i].bottom_nm * 1e-9, hi = zs[i].top_nm * 1e-9;
        if (hi - lo <= 0) continue;
        if (l.type == LayerType::Copper) {
            seen = true;
            slabs.push_back({lo, hi, l.copper_index, 0, 0});
            for (int k = 0; k <= es.copper_divisions; ++k) levels.push_back(lo + (hi - lo) * k / es.copper_divisions);
            continue;
        }
        const bool mask = l.type == LayerType::SolderMask;
        double kxy = mask ? ts.mask_k : ts.dielectric_k_xy, kz = mask ? ts.mask_k : ts.dielectric_k_z;
        if (auto it = ts.layer_k.find(l.name); it != ts.layer_k.end()) {
            kxy = it->second[0];
            kz = it->second[1];
        }
        if (mask) {
            (seen ? r_bottom : r_top) += (hi - lo) / kz;
            continue;
        }
        slabs.push_back({lo, hi, 0, kxy, kz});
        const int k = std::max(1, static_cast<int>(std::ceil((hi - lo) / (2 * es.min_size_nm * 1e-9))));
        for (int j = 0; j <= k; ++j) levels.push_back(lo + (hi - lo) * j / k);
    }
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end(), [](double a, double c) { return std::abs(a - c) < 1e-12; }), levels.end());
    auto level_of = [&](double z) {
        return static_cast<std::size_t>(std::lower_bound(levels.begin(), levels.end(), z - 1e-12) - levels.begin());
    };
    std::vector<double> zb(static_cast<std::size_t>(n) + 1), zt(static_cast<std::size_t>(n) + 1);
    std::vector<std::size_t> copper_slab(static_cast<std::size_t>(n) + 1);
    for (std::size_t i = 0; i < slabs.size(); ++i)
        if (slabs[i].copper) {
            zb[static_cast<std::size_t>(slabs[i].copper)] = slabs[i].lo;
            zt[static_cast<std::size_t>(slabs[i].copper)] = slabs[i].hi;
            copper_slab[static_cast<std::size_t>(slabs[i].copper)] = i;
        }

    // Bodies: one per material, a source's copper its own body with the source's heat.
    struct Body {
        double kxy, kz;
        int source;  // or -1
    };
    std::vector<Body> bodies{{}};  // from 1
    std::map<std::tuple<long long, long long, int>, int> body_of;
    auto body = [&](double kxy, double kz, int source) {
        const auto key = std::make_tuple(std::llround(kxy * 1e6), std::llround(kz * 1e6), source);
        auto it = body_of.find(key);
        if (it == body_of.end()) {
            bodies.push_back({kxy, kz, source});
            it = body_of.emplace(key, int(bodies.size()) - 1).first;
        }
        return it->second;
    };
    const double kcu = ts.copper_k, tp = ts.plating_thickness_nm * 1e-9;
    // Material between the copper of a layer: inner layers are embedded in the dielectric
    // below them; outer ones in air (no element) or mask.
    std::vector<int> fill(static_cast<std::size_t>(n) + 1, 0);
    for (int l = 1; l <= n; ++l) {
        const std::size_t i = copper_slab[static_cast<std::size_t>(l)];
        if (l == 1 || l == n) {
            if (ts.outer_gaps_mask) fill[static_cast<std::size_t>(l)] = body(ts.mask_k, ts.mask_k, -1);
        } else if (i + 1 < slabs.size() && !slabs[i + 1].copper) {
            fill[static_cast<std::size_t>(l)] = body(slabs[i + 1].kxy, slabs[i + 1].kz, -1);
        }
    }

    std::vector<Region> in_copper;
    in_copper.reserve(static_cast<std::size_t>(n) + 1);
    for (int l = 0; l <= n; ++l) in_copper.emplace_back(copper[static_cast<std::size_t>(l)]);
    std::vector<std::vector<std::pair<int, Region>>> in_source(source_area.size());
    for (std::size_t s_ = 0; s_ < source_area.size(); ++s_)
        for (const auto& [layer, poly] : source_area[s_]) in_source[s_].emplace_back(layer, Region(poly));
    // Holes on a 1 mm grid of their centres.
    std::map<std::pair<long long, long long>, std::vector<std::size_t>> hole_cells;
    const double hcell = 1e6;
    for (std::size_t h = 0; h < holes.size(); ++h)
        hole_cells[{std::llround(std::floor(holes[h].d->at.x / hcell)), std::llround(std::floor(holes[h].d->at.y / hcell))}].push_back(h);
    std::vector<Region> in_hole;
    for (const auto& h : holes) in_hole.emplace_back(h.disk);

    std::vector<Wedge> wedges;
    std::vector<Face> faces;
    std::vector<double> source_volume(source_area.size(), 0.0);
    std::vector<std::vector<std::size_t>> source_wedges(source_area.size());
    for (std::size_t m = 0; m < meshes.size(); ++m) {
        const auto& tm = meshes[m];
        out.triangles += tm.triangles.size();
        for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
            const auto& tr = tm.triangles[t];
            const double cx = (tm.x[tr[0]] + tm.x[tr[1]] + tm.x[tr[2]]) / 3, cy = (tm.y[tr[0]] + tm.y[tr[1]] + tm.y[tr[2]]) / 3;
            const double area = triangle_area(tm, t) * 1e-18;
            int hole = -1;
            const long long hx = std::llround(std::floor(cx / hcell)), hy = std::llround(std::floor(cy / hcell));
            for (long long ox = -1; ox <= 1 && hole < 0; ++ox)
                for (long long oy = -1; oy <= 1 && hole < 0; ++oy) {
                    const auto it = hole_cells.find({hx + ox, hy + oy});
                    if (it == hole_cells.end()) continue;
                    for (auto h : it->second) {
                        const auto& d = *holes[h].d;
                        if (std::hypot(cx - double(d.at.x), cy - double(d.at.y)) < 0.5 * d.diameter_nm && in_hole[h].inside(cx, cy)) {
                            hole = int(h);
                            break;
                        }
                    }
                }
            std::vector<int> layer_body(static_cast<std::size_t>(n) + 1, 0);
            for (int l = 1; l <= n; ++l) {
                if (!in_copper[static_cast<std::size_t>(l)].inside(cx, cy)) {
                    layer_body[static_cast<std::size_t>(l)] = fill[static_cast<std::size_t>(l)];
                    continue;
                }
                int src = -1;
                for (std::size_t s_ = 0; s_ < in_source.size() && src < 0; ++s_)
                    for (const auto& [layer, reg] : in_source[s_])
                        if (layer == l && reg.inside(cx, cy)) src = int(s_);
                layer_body[static_cast<std::size_t>(l)] = body(kcu, kcu, src);
            }
            const std::size_t first = wedges.size();
            for (const auto& sl : slabs) {
                const std::size_t lo = level_of(sl.lo), hi = level_of(sl.hi);
                int bd;
                const Drill* d = hole >= 0 ? holes[static_cast<std::size_t>(hole)].d : nullptr;
                if (d && sl.lo >= zb[static_cast<std::size_t>(d->to_layer)] - 1e-12 && sl.hi <= zt[static_cast<std::size_t>(d->from_layer)] + 1e-12) {
                    // A plated hole: the wall conducts around the hole as copper, along it as
                    // the plating spread over the hole.
                    const double dm = d->diameter_nm * 1e-9;
                    bd = body(kcu, kcu * kPi * tp * (dm + tp) / (kPi * dm * dm / 4), -1);
                } else {
                    bd = sl.copper ? layer_body[static_cast<std::size_t>(sl.copper)] : body(sl.kxy, sl.kz, -1);
                }
                if (!bd) continue;
                for (std::size_t k = lo; k < hi; ++k) {
                    wedges.push_back({bd, m, t, k, k + 1});
                    if (const int src = bodies[static_cast<std::size_t>(bd)].source; src >= 0) {
                        source_volume[static_cast<std::size_t>(src)] += area * (levels[k + 1] - levels[k]);
                        source_wedges[static_cast<std::size_t>(src)].push_back(wedges.size() - 1);
                    }
                }
            }
            if (wedges.size() == first) continue;
            // Faces: the column's top and bottom; boundary 1 the top, 2 the bottom.
            faces.push_back({1, wedges.size() - 1, true});
            faces.push_back({2, first, false});
        }
    }
    if (wedges.empty()) throw std::runtime_error("elmer: no elements");
    const auto elem = write_mesh(dir, meshes, levels, wedges, faces, out.nodes);
    out.elements = wedges.size();
    out.boundary_elements = faces.size();

    // The probes: each source's nodes with their share of its volume, each copper layer's
    // nodes.
    {
        std::ofstream f(fs::path(dir) / "probes.txt");
        f.precision(12);
        for (std::size_t s_ = 0; s_ < source_area.size(); ++s_) {
            std::map<std::size_t, double> w;
            for (auto k : source_wedges[s_]) {
                const auto& W = wedges[k];
                const double v = triangle_area(meshes[W.mesh], W.tri) * 1e-18 * (levels[W.hi] - levels[W.lo]);
                for (auto node : elem[k]) w[node] += v / 6;
            }
            f << "source " << ts.sources[s_].area.name << ' ' << w.size() << '\n';
            for (const auto& [node, v] : w) f << node << ' ' << v / source_volume[s_] << '\n';
        }
        for (int l = 1; l <= n; ++l) {
            std::set<std::size_t> nodes;
            const std::size_t lo = level_of(zb[static_cast<std::size_t>(l)]), hi = level_of(zt[static_cast<std::size_t>(l)]);
            for (std::size_t k = 0; k < wedges.size(); ++k)
                if (wedges[k].lo >= lo && wedges[k].hi <= hi)
                    for (auto node : elem[k]) nodes.insert(node);
            f << "layer " << l << ' ' << nodes.size() << '\n';
            for (auto node : nodes) f << node << '\n';
        }
    }

    // The case: steady conduction, convection through the mask on both faces.
    std::ofstream f(fs::path(dir) / "thermal.sif");
    char num[200];
    f << "! kairoq: steady temperature of " << b.name << "\n";
    f << "Header\n  Mesh DB \".\" \"mesh\"\n  Results Directory \".\"\nEnd\n\n";
    f << "Simulation\n  Max Output Level = 3\n  Coordinate System = Cartesian\n  Simulation Type = Steady state\n"
         "  Steady State Max Iterations = 1\n  Output Intervals = 0\nEnd\n\n";
    for (std::size_t k = 1; k < bodies.size(); ++k) {
        f << "Body " << k << "\n  Equation = 1\n  Material = " << k << "\n";
        if (bodies[k].source >= 0) f << "  Body Force = " << bodies[k].source + 1 << "\n";
        f << "End\n";
        std::snprintf(num, sizeof num, "Real %.10g 0 0  0 %.10g 0  0 0 %.10g", bodies[k].kxy, bodies[k].kxy, bodies[k].kz);
        f << "Material " << k << "\n  Density = 1.0\n  Heat Conductivity(3,3) = " << num << "\nEnd\n";
    }
    for (std::size_t s_ = 0; s_ < ts.sources.size(); ++s_) {
        // Heat per mass at a density of 1: the power over the source's copper volume.
        std::snprintf(num, sizeof num, "%.10g", *ts.sources[s_].power_w / source_volume[s_]);
        f << "Body Force " << s_ + 1 << "\n  Name = \"" << ts.sources[s_].area.name << "\"\n  Heat Source = " << num << "\nEnd\n";
    }
    f << "\nEquation 1\n  Active Solvers(1) = 1\nEnd\n\n";
    // IDR(s) with ILU1: on these stiff systems (thin copper, conductivities 1e3 apart) about
    // fifteen times faster than conjugate gradients with the same preconditioner.
    f << "Solver 1\n  Equation = Heat Equation\n  Procedure = \"HeatSolve\" \"HeatSolver\"\n  Variable = Temperature\n"
         "  Variable DOFs = 1\n  Nonlinear System Max Iterations = 1\n"
         "  Linear System Solver = Iterative\n  Linear System Iterative Method = IDRS\n  Linear System Preconditioning = ILU1\n"
         "  Linear System Max Iterations = 20000\n  Linear System Convergence Tolerance = 1.0e-12\n"
         "  Linear System Abort Not Converged = True\nEnd\n\n";
    f << "Solver 2\n  Exec Solver = After Simulation\n  Equation = ResultOutput\n  Procedure = \"ResultOutputSolve\" \"ResultOutputSolver\"\n"
         "  Output File Name = \"thermal\"\n  Output Format = String \"gmsh\"\n  Scalar Field 1 = Temperature\nEnd\n\n";
    const SurfaceSpec* side[2] = {&ts.top, &ts.bottom};
    const double r_mask[2] = {r_top, r_bottom};
    for (int k = 0; k < 2; ++k) {
        const double h = side[k]->convection_w_m2k > 0 ? 1.0 / (1.0 / side[k]->convection_w_m2k + r_mask[k]) : 0.0;
        std::snprintf(num, sizeof num, "  Heat Transfer Coefficient = %.10g\n  External Temperature = %.10g\n", h, ts.ambient_c + 273.15);
        f << "Boundary Condition " << k + 1 << "\n  Name = \"" << (k ? "bottom" : "top") << "\"\n  Target Boundaries(1) = " << k + 1 << "\n" << num
          << "End\n";
    }
    return out;
}

ElmerThermalResult run_elmer_thermal(const std::string& dir, const std::string& solver) {
    const auto t0 = Clock::now();
    const std::string cmd = "cd '" + dir + "' && '" + solver + "' thermal.sif > thermal.log 2>&1";
    if (std::system(cmd.c_str()) != 0) throw std::runtime_error("elmer: " + solver + " failed on thermal.sif (see " + dir + "/thermal.log)");
    // The temperatures (K) at the nodes, from the gmsh file.
    fs::path msh;
    for (const auto& e : fs::recursive_directory_iterator(dir))
        if (e.path().extension() == ".msh" && e.path().stem().string().rfind("thermal", 0) == 0) msh = e.path();
    if (msh.empty()) throw std::runtime_error("elmer: no thermal*.msh in " + dir);
    std::ifstream g(msh);
    std::string line;
    std::vector<double> T;
    while (std::getline(g, line)) {
        if (line != "$NodeData") continue;
        int tags = 0;
        g >> tags;
        std::getline(g, line);
        for (int i = 0; i < tags; ++i) std::getline(g, line);
        g >> tags;
        double x;
        for (int i = 0; i < tags; ++i) g >> x;
        g >> tags;
        std::size_t count = 0;
        for (int i = 0; i < tags; ++i) {
            long v;
            g >> v;
            if (i == 2) count = static_cast<std::size_t>(v);
        }
        T.assign(count + 1, std::nan(""));
        for (std::size_t i = 0; i < count; ++i) {
            std::size_t id;
            double v;
            g >> id >> v;
            if (id < T.size()) T[id] = v - 273.15;
        }
        break;
    }
    if (T.size() < 2) throw std::runtime_error("elmer: no temperatures in " + msh.string());
    ElmerThermalResult r;
    std::ifstream p(fs::path(dir) / "probes.txt");
    std::string kind, name;
    std::size_t count;
    while (p >> kind) {
        if (kind == "source") {
            p >> name >> count;
            ElmerThermalResult::Source s;
            s.name = name;
            s.max_c = -INFINITY;
            for (std::size_t i = 0; i < count; ++i) {
                std::size_t node;
                double w;
                p >> node >> w;
                if (node >= T.size()) throw std::runtime_error("elmer: probe node " + std::to_string(node) + " has no temperature");
                s.mean_c += w * T[node];
                s.max_c = std::max(s.max_c, T[node]);
            }
            r.sources.push_back(s);
        } else if (kind == "layer") {
            int l;
            p >> l >> count;
            if (r.layer_max_c.size() <= static_cast<std::size_t>(l)) r.layer_max_c.resize(static_cast<std::size_t>(l) + 1, -INFINITY);
            for (std::size_t i = 0; i < count; ++i) {
                std::size_t node;
                p >> node;
                if (node < T.size()) r.layer_max_c[static_cast<std::size_t>(l)] = std::max(r.layer_max_c[static_cast<std::size_t>(l)], T[node]);
            }
        }
    }
    r.seconds = since(t0);
    return r;
}

std::string elmer_thermal_summary(const ElmerExport& e, const ElmerThermalResult* r, const ThermalResult* native) {
    std::ostringstream o;
    char line[256];
    std::snprintf(line, sizeof line, "elmer model: %zu triangles, %zu nodes, %zu wedge elements, %zu surface faces\n", e.triangles, e.nodes,
                  e.elements, e.boundary_elements);
    o << line;
    if (r) {
        o << "source                   mean C (Elmer)   max C (Elmer)";
        if (native) o << "   mean C (kairoq)   max C (kairoq)";
        o << "\n";
        for (const auto& s : r->sources) {
            std::snprintf(line, sizeof line, "%-24s %14.4f   %13.4f", s.name.c_str(), s.mean_c, s.max_c);
            o << line;
            if (native)
                for (const auto& q : native->sources)
                    if (q.name == s.name) {
                        std::snprintf(line, sizeof line, "   %17.4f   %16.4f", q.mean_c, q.max_c);
                        o << line;
                    }
            o << "\n";
        }
        o << "copper layer             max C (Elmer)";
        if (native) o << "   max C (kairoq)";
        o << "\n";
        for (std::size_t l = 1; l < r->layer_max_c.size(); ++l) {
            std::snprintf(line, sizeof line, "L%-23zu %13.4f", l, r->layer_max_c[l]);
            o << line;
            if (native)
                for (const auto& lf : native->layers)
                    if (lf.layer == int(l) && !lf.temperature_c.empty()) {
                        std::snprintf(line, sizeof line, "   %16.4f", *std::max_element(lf.temperature_c.begin(), lf.temperature_c.end()));
                        o << line;
                    }
            o << "\n";
        }
        std::snprintf(line, sizeof line, "elmer solved in %.1f s\n", r->seconds);
        o << line;
    }
    for (const auto& w : e.warnings) o << "warning: " << w << "\n";
    return o.str();
}

}  // namespace kairoq
