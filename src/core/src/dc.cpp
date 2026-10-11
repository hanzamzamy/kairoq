#include "kairoq/dc.hpp"

#include "fem_util.hpp"


#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace kairoq {

namespace {

using namespace detail;

// A connected piece of copper of one net on one layer.
struct Island {
    int layer = 0;  // 1-based
    int net = 0;
    const Polygon* polygon = nullptr;
    Rect64 box;
};

}  // namespace

DcResult solve_dc(const Board& b, const DcSettings& s) {
    DcResult res;
    if (s.ports.empty()) throw std::runtime_error("no ports given");
    const int nlayers = static_cast<int>(b.copper.size());
    const double tol = b.source.arc_tolerance_nm > 0 ? b.source.arc_tolerance_nm : 1000.0;
    const double thermal = 1.0 + s.temperature_coefficient * (s.temperature_c - 20.0);

    // Layer properties.
    const auto zs = stackup_z(b);
    std::vector<double> thickness(nlayers + 1, 0.0), sigma(nlayers + 1, 0.0), zmid(nlayers + 1, 0.0);
    for (std::size_t i = 0; i < b.stackup.size(); ++i) {
        const auto& l = b.stackup[i];
        if (l.type != LayerType::Copper || l.copper_index < 1 || l.copper_index > nlayers) continue;
        thickness[l.copper_index] = l.thickness_nm;
        const double s20 = s.copper_conductivity_s_per_m ? *s.copper_conductivity_s_per_m
                                                          : l.conductivity_s_per_m.value_or(s.default_conductivity_s_per_m);
        sigma[l.copper_index] = s20 / thermal;
        zmid[l.copper_index] = 0.5 * (zs[i].bottom_nm + zs[i].top_nm);
    }
    for (int l = 1; l <= nlayers; ++l)
        if (thickness[l] <= 0) throw std::runtime_error("copper layer " + std::to_string(l) + " has no thickness in the stackup");

    // Islands and their grid index per layer.
    std::vector<Island> islands;
    std::vector<std::vector<std::size_t>> layer_islands(nlayers + 1);
    for (const auto& cl : b.copper)
        for (const auto& nc : cl.nets)
            for (const auto& p : nc.polygons) {
                layer_islands[cl.index].push_back(islands.size());
                islands.push_back({cl.index, nc.net, &p, bounds(p)});
            }
    std::vector<std::unique_ptr<GridIndex>> index(nlayers + 1);
    for (int l = 1; l <= nlayers; ++l) {
        std::vector<Rect64> boxes;
        for (auto i : layer_islands[l]) boxes.push_back(islands[i].box);
        index[l] = std::make_unique<GridIndex>(boxes, 1.0e6);
    }
    // Island at a point; a net other than 0 restricts the search to that net.
    auto island_at = [&](int layer, Point64 p, int net) -> long {
        for (auto k : index[layer]->query(p)) {
            const auto i = layer_islands[layer][k];
            if ((net == 0 || islands[i].net == net) && contains(*islands[i].polygon, p)) return static_cast<long>(i);
        }
        return -1;
    };

    // Ports.
    std::vector<PortShape> port_shapes;
    std::vector<std::vector<std::size_t>> port_islands(s.ports.size());
    for (std::size_t p = 0; p < s.ports.size(); ++p) {
        const auto& spec = s.ports[p];
        if (spec.voltage_v.has_value() == spec.current_a.has_value())
            throw std::runtime_error("port " + spec.name + ": give either voltage_v or current_a");
        port_shapes.push_back(resolve_port(b, spec, tol));
        for (const auto& [layer, shape] : port_shapes.back().by_layer) {
            const Rect64 box = bounds(shape);
            for (auto i : layer_islands[layer])
                if (boxes_overlap(islands[i].box, box) && !intersection({*islands[i].polygon}, shape).empty())
                    port_islands[p].push_back(i);
        }
        if (port_islands[p].empty()) throw std::runtime_error("port " + spec.name + " touches no copper");
    }

    // Plated holes: the islands they connect, per layer.
    struct HoleLink {
        std::size_t drill;
        int layer;
        std::size_t island;
    };
    std::vector<HoleLink> links;
    std::vector<std::vector<std::size_t>> drill_links(b.drills.size());
    for (std::size_t d = 0; d < b.drills.size(); ++d) {
        const auto& dr = b.drills[d];
        if (!dr.plated) continue;
        const int from = std::clamp(dr.from_layer, 1, nlayers), to = std::clamp(dr.to_layer, 1, nlayers);
        for (int l = from; l <= to; ++l) {
            // A hole with a net joins copper of that net only; on a board with shorts it
            // must not couple two nets.
            const long i = island_at(l, dr.at, dr.net);
            if (i < 0) continue;
            drill_links[d].push_back(links.size());
            links.push_back({d, l, static_cast<std::size_t>(i)});
        }
    }

    // Islands connected to a port.
    UnionFind uf_islands(islands.size());
    for (const auto& dl : drill_links)
        for (std::size_t k = 1; k < dl.size(); ++k) uf_islands.unite(links[dl[k - 1]].island, links[dl[k]].island);
    for (const auto& pi : port_islands)
        for (std::size_t k = 1; k < pi.size(); ++k) uf_islands.unite(pi[k - 1], pi[k]);
    std::set<std::size_t> active_roots;
    for (const auto& pi : port_islands) active_roots.insert(uf_islands.find(pi[0]));
    std::vector<std::size_t> active;
    for (std::size_t i = 0; i < islands.size(); ++i)
        if (active_roots.count(uf_islands.find(i))) active.push_back(i);
    res.islands = active.size();

    // Mesh the active islands. Plated holes and ports inside an island are constraints
    // (their areas become equipotential); non-plated holes are removed.
    std::vector<PolygonSet> port_by_island(islands.size());
    for (std::size_t p = 0; p < s.ports.size(); ++p)
        for (auto i : port_islands[p])
            port_by_island[i] = unite(port_by_island[i], port_shapes[p].by_layer.at(islands[i].layer));
    std::vector<std::vector<std::size_t>> island_links(islands.size());
    for (std::size_t k = 0; k < links.size(); ++k) island_links[links[k].island].push_back(k);

    std::vector<PolygonSet> hole_shapes(links.size());
    for (std::size_t k = 0; k < links.size(); ++k) hole_shapes[k] = union_polygons({hole_outline(b.drills[links[k].drill], tol)});
    std::vector<PolygonSet> regions(islands.size()), constraints(islands.size());
    parallel_for(active.size(), [&](std::size_t a) {
        const std::size_t i = active[a];
        const Island& isl = islands[i];
        regions[i] = PolygonSet{*isl.polygon};
        Paths64 npth;
        for (const auto& d : b.drills)
            if (!d.plated) {
                const Path64 h = hole_outline(d, tol);
                if (boxes_overlap(Clipper2Lib::GetBounds(h), isl.box)) npth.push_back(h);
            }
        if (!npth.empty()) regions[i] = difference(regions[i], union_polygons(npth));
        constraints[i] = port_by_island[i];
        for (auto k : island_links[i]) constraints[i].insert(constraints[i].end(), hole_shapes[k].begin(), hole_shapes[k].end());
    });

    // Adaptive passes: mesh, solve, estimate the error, refine where it is large.
    std::vector<std::vector<std::array<double, 3>>> seeds(islands.size());
    if (s.seeds.size() == islands.size()) seeds = s.seeds;
    const std::vector<std::string> setup_warnings = res.warnings;
    for (int pass = 0;; ++pass) {
        res.warnings = setup_warnings;
        const auto t_mesh = Clock::now();
        std::vector<std::vector<TriMesh>> meshes(islands.size());
        std::atomic<std::size_t> unrefined{0};
        parallel_for(active.size(), [&](std::size_t a) {
            const std::size_t i = active[a];
            for (const auto& part : regions[i]) {
                meshes[i].push_back(mesh_region(part, constraints[i], s.mesh, seeds[i]));
                unrefined += meshes[i].back().unrefined;
            }
        });
        if (unrefined) res.warnings.push_back(std::to_string(unrefined.load()) + " mesh refinements could not be performed");
        res.seconds_mesh += since(t_mesh);
        const auto t_assemble = Clock::now();

        // Global nodes and the equipotential groups (ports, hole walls).
        struct Part {
            std::size_t island;
            const TriMesh* mesh;
            std::size_t offset;
        };
        std::vector<Part> parts;
        std::size_t nnodes = 0;
        for (auto i : active)
            for (const auto& m : meshes[i]) {
                parts.push_back({i, &m, nnodes});
                nnodes += m.x.size();
            }
        // Sheet conductance of each triangle (S); with a temperature field, the copper
        // conductivity at the temperature of the triangle's centroid.
        const auto sigma_at = [&](int layer, double x, double y) {
            if (!s.temperature_field) return sigma[layer];
            const double s20 = sigma[layer] * thermal;
            return s20 / (1.0 + s.temperature_coefficient * (s.temperature_field(layer, x, y) - 20.0));
        };
        std::vector<std::vector<double>> gsheet(parts.size());
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const TriMesh& m = *parts[k].mesh;
            const int layer = islands[parts[k].island].layer;
            gsheet[k].resize(m.triangles.size());
            for (std::size_t t = 0; t < m.triangles.size(); ++t) {
                const auto& tr = m.triangles[t];
                const double cx = (m.x[tr[0]] + m.x[tr[1]] + m.x[tr[2]]) / 3, cy = (m.y[tr[0]] + m.y[tr[1]] + m.y[tr[2]]) / 3;
                gsheet[k][t] = sigma_at(layer, cx, cy) * thickness[layer] * 1e-9;
            }
        }
        UnionFind uf(nnodes);
        std::vector<long> port_node(s.ports.size(), -1);
        std::vector<std::size_t> port_count(s.ports.size(), 0);
        std::vector<long> link_node(links.size(), -1);
        {
            std::vector<std::vector<std::size_t>> parts_of(islands.size());
            for (std::size_t k = 0; k < parts.size(); ++k) parts_of[parts[k].island].push_back(k);
            for (std::size_t k = 0; k < parts.size(); ++k) {
                const Part& pt = parts[k];
                NodeGrid grid(*pt.mesh, 2.0e5);
                for (std::size_t p = 0; p < s.ports.size(); ++p) {
                    if (std::find(port_islands[p].begin(), port_islands[p].end(), pt.island) == port_islands[p].end()) continue;
                    for (auto n : grid.inside(port_shapes[p].by_layer.at(islands[pt.island].layer))) {
                        const std::size_t g = pt.offset + n;
                        if (port_node[p] < 0) port_node[p] = static_cast<long>(g);
                        else uf.unite(g, static_cast<std::size_t>(port_node[p]));
                        ++port_count[p];
                    }
                }
                for (auto l : island_links[pt.island]) {
                    for (auto n : grid.inside(hole_shapes[l])) {
                        const std::size_t g = pt.offset + n;
                        if (link_node[l] < 0) link_node[l] = static_cast<long>(g);
                        else uf.unite(g, static_cast<std::size_t>(link_node[l]));
                    }
                    if (link_node[l] < 0) {
                        // No node on the hole wall: the nearest node within two radii.
                        const auto& dr = b.drills[links[l].drill];
                        double best = dr.diameter_nm;
                        const TriMesh& m = *pt.mesh;
                        for (std::size_t n = 0; n < m.x.size(); ++n) {
                            const double dd = std::hypot(m.x[n] - double(dr.at.x), m.y[n] - double(dr.at.y));
                            if (dd <= best) {
                                best = dd;
                                link_node[l] = static_cast<long>(pt.offset + n);
                            }
                        }
                    }
                }
            }
        }
        for (std::size_t p = 0; p < s.ports.size(); ++p)
            if (port_node[p] < 0) throw std::runtime_error("port " + s.ports[p].name + " covers no mesh node");

        // Unknowns: one per group of tied nodes.
        std::vector<std::size_t> dof(nnodes);
        std::unordered_map<std::size_t, std::size_t> root_dof;
        for (std::size_t g = 0; g < nnodes; ++g) {
            const auto r = uf.find(g);
            auto [it, inserted] = root_dof.emplace(r, root_dof.size());
            dof[g] = it->second;
        }
        const std::size_t ndof = root_dof.size();

        // Element and barrel conductances.
        std::vector<Triplet> trip;
        std::size_t ntri = 0;
        for (const auto& pt : parts) ntri += pt.mesh->triangles.size();
        trip.reserve(9 * ntri + 4 * links.size());
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const Part& pt = parts[k];
            const TriMesh& m = *pt.mesh;
            for (std::size_t ti = 0; ti < m.triangles.size(); ++ti) {
                const auto& t = m.triangles[ti];
                const double gs = gsheet[k][ti];
                const double x1 = m.x[t[0]], y1 = m.y[t[0]], x2 = m.x[t[1]], y2 = m.y[t[1]], x3 = m.x[t[2]], y3 = m.y[t[2]];
                const double bb[3] = {y2 - y3, y3 - y1, y1 - y2};
                const double cc[3] = {x3 - x2, x1 - x3, x2 - x1};
                const double area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
                if (area2 <= 0) continue;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        trip.push_back({dof[pt.offset + t[i]], dof[pt.offset + t[j]], gs * (bb[i] * bb[j] + cc[i] * cc[j]) / (2.0 * area2)});
            }
        }
        struct Barrel {
            std::size_t a, b;  // dofs
            double g;
            std::size_t drill;
            int la, lb;
            int net;
        };
        std::vector<Barrel> barrels;
        std::size_t unconnected_links = 0;
        const double sigma_p = s.plating_conductivity_s_per_m / thermal;
        const double tp = s.plating_thickness_nm;
        for (std::size_t d = 0; d < b.drills.size(); ++d) {
            std::vector<std::size_t> used;
            for (auto k : drill_links[d]) {
                if (!active_roots.count(uf_islands.find(links[k].island))) continue;
                if (link_node[k] < 0) {
                    ++unconnected_links;
                    continue;
                }
                used.push_back(k);
            }
            const auto& dr = b.drills[d];
            // Plating annulus on the hole wall: the finished hole diameter is the inner one.
            const double slot = dr.slot_end ? std::hypot(double(dr.slot_end->x - dr.at.x), double(dr.slot_end->y - dr.at.y)) : 0.0;
            const double area_nm2 = tp * (kPi * (dr.diameter_nm + tp) + 2.0 * slot);
            for (std::size_t k = 1; k < used.size(); ++k) {
                const auto& la = links[used[k - 1]];
                const auto& lb = links[used[k]];
                const double length_nm = std::abs(zmid[la.layer] - zmid[lb.layer]);
                double sp = sigma_p;
                if (s.temperature_field) {
                    const double tm = 0.5 * (s.temperature_field(la.layer, double(dr.at.x), double(dr.at.y)) +
                                             s.temperature_field(lb.layer, double(dr.at.x), double(dr.at.y)));
                    sp = s.plating_conductivity_s_per_m / (1.0 + s.temperature_coefficient * (tm - 20.0));
                }
                const double g = sp * area_nm2 * 1e-9 / length_nm;  // S
                const std::size_t da = dof[static_cast<std::size_t>(link_node[used[k - 1]])];
                const std::size_t db = dof[static_cast<std::size_t>(link_node[used[k]])];
                barrels.push_back({da, db, g, d, la.layer, lb.layer, islands[la.island].net});
                trip.push_back({da, da, g});
                trip.push_back({db, db, g});
                trip.push_back({da, db, -g});
                trip.push_back({db, da, -g});
            }
        }
        if (unconnected_links)
            res.warnings.push_back(std::to_string(unconnected_links) + " hole walls have no mesh nodes and were left unconnected");
        Csr K = to_csr(ndof, trip);
        trip.clear();
        trip.shrink_to_fit();

        // Boundary conditions. Every port is an area of one potential: set-voltage ports
        // at their voltage, set-current ports a free unknown loaded with their current (one
        // solution; its row couples the port's whole area). Copper connected to ports
        // without a set voltage takes its first port as the 0 V reference.
        std::vector<double> fixed(ndof, NAN);
        std::vector<std::size_t> port_dof(s.ports.size());
        std::map<std::size_t, double> set_current;
        for (std::size_t p = 0; p < s.ports.size(); ++p) {
            port_dof[p] = dof[static_cast<std::size_t>(port_node[p])];
            if (s.ports[p].voltage_v) {
                const double v = *s.ports[p].voltage_v;
                if (!std::isnan(fixed[port_dof[p]]) && fixed[port_dof[p]] != v)
                    throw std::runtime_error("port " + s.ports[p].name + " shares copper with another port at a different voltage");
                fixed[port_dof[p]] = v;
            } else {
                set_current[port_dof[p]] += *s.ports[p].current_a;
            }
        }
        {
            UnionFind uc(ndof);
            for (std::size_t r = 0; r < ndof; ++r)
                for (auto k = K.ptr[r]; k < K.ptr[r + 1]; ++k) uc.unite(r, static_cast<std::size_t>(K.col[k]));
            std::map<std::size_t, std::vector<std::size_t>> ports_of;
            for (std::size_t p = 0; p < s.ports.size(); ++p) ports_of[uc.find(port_dof[p])].push_back(p);
            std::set<std::size_t> referenced;
            for (std::size_t d = 0; d < ndof; ++d)
                if (!std::isnan(fixed[d])) referenced.insert(uc.find(d));
            for (const auto& [root, ps] : ports_of) {
                if (!referenced.count(root)) {
                    double sum = 0, mag = 0;
                    for (auto p : ps) {
                        sum += *s.ports[p].current_a;
                        mag += std::abs(*s.ports[p].current_a);
                    }
                    if (std::abs(sum) > 1e-9 * std::max(mag, 1e-12)) {
                        std::string names;
                        for (auto p : ps) names += (names.empty() ? "" : ", ") + s.ports[p].name;
                        throw std::runtime_error("ports " + names +
                                                 " are connected only to each other and their currents do not sum to zero");
                    }
                    fixed[port_dof[ps[0]]] = 0.0;
                    res.warnings.push_back("ports connected without a set voltage; potentials are relative to port " +
                                           s.ports[ps[0]].name);
                    referenced.insert(root);
                }
            }
            // Copper without a port anywhere in its group carries no current.
            for (std::size_t d = 0; d < ndof; ++d)
                if (!referenced.count(uc.find(d))) {
                    fixed[d] = 0.0;
                    referenced.insert(uc.find(d));
                }
        }
        // Reduced system over the free unknowns; its coupling to the fixed potentials and
        // the set currents form the right-hand side.
        std::vector<long> free_index(ndof, -1);
        std::size_t nfree = 0;
        for (std::size_t d = 0; d < ndof; ++d)
            if (std::isnan(fixed[d])) free_index[d] = static_cast<long>(nfree++);
        std::vector<std::ptrdiff_t> ptr(nfree + 1, 0), col;
        std::vector<double> val;
        std::vector<std::vector<double>> rhs(1, std::vector<double>(nfree, 0.0));
        col.reserve(K.col.size());
        val.reserve(K.val.size());
        for (std::size_t r = 0; r < ndof; ++r) {
            if (free_index[r] < 0) continue;
            const auto fr = static_cast<std::size_t>(free_index[r]);
            for (auto k = K.ptr[r]; k < K.ptr[r + 1]; ++k) {
                const auto c = static_cast<std::size_t>(K.col[k]);
                if (free_index[c] >= 0) {
                    col.push_back(free_index[c]);
                    val.push_back(K.val[k]);
                } else {
                    rhs[0][fr] -= K.val[k] * fixed[c];
                }
            }
            ptr[fr + 1] = static_cast<std::ptrdiff_t>(col.size());
        }
        for (const auto& [d, current] : set_current)
            if (free_index[d] >= 0) rhs[0][static_cast<std::size_t>(free_index[d])] += current;
        res.unknowns = nfree;
        res.triangles = ntri;
        res.seconds_assemble += since(t_assemble);

        // Solve.
        const auto t_solve = Clock::now();
        std::vector<std::vector<double>> sol;
        {
            const SolveStats st = solve_spd(nfree, ptr, col, val, rhs, sol, s.tolerance, s.max_iterations);
            res.iterations = st.iterations;
            res.residual = st.residual;
            if (nfree > 0 && !(st.residual <= s.tolerance * 10))
                res.warnings.push_back("linear solver did not reach the tolerance: residual " + std::to_string(st.residual));
            if (nfree == 0) sol.assign(1, {});
        }
        res.seconds_solve += since(t_solve);

        std::vector<double> V(ndof);
        for (std::size_t d = 0; d < ndof; ++d)
            V[d] = free_index[d] >= 0 ? sol[0][static_cast<std::size_t>(free_index[d])] : fixed[d];

        // Error estimate by gradient recovery. Triangles inside an equipotential area
        // (port, hole wall) are boundary conditions, not field: they stay out.
        std::vector<std::vector<double>> eta2(parts.size());  // per triangle, W
        double eta2_total = 0, energy = 0, area_total = 0;
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const Part& pt = parts[k];
            const TriMesh& m = *pt.mesh;
            std::vector<double> u(m.x.size());
            for (std::size_t n = 0; n < m.x.size(); ++n) u[n] = V[dof[pt.offset + n]];
            std::vector<char> skip(m.triangles.size(), 0);
            for (std::size_t t = 0; t < m.triangles.size(); ++t) {
                const auto& tr = m.triangles[t];
                const auto d0 = dof[pt.offset + tr[0]];
                skip[t] = d0 == dof[pt.offset + tr[1]] && d0 == dof[pt.offset + tr[2]];
            }
            eta2[k] = zz_estimate(m, u, gsheet[k], skip, eta2_total, energy, area_total);
        }
        res.estimated_error = energy > 0 ? std::sqrt(eta2_total / energy) : 0.0;
        res.passes.push_back({ntri, res.estimated_error, energy});
        if (res.estimated_error <= s.target_error || pass + 1 >= s.max_passes) {
            if (res.estimated_error > s.target_error)
                res.warnings.push_back("estimated error " + std::to_string(res.estimated_error) + " above the target after " +
                                       std::to_string(pass + 1) + " passes");
        } else {
            std::vector<RefinePart> rp;
            for (std::size_t k = 0; k < parts.size(); ++k) rp.push_back({parts[k].mesh, &eta2[k], parts[k].island});
            bool budget_reached = false;
            // A margin of 0.9: on large boards it reaches the target in one pass fewer
            // with fewer triangles (docs/validation/dc.md).
            const std::size_t marked = mark_and_seed(rp, eta2_total, energy, s.target_error, ntri, s.max_triangles,
                                                     s.mesh.min_edge_nm, seeds, budget_reached, 0.9);
            if (budget_reached)
                res.warnings.push_back("refinement limited by the triangle budget (" + std::to_string(s.max_triangles) + ")");
            if (marked > 0) continue;
        }

        // Mesh of this solution, for a later solution on the same mesh.
        res.mesh_seeds.assign(islands.size(), {});
        for (const auto& pt : parts)
            for (std::size_t n = 0; n < pt.mesh->x.size(); ++n) res.mesh_seeds[pt.island].push_back({pt.mesh->x[n], pt.mesh->y[n], 0.0});

        // Port voltages and currents. The current into the copper at an unknown is (K V);
        // a set-voltage port supplies what the set-current ports on the same unknown do not.
        std::map<std::size_t, std::vector<std::size_t>> ports_on;
        for (std::size_t p = 0; p < s.ports.size(); ++p) ports_on[port_dof[p]].push_back(p);
        res.ports.resize(s.ports.size());
        for (const auto& [d, ps] : ports_on) {
            double kv = 0, set_current = 0;
            for (auto k = K.ptr[d]; k < K.ptr[d + 1]; ++k) kv += K.val[k] * V[static_cast<std::size_t>(K.col[k])];
            std::vector<std::size_t> voltage_ports;
            for (auto p : ps) {
                if (s.ports[p].current_a) set_current += *s.ports[p].current_a;
                else voltage_ports.push_back(p);
            }
            for (auto p : ps) {
                DcPortResult& pr = res.ports[p];
                pr.name = s.ports[p].name;
                pr.voltage_v = V[d];
                pr.voltage_set = s.ports[p].voltage_v.has_value();
                pr.nodes = port_count[p];
                pr.current_a = pr.voltage_set ? (kv - set_current) / static_cast<double>(voltage_ports.size()) : *s.ports[p].current_a;
            }
            if (ps.size() > 1) {
                std::string names;
                for (auto p : ps) names += (names.empty() ? "" : ", ") + s.ports[p].name;
                res.warnings.push_back("ports " + names + " share copper and form one equipotential area");
            }
            res.power_ports_w += V[d] * kv;
        }

        // Fields per layer.
        std::map<int, std::size_t> layer_slot;
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const Part& pt = parts[k];
            const int layer = islands[pt.island].layer;
            if (!layer_slot.count(layer)) {
                layer_slot[layer] = res.layers.size();
                DcLayerField f;
                f.layer = layer;
                f.thickness_nm = thickness[layer];
                f.conductivity_s_per_m = sigma[layer];
                res.layers.push_back(std::move(f));
            }
            DcLayerField& f = res.layers[layer_slot[layer]];
            const TriMesh& m = *pt.mesh;
            const auto base = static_cast<std::uint32_t>(f.x.size());
            for (std::size_t n = 0; n < m.x.size(); ++n) {
                f.x.push_back(m.x[n]);
                f.y.push_back(m.y[n]);
                f.voltage.push_back(V[dof[pt.offset + n]]);
            }
            for (std::size_t ti = 0; ti < m.triangles.size(); ++ti) {
                const auto& t = m.triangles[ti];
                const double gs = gsheet[k][ti];
                const double x1 = m.x[t[0]], y1 = m.y[t[0]], x2 = m.x[t[1]], y2 = m.y[t[1]], x3 = m.x[t[2]], y3 = m.y[t[2]];
                const double area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
                const double v1 = V[dof[pt.offset + t[0]]], v2 = V[dof[pt.offset + t[1]]], v3 = V[dof[pt.offset + t[2]]];
                // Gradient in V/nm, then V/m.
                const double gx = ((y2 - y3) * v1 + (y3 - y1) * v2 + (y1 - y2) * v3) / area2 * 1e9;
                const double gy = ((x3 - x2) * v1 + (x1 - x3) * v2 + (x2 - x1) * v3) / area2 * 1e9;
                f.triangles.push_back({base + t[0], base + t[1], base + t[2]});
                f.jx.push_back(-gs * gx);
                f.jy.push_back(-gs * gy);
                const double p = gs * (gx * gx + gy * gy) * 0.5 * area2 * 1e-18;
                f.power_w.push_back(p);
                f.net.push_back(islands[pt.island].net);
                res.power_copper_w += p;
            }
        }
        std::sort(res.layers.begin(), res.layers.end(), [](const auto& a, const auto& c) { return a.layer < c.layer; });

        for (const auto& br : barrels) {
            DcBarrelResult r;
            r.drill = br.drill;
            r.net = br.net;
            r.layer_a = br.la;
            r.layer_b = br.lb;
            r.resistance_ohm = 1.0 / br.g;
            r.current_a = br.g * (V[br.a] - V[br.b]);
            res.power_barrels_w += r.current_a * r.current_a * r.resistance_ohm;
            res.barrels.push_back(r);
        }

        // Hot spots: largest volume current densities, at least 1 mm apart.
        struct Cand {
            double density;
            std::size_t layer, tri;
        };
        std::vector<Cand> cands;
        for (std::size_t li = 0; li < res.layers.size(); ++li) {
            const auto& f = res.layers[li];
            for (std::size_t t = 0; t < f.jx.size(); ++t)
                cands.push_back({std::hypot(f.jx[t], f.jy[t]) / (f.thickness_nm * 1e-9) * 1e-6, li, t});
        }
        const std::size_t keep = std::min<std::size_t>(cands.size(), 5000);
        std::partial_sort(cands.begin(), cands.begin() + static_cast<std::ptrdiff_t>(keep), cands.end(),
                          [](const Cand& a, const Cand& c) { return a.density > c.density; });
        for (std::size_t k = 0; k < keep && res.hot_spots.size() < 20; ++k) {
            const auto& f = res.layers[cands[k].layer];
            const auto& t = f.triangles[cands[k].tri];
            const double cx = (f.x[t[0]] + f.x[t[1]] + f.x[t[2]]) / 3, cy = (f.y[t[0]] + f.y[t[1]] + f.y[t[2]]) / 3;
            bool near = false;
            for (const auto& h : res.hot_spots)
                if (h.layer == f.layer && std::hypot(h.x_nm - cx, h.y_nm - cy) < 1e6) near = true;
            if (near) continue;
            res.hot_spots.push_back({f.layer, f.net[cands[k].tri], cx, cy,
                                     std::hypot(f.jx[cands[k].tri], f.jy[cands[k].tri]) * 1e-3, cands[k].density});
        }
        return res;
    }  // adaptive passes
}

}  // namespace kairoq
