// Capacitance of a switch node to a plate by the method of moments (emi_capacitance.hpp).
//
// The plate is a conducting plane under the stack it faces (gap, dielectrics, copper
// layers filled with the dielectric beside them, mask), with air beyond the board. The
// copper of every net within a window around the node is meshed into triangles carrying a
// uniform charge each. With the potential of a unit charge in that medium (layered_green.hpp),
// the Galerkin system P q = v gives the charges for the node at 1 V and all other copper at
// 0 V: the node's capacitance to each net is minus that net's charge, and to the plate the
// rest of the node's charge.

#include "emi_capacitance.hpp"

#include "fem_util.hpp"
#include "layered_green.hpp"
#include "peec.hpp"
#include "kairoq/mesh2d.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>

namespace kairoq::detail {

namespace {

struct Panel {
    std::array<Vec2, 3> v;  // m
    Vec2 c;
    double z = 0, area = 0, size = 0;
    int layer = 0, group = 0;  // group 0: the node; others index the nets
};

double eps_of(const StackupLayer& l, const EmiSettings& s) {
    if (l.epsilon_r) return *l.epsilon_r;
    return l.type == LayerType::SolderMask ? s.mask_epsilon_r : s.dielectric_epsilon_r;
}

}  // namespace

NodeCapacitance node_capacitance(const Board& b, int net, const PlateSpec& plate, const EmiSettings& s) {
    const auto t0 = Clock::now();
    NodeCapacitance res;
    // The stack from the plate's face: the gap, then the board's layers away from the plate.
    Stratified st;
    st.ground = true;
    if (plate.gap_nm > 0) {
        st.thickness.push_back(plate.gap_nm * 1e-9);
        st.epsilon_r.push_back(plate.gap_epsilon_r);
    }
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < b.stackup.size(); ++i) order.push_back(i);
    if (plate.bottom) std::reverse(order.begin(), order.end());
    // A copper layer is a sheet at its middle: half its thickness joins the slab on each
    // side, so that the sheet lies on the interface between them and no slab is thinner than
    // the dielectrics.
    std::map<int, double> zcu;  // copper layer -> height of its middle above the plate (m)
    double z = st.thickness.empty() ? 0.0 : st.thickness.front();
    double carry = 0;
    for (std::size_t q = 0; q < order.size(); ++q) {
        const auto& l = b.stackup[order[q]];
        if (!(l.thickness_nm > 0)) continue;
        if (l.type == LayerType::Silkscreen || l.type == LayerType::SolderPaste) continue;
        const double t = l.thickness_nm * 1e-9;
        if (l.type == LayerType::Copper) {
            if (!st.thickness.empty()) {
                st.thickness.back() += 0.5 * t;
                carry += 0.5 * t;
            } else {
                carry += t;  // on the plate's face: the whole layer goes to the slab above
            }
            zcu[l.copper_index] = z + (st.thickness.empty() ? 0.0 : 0.5 * t);
        } else {
            st.thickness.push_back(t + carry);
            st.epsilon_r.push_back(eps_of(l, s));
            carry = 0;
        }
        z += t;
    }

    // Window: the node's copper over the plate, with a margin.
    const PolygonSet plate_area = plate.area.empty() ? b.outline : intersection(plate.area, b.outline);
    PolygonSet node_all;
    double zmax = 0;
    for (const auto& cl : b.copper)
        for (const auto& nc : cl.nets)
            if (nc.net == net) {
                const PolygonSet r = intersection(nc.polygons, plate_area);
                if (r.empty()) continue;
                node_all.insert(node_all.end(), r.begin(), r.end());
                zmax = std::max(zmax, zcu.count(cl.index) ? zcu[cl.index] : 0.0);
            }
    if (node_all.empty()) return res;
    const double margin_m = s.coupling.window_margin_nm > 0 ? s.coupling.window_margin_nm * 1e-9 : std::clamp(3.0 * zmax, 2e-3, 30e-3);
    res.window_margin_m = margin_m;
    const Rect64 nb = bounds(node_all);
    const double mg = margin_m * 1e9;
    const Rect64 wr{static_cast<int64_t>(nb.left - mg), static_cast<int64_t>(nb.top - mg), static_cast<int64_t>(nb.right + mg),
                    static_cast<int64_t>(nb.bottom + mg)};
    const PolygonSet window{Polygon{Path64{{wr.left, wr.top}, {wr.right, wr.top}, {wr.right, wr.bottom}, {wr.left, wr.bottom}}, {}}};
    if (!plate.area.empty() && area_nm2(difference(intersection(window, b.outline), plate_area)) > 1e-3 * area_nm2(window))
        res.warnings.push_back("the plate does not cover the window around " + b.nets[net] + "; it is taken as an infinite plane");

    // Conductors: per layer and net, the copper within the window.
    struct Piece {
        int layer, group;
        PolygonSet shape;
    };
    std::vector<Piece> pieces;
    std::vector<int> group_net{net};
    std::map<int, int> group_of{{net, 0}};
    for (const auto& cl : b.copper) {
        if (!zcu.count(cl.index)) continue;
        for (const auto& nc : cl.nets) {
            PolygonSet r = intersection(nc.polygons, window);
            if (nc.net == net) r = intersection(r, plate_area);
            if (r.empty()) continue;
            if (!group_of.count(nc.net)) {
                group_of[nc.net] = static_cast<int>(group_net.size());
                group_net.push_back(nc.net);
            }
            pieces.push_back({cl.index, group_of[nc.net], std::move(r)});
        }
    }
    // Element size at the node's edges: a quarter of its distance to the plate or to the
    // nearest other copper layer (twice that on other copper), grown to the budget.
    double dmin = INFINITY;
    for (const auto& p : pieces)
        if (p.group == 0) {
            dmin = std::min(dmin, zcu[p.layer]);
            for (const auto& q : pieces)
                if (q.layer != p.layer) dmin = std::min(dmin, std::abs(zcu[q.layer] - zcu[p.layer]));
        }
    const double h0 = s.coupling.min_size_nm > 0 ? s.coupling.min_size_nm * 1e-9 : std::clamp(dmin / 4.0, 20e-6, 0.5e-3);
    // Away from the node the size grows by 0.35 per unit distance from its bounding box, up
    // to a third of the margin; there other copper is simplified to the local size: holes
    // smaller than half of it (antipads in planes) are closed, outlines simplified to 0.15 of
    // it.
    // Large elements of copper stacked closely on two layers make the system ill-conditioned:
    // the largest is at most eight times the smallest spacing of copper layers.
    double layer_gap = INFINITY;
    for (const auto& p : pieces)
        for (const auto& q : pieces)
            if (q.layer != p.layer) layer_gap = std::min(layer_gap, std::abs(zcu[q.layer] - zcu[p.layer]));
    const double hmax = std::min(margin_m / 3, std::isfinite(layer_gap) ? 8 * layer_gap : INFINITY);
    auto distance = [&](double x, double y) {  // nm, to the node's bounding box
        const double dx = std::max({0.0, double(nb.left) - x, x - double(nb.right)});
        const double dy = std::max({0.0, double(nb.top) - y, y - double(nb.bottom)});
        return std::hypot(dx, dy);
    };
    // Pieces in order: the node's, then the others by decreasing area. A simplified piece
    // gives up what the pieces before it on its layer already cover, so that pads within a
    // plane's closed antipads join the plane (all copper but the node's is at 0 V).
    std::vector<std::size_t> rank(pieces.size());
    std::vector<double> piece_area(pieces.size());
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        rank[i] = i;
        piece_area[i] = area_nm2(pieces[i].shape);
    }
    std::stable_sort(rank.begin(), rank.end(), [&](std::size_t x, std::size_t y) {
        if ((pieces[x].group == 0) != (pieces[y].group == 0)) return pieces[x].group == 0;
        return piece_area[x] > piece_area[y];
    });
    auto mesh = [&](double scale) {
        std::vector<Panel> panels;
        std::map<int, Paths64> placed;  // per layer
        for (const std::size_t pi : rank) {
            const Piece& p = pieces[pi];
            const double h = h0 * scale * (p.group == 0 ? 1.0 : 2.0);
            auto local = [&](double x, double y) { return std::min(std::max(hmax, h), h + 0.35 * distance(x, y) * 1e-9); };
            MeshSettings ms;
            ms.min_size_nm = h * 1e9;
            ms.boundary_tolerance_nm = std::min(2000.0, 0.1 * h0 * 1e9);
            ms.grading = 0.3;
            ms.max_edge_nm = std::max(hmax, h) * 1e9;
            ms.boundary_size = [&](double x, double y) { return local(x, y) * 1e9; };
            PolygonSet shape = p.shape;
            if (p.group != 0) {
                PolygonSet simplified;
                for (const auto& isl : p.shape) {
                    auto size_of = [&](const Path64& r) {
                        const Rect64 bb = bounds(Polygon{r, {}});
                        const double cx = 0.5 * (bb.left + bb.right), cy = 0.5 * (bb.top + bb.bottom);
                        const double extent = 0.5 * std::hypot(double(bb.right - bb.left), double(bb.bottom - bb.top));
                        const double d = std::max(0.0, distance(cx, cy) - extent) * 1e-9;
                        return std::make_pair(std::min(std::max(hmax, h), h + 0.35 * d) * 1e9, extent);
                    };
                    Paths64 paths;
                    const auto [ho, eo] = size_of(isl.outer);
                    paths.push_back(Clipper2Lib::SimplifyPath(isl.outer, std::max(1000.0, 0.15 * ho), true));
                    for (const auto& hole : isl.holes) {
                        const auto [hh, eh] = size_of(hole);
                        if (2 * eh < 0.5 * hh) continue;
                        paths.push_back(Clipper2Lib::SimplifyPath(hole, std::max(1000.0, 0.15 * hh), true));
                    }
                    for (auto& q : union_polygons(paths)) simplified.push_back(std::move(q));
                }
                auto& before = placed[p.layer];
                shape = before.empty() ? std::move(simplified) : difference(simplified, union_polygons(before));
            }
            {
                const auto sp = to_paths(shape);
                auto& before = placed[p.layer];
                before.insert(before.end(), sp.begin(), sp.end());
            }
            // Other copper carries the node's image under and around it: vertices at the local
            // size over the node's bounding box and in rings around it.
            std::vector<std::array<double, 3>> seeds;
            if (p.group != 0) {
                const double hs = 2 * h0 * scale * 1e9;  // nm
                for (double x = nb.left; x <= nb.right; x += hs)
                    for (double y = nb.top; y <= nb.bottom; y += hs) seeds.push_back({x, y, hs});
                for (double d = hs, step = hs; d < margin_m * 1e9; d += step) {
                    step = std::min(std::max(hmax, h), h + 0.35 * d * 1e-9) * 1e9;
                    const double x0 = nb.left - d, x1 = nb.right + d, y0 = nb.top - d, y1 = nb.bottom + d;
                    for (double x = x0; x <= x1; x += step) {
                        seeds.push_back({x, y0, step});
                        seeds.push_back({x, y1, step});
                    }
                    for (double y = y0 + step; y < y1; y += step) {
                        seeds.push_back({x0, y, step});
                        seeds.push_back({x1, y, step});
                    }
                }
            }
            for (const auto& poly : shape) {
                const TriMesh m = mesh_region(poly, {}, ms, seeds);
                for (const auto& t : m.triangles) {
                    Panel pn;
                    for (int k = 0; k < 3; ++k) pn.v[k] = {m.x[t[k]] * 1e-9, m.y[t[k]] * 1e-9};
                    pn.c = {(pn.v[0].x + pn.v[1].x + pn.v[2].x) / 3, (pn.v[0].y + pn.v[1].y + pn.v[2].y) / 3};
                    pn.area = 0.5 * std::abs((pn.v[1].x - pn.v[0].x) * (pn.v[2].y - pn.v[0].y) - (pn.v[2].x - pn.v[0].x) * (pn.v[1].y - pn.v[0].y));
                    for (int k = 0; k < 3; ++k)
                        pn.size = std::max(pn.size, std::hypot(pn.v[(k + 1) % 3].x - pn.v[k].x, pn.v[(k + 1) % 3].y - pn.v[k].y));
                    pn.z = zcu[p.layer];
                    pn.layer = p.layer;
                    pn.group = p.group;
                    if (pn.area > 0) panels.push_back(pn);
                }
            }
        }
        return panels;
    };
    double scale = 1.0;
    std::vector<Panel> panels = mesh(scale);
    for (int attempt = 0; attempt < 8 && panels.size() > s.coupling.max_panels; ++attempt) {
        scale *= std::clamp(std::sqrt(double(panels.size()) / (0.8 * s.coupling.max_panels)), 1.1, 3.0);
        panels = mesh(scale);
    }
    if (panels.size() > 2 * s.coupling.max_panels) {
        res.warnings.push_back("the copper around " + b.nets[net] + " needs " + std::to_string(panels.size()) +
                               " panels at the coarsest size, more than twice the budget; its capacitance with fringing is not computed");
        res.panels = panels.size();
        return res;
    }
    if (panels.size() > s.coupling.max_panels)
        res.warnings.push_back("the copper around " + b.nets[net] + " needs " + std::to_string(panels.size()) + " panels, above the budget");
    res.panels = panels.size();

    // The medium's potential for each pair of copper layers.
    const double rho_max = std::hypot(double(wr.right - wr.left), double(wr.bottom - wr.top)) * 1e-9;
    std::map<std::pair<int, int>, std::unique_ptr<LayeredGreen>> green;
    for (const auto& [la, za] : zcu)
        for (const auto& [lb, zb] : zcu) {
            if (lb < la) continue;
            bool used_a = false, used_b = false;
            for (const auto& p : pieces) {
                used_a |= p.layer == la;
                used_b |= p.layer == lb;
            }
            if (used_a && used_b) green[{la, lb}] = std::make_unique<LayeredGreen>(st, za, zb, rho_max);
        }
    auto G = [&](int la, int lb) -> const LayeredGreen& { return *green.at({std::min(la, lb), std::max(la, lb)}); };

    // Galerkin system: the mean potential over panel i of a unit charge spread over panel j.
    // Returns the charge of each group with the node at 1 V.
    const TriangleRule& rule = triangle_rule_7();
    const TriangleRule& fine_rule = triangle_rule_28();
    auto charges = [&](const std::vector<Panel>& panels) {
        const std::size_t n = panels.size();
        std::vector<double> P(n * n, 0.0);
        parallel_for(n, [&](std::size_t i) {
            const Panel& a = panels[i];
            std::vector<Vec2> pa(rule.weight.size());
            for (std::size_t q = 0; q < pa.size(); ++q) {
                const auto& w = rule.bary[q];
                pa[q] = {w[0] * a.v[0].x + w[1] * a.v[1].x + w[2] * a.v[2].x, w[0] * a.v[0].y + w[1] * a.v[1].y + w[2] * a.v[2].y};
            }
            for (std::size_t j = i; j < n; ++j) {
                const Panel& c = panels[j];
                const LayeredGreen& g = G(a.layer, c.layer);
                const double d = std::hypot(a.c.x - c.c.x, a.c.y - c.c.y);
                double v = 0;
                if (d > 2.5 * (a.size + c.size)) {
                    v = g(d);
                } else {
                    // The direct term in closed form over panel j, the remainder by quadrature;
                    // over panel i a finer rule where the panels are large against their
                    // vertical distance.
                    const double dz = std::abs(a.z - c.z);
                    const bool fine = dz > 0 && std::max(a.size, c.size) > 2 * dz;
                    const TriangleRule& outer = fine ? fine_rule : rule;
                    std::vector<Vec2> po(outer.weight.size());
                    for (std::size_t q = 0; q < po.size(); ++q) {
                        const auto& w = outer.bary[q];
                        po[q] = {w[0] * a.v[0].x + w[1] * a.v[1].x + w[2] * a.v[2].x, w[0] * a.v[0].y + w[1] * a.v[1].y + w[2] * a.v[2].y};
                    }
                    for (std::size_t q = 0; q < po.size(); ++q) {
                        double rem = 0;
                        for (std::size_t r = 0; r < pa.size(); ++r) {
                            const auto& w = rule.bary[r];
                            const double x = w[0] * c.v[0].x + w[1] * c.v[1].x + w[2] * c.v[2].x;
                            const double y = w[0] * c.v[0].y + w[1] * c.v[1].y + w[2] * c.v[2].y;
                            rem += rule.weight[r] * g.remainder(std::hypot(po[q].x - x, po[q].y - y));
                        }
                        double direct = g.direct_coefficient() * triangle_potential(c.v, po[q], dz).scalar;
                        for (const auto& im : g.images()) direct += im.coefficient * triangle_potential(c.v, po[q], dz + im.depth).scalar;
                        v += outer.weight[q] * (direct / c.area + rem);
                    }
                }
                P[i * n + j] = v;
                P[j * n + i] = v;
            }
        });
        if (!cholesky(P, n)) throw std::runtime_error("switch node " + b.nets[net] + ": the potential matrix is not positive definite");
        std::vector<double> q(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) q[i] = panels[i].group == 0 ? 1.0 : 0.0;
        forward_solve(P, n, q, 1);
        backward_solve(P, n, q, 1);
        std::vector<double> out(group_net.size(), 0.0);
        for (std::size_t i = 0; i < n; ++i) out[panels[i].group] += q[i];
        return out;
    };
    // The Galerkin charges converge from below, in proportion to the element size at the
    // edges: a second solution with elements twice as large extrapolates them.
    const std::vector<double> fine = charges(panels);
    std::vector<double> best = fine;
    const std::vector<Panel> coarse_panels = mesh(2 * scale);
    if (coarse_panels.size() < panels.size()) {
        const std::vector<double> coarse = charges(coarse_panels);
        for (double v : coarse) res.coarse_to_plate_f += v;
        for (std::size_t g = 0; g < fine.size(); ++g) best[g] = 2 * fine[g] - coarse[g];
    }
    res.total_f = best[0];
    res.to_plate_f = best[0];
    for (std::size_t g = 1; g < group_net.size(); ++g) {
        res.to_plate_f += best[g];
        res.to_nets.push_back({group_net[g], -best[g]});
    }
    res.unextrapolated_to_plate_f = 0;
    for (double v : fine) res.unextrapolated_to_plate_f += v;
    std::sort(res.to_nets.begin(), res.to_nets.end(), [](const auto& x, const auto& y) { return x.second > y.second; });
    res.solved = true;
    res.seconds = since(t0);
    return res;
}

}  // namespace kairoq::detail
