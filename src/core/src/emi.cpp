// Return-path checks, switch-node coupling, and the EMI analysis as a whole.
#include "kairoq/emi.hpp"

#include "emi_capacitance.hpp"
#include "fem_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

using namespace detail;

constexpr double kEps0 = 8.8541878128e-12;  // F/m

PolygonSet grow(const PolygonSet& p, double d_nm) {
    if (p.empty()) return {};
    return union_polygons(
        Clipper2Lib::InflatePaths(to_paths(p), d_nm, Clipper2Lib::JoinType::Round, Clipper2Lib::EndType::Polygon, 2.0, 50.0));
}

double perimeter_nm(const Polygon& p) {
    auto ring = [](const Path64& r) {
        double l = 0;
        for (std::size_t k = 0; k < r.size(); ++k) {
            const auto& a = r[k];
            const auto& b = r[(k + 1) % r.size()];
            l += std::hypot(double(b.x - a.x), double(b.y - a.y));
        }
        return l;
    };
    double l = ring(p.outer);
    for (const auto& h : p.holes) l += ring(h);
    return l;
}

const PolygonSet* net_copper(const Board& b, int layer, int net) {
    for (const auto& cl : b.copper)
        if (cl.index == layer)
            for (const auto& nc : cl.nets)
                if (nc.net == net) return &nc.polygons;
    return nullptr;
}

}  // namespace

std::vector<ReturnPathFinding> check_return_paths(const Board& b, const ReturnPathSettings& s, std::vector<std::string>& warnings) {
    const double board_area = area_nm2(b.outline);
    const int n = static_cast<int>(b.copper.size());
    // Reference nets: given, or those that fill a share of some layer.
    std::set<int> refs;
    for (const auto& name : s.reference_nets) {
        const int id = b.find_net(name);
        if (id < 0) throw std::runtime_error("return path: no net '" + name + "'");
        refs.insert(id);
    }
    if (refs.empty())
        for (const auto& cl : b.copper)
            for (const auto& nc : cl.nets)
                if (nc.net != 0 && area_nm2(nc.polygons) >= s.plane_fraction * board_area) refs.insert(nc.net);
    if (refs.empty()) {
        warnings.push_back("return path: no reference plane found");
        return {};
    }
    // Plane layers: those on which the reference nets fill the share.
    std::vector<PolygonSet> ref_union(n + 1);
    std::vector<bool> plane(n + 1, false);
    for (const auto& cl : b.copper) {
        Paths64 all;
        for (const auto& nc : cl.nets)
            if (refs.count(nc.net)) {
                const auto p = to_paths(nc.polygons);
                all.insert(all.end(), p.begin(), p.end());
            }
        ref_union[cl.index] = union_polygons(all);
        plane[cl.index] = area_nm2(ref_union[cl.index]) >= s.plane_fraction * board_area;
    }
    std::set<int> checked;
    for (const auto& name : s.nets) {
        const int id = b.find_net(name);
        if (id < 0) throw std::runtime_error("return path: no net '" + name + "'");
        checked.insert(id);
    }
    if (checked.empty())
        for (std::size_t id = 1; id < b.nets.size(); ++id)
            if (!refs.count(static_cast<int>(id))) checked.insert(static_cast<int>(id));
    const double min_area = s.min_area_mm2 * 1e12;
    std::vector<ReturnPathFinding> out;
    for (const auto& cl : b.copper) {
        const int L = cl.index;
        // Nearest plane layers above and below.
        std::vector<int> rlayers;
        for (int k = L - 1; k >= 1; --k)
            if (plane[k]) {
                rlayers.push_back(k);
                break;
            }
        for (int k = L + 1; k <= n; ++k)
            if (plane[k]) {
                rlayers.push_back(k);
                break;
            }
        if (rlayers.empty()) continue;
        for (const auto& nc : cl.nets) {
            if (!checked.count(nc.net)) continue;
            // The trace without its pads and holes (where it changes layer or meets a part).
            Paths64 excl;
            for (const auto& pad : b.pads)
                if (pad.layer == L && pad.net == nc.net) {
                    const auto p = to_paths(pad.shape);
                    excl.insert(excl.end(), p.begin(), p.end());
                }
            PolygonSet exclusion = grow(union_polygons(excl), s.pad_clearance_nm);
            Paths64 holes;
            for (const auto& d : b.drills)
                if (d.plated && L >= d.from_layer && L <= d.to_layer && contains(nc.polygons, d.at))
                    holes.push_back(circle_path(d.at, 0.5 * d.diameter_nm + s.pad_clearance_nm, 2000.0));
            if (!holes.empty()) exclusion = unite(exclusion, union_polygons(holes));
            const PolygonSet trace = intersection(difference(nc.polygons, exclusion), b.outline);
            for (int R : rlayers) {
                for (const auto& piece : trace) {
                    const PolygonSet P{piece};
                    const double pa = area_nm2(piece);
                    if (pa < min_area) continue;
                    const PolygonSet gaps = difference(P, ref_union[R]);
                    if (gaps.empty()) continue;
                    // Parts over each reference net.
                    std::vector<std::pair<int, PolygonSet>> covered;
                    for (int r : refs)
                        if (const PolygonSet* rc = net_copper(b, R, r))
                            for (const auto& c : intersection(P, *rc)) covered.push_back({r, PolygonSet{c}});
                    const double width = 2.0 * pa / std::max(1.0, perimeter_nm(piece));
                    for (const auto& g : gaps) {
                        const double ga = area_nm2(g);
                        if (ga < min_area) continue;
                        const PolygonSet near = grow(PolygonSet{g}, 2000.0);
                        std::vector<int> touching;
                        for (const auto& [net, c] : covered)
                            if (!intersection(near, c).empty()) touching.push_back(net);
                        ReturnPathFinding f;
                        f.net = b.nets[nc.net];
                        f.layer = L;
                        f.reference_layer = R;
                        std::set<int> tn(touching.begin(), touching.end());
                        for (int t : tn) f.reference_nets.push_back(b.nets[t]);
                        f.kind = touching.size() >= 2 ? (tn.size() >= 2 ? "split" : "crossing") : "partial";
                        // Centroid by the polygon's vertices weighted by area of the outer ring.
                        double cx = 0, cy = 0, a2 = 0;
                        const auto& r = g.outer;
                        for (std::size_t k = 0; k < r.size(); ++k) {
                            const auto& p0 = r[k];
                            const auto& p1 = r[(k + 1) % r.size()];
                            const double cr = double(p0.x) * double(p1.y) - double(p1.x) * double(p0.y);
                            a2 += cr;
                            cx += (double(p0.x) + double(p1.x)) * cr;
                            cy += (double(p0.y) + double(p1.y)) * cr;
                        }
                        f.x_nm = a2 != 0 ? cx / (3 * a2) : double(r[0].x);
                        f.y_nm = a2 != 0 ? cy / (3 * a2) : double(r[0].y);
                        f.area_mm2 = ga * 1e-12;
                        f.span_mm = ga / std::max(1.0, width) * 1e-6;
                        f.piece = PolygonSet{g};
                        out.push_back(std::move(f));
                    }
                }
            }
        }
    }
    // Crossings and splits first, larger spans first.
    auto rank = [](const ReturnPathFinding& f) { return f.kind == "partial" ? 1 : 0; };
    std::stable_sort(out.begin(), out.end(), [&](const auto& a, const auto& c) {
        if (rank(a) != rank(c)) return rank(a) < rank(c);
        return a.span_mm > c.span_mm;
    });
    return out;
}

CouplingResult switch_node_coupling(const Board& b, const SwitchNodeSpec& node, const PlateSpec& plate, const EmiSettings& s) {
    CouplingResult res;
    res.net = node.net;
    res.plate = plate.name;
    const int net = b.find_net(node.net);
    if (net < 0) throw std::runtime_error("switch node: no net '" + node.net + "'");
    const int n = static_cast<int>(b.copper.size());
    const PolygonSet plate_area = plate.area.empty() ? b.outline : intersection(plate.area, b.outline);
    // Stackup positions of the copper layers.
    std::vector<std::size_t> pos(n + 1, 0);
    for (std::size_t i = 0; i < b.stackup.size(); ++i)
        if (b.stackup[i].type == LayerType::Copper) pos[b.stackup[i].copper_index] = i;
    auto eps_of = [&](const StackupLayer& l) {
        if (l.epsilon_r) return *l.epsilon_r;
        return l.type == LayerType::SolderMask ? s.mask_epsilon_r : s.dielectric_epsilon_r;
    };
    // Electrical thickness (t / eps_r, m) from the bottom of a copper layer to the plate's
    // face, through the dielectrics, the copper layers in between (filled with resin where
    // they have no copper), the mask and the gap.
    auto electrical = [&](int L) {
        double d = 0;
        const std::size_t i0 = pos[L];
        if (plate.bottom) {
            for (std::size_t i = i0 + 1; i < b.stackup.size(); ++i) {
                const auto& l = b.stackup[i];
                double eps = eps_of(l);
                if (l.type == LayerType::Copper) eps = i + 1 < b.stackup.size() ? eps_of(b.stackup[i + 1]) : s.dielectric_epsilon_r;
                d += l.thickness_nm * 1e-9 / eps;
            }
        } else {
            for (std::size_t i = i0; i-- > 0;) {
                const auto& l = b.stackup[i];
                double eps = eps_of(l);
                if (l.type == LayerType::Copper) eps = i > 0 ? eps_of(b.stackup[i - 1]) : s.dielectric_epsilon_r;
                d += l.thickness_nm * 1e-9 / eps;
            }
        }
        return d + plate.gap_nm * 1e-9 / plate.gap_epsilon_r;
    };
    for (int L = 1; L <= n; ++L) {
        const PolygonSet* own = net_copper(b, L, net);
        if (!own) continue;
        const PolygonSet region = intersection(*own, plate_area);
        if (region.empty()) continue;
        // Copper of any net on the layers between shields the plate.
        Paths64 shield;
        const int a = plate.bottom ? L + 1 : 1, z = plate.bottom ? n : L - 1;
        for (int k = a; k <= z; ++k)
            for (const auto& cl : b.copper)
                if (cl.index == k)
                    for (const auto& nc : cl.nets) {
                        const auto p = to_paths(nc.polygons);
                        shield.insert(shield.end(), p.begin(), p.end());
                    }
        const PolygonSet open = shield.empty() ? region : difference(region, union_polygons(shield));
        const double ao = area_nm2(open) * 1e-18, ar = area_nm2(region) * 1e-18;
        res.unshielded_area_mm2 += ao * 1e6;
        res.shielded_area_mm2 += (ar - ao) * 1e6;
        const double d = electrical(L);
        const double c = d > 0 ? kEps0 * ao / d : 0.0;
        if (ao > 0) res.layers.push_back({L, c * 1e12});
        res.parallel_plate_pf += c * 1e12;
    }
    // With fringing: the method of moments over the copper around the node.
    const NodeCapacitance mom = node_capacitance(b, net, plate, s);
    // Without the solution (too many panels), the parallel-plate value.
    res.capacitance_pf = mom.solved ? mom.to_plate_f * 1e12 : res.parallel_plate_pf;
    res.total_pf = mom.total_f * 1e12;
    res.mesh_pf = mom.unextrapolated_to_plate_f * 1e12;
    res.coarse_pf = mom.coarse_to_plate_f * 1e12;
    for (const auto& [k, c] : mom.to_nets) res.nets_pf.push_back({b.nets[k], c * 1e12});
    res.panels = mom.panels;
    res.window_margin_mm = mom.window_margin_m * 1e3;
    res.warnings = mom.warnings;
    res.common_mode_current_ma = res.capacitance_pf * node.dv_dt_v_per_ns;  // pF * V/ns = mA
    // Overlap with other nets on the neighbouring copper layers.
    std::map<std::pair<int, int>, double> nb;  // (net, layer) -> pF
    for (int L = 1; L <= n; ++L) {
        const PolygonSet* own = net_copper(b, L, net);
        if (!own) continue;
        for (int K : {L - 1, L + 1}) {
            if (K < 1 || K > n) continue;
            double d = 0;
            const std::size_t i0 = std::min(pos[L], pos[K]), i1 = std::max(pos[L], pos[K]);
            for (std::size_t i = i0 + 1; i < i1; ++i) d += b.stackup[i].thickness_nm * 1e-9 / eps_of(b.stackup[i]);
            if (d <= 0) continue;
            for (const auto& cl : b.copper)
                if (cl.index == K)
                    for (const auto& nc : cl.nets) {
                        if (nc.net == net) continue;
                        const double ov = area_nm2(intersection(*own, nc.polygons)) * 1e-18;
                        if (ov > 0) nb[{nc.net, K}] += kEps0 * ov / d * 1e12;
                    }
        }
    }
    for (const auto& [key, c] : nb) res.neighbours.push_back({b.nets[key.first], key.second, c});
    std::sort(res.neighbours.begin(), res.neighbours.end(), [](const auto& a, const auto& c) { return std::get<2>(a) > std::get<2>(c); });
    if (res.neighbours.size() > 5) res.neighbours.resize(5);
    return res;
}

EmiResult analyse_emi(const Board& b, const EmiSettings& s) {
    EmiResult r;
    for (const auto& loop : s.loops) r.loops.push_back(solve_loop(b, loop, s.loop));
    if (s.return_path) r.return_path = check_return_paths(b, *s.return_path, r.warnings);
    for (const auto& sw : s.switch_nodes) {
        const PlateSpec* plate = nullptr;
        for (const auto& p : s.plates)
            if (p.name == sw.plate) plate = &p;
        if (!plate) throw std::runtime_error("switch node " + sw.net + ": no plate '" + sw.plate + "'");
        r.coupling.push_back(switch_node_coupling(b, sw, *plate, s));
        r.warnings.insert(r.warnings.end(), r.coupling.back().warnings.begin(), r.coupling.back().warnings.end());
    }
    return r;
}

std::string emi_summary(const Board& b, const EmiResult& r) {
    std::ostringstream s;
    char line[320];
    auto dbuv = [](double e) { return e > 0 ? 20.0 * std::log10(e / 1e-6) : -INFINITY; };
    for (const auto& l : r.loops) {
        const double m = std::sqrt(l.moment_m2[0] * l.moment_m2[0] + l.moment_m2[1] * l.moment_m2[1] + l.moment_m2[2] * l.moment_m2[2]);
        std::snprintf(line, sizeof line, "loop %-12s %8.3f nH", l.name.c_str(), l.inductance_nh);
        s << line;
        if (l.inductance_coarse_nh > 0) {
            std::snprintf(line, sizeof line, " (%.3f nH with about half the triangles; it converges from above, and the error is at least the difference)",
                          l.inductance_coarse_nh);
            s << line;
        }
        if (l.element_inductance_nh > 0) {
            std::snprintf(line, sizeof line, ", of which %.3f nH inside the parts", l.element_inductance_nh);
            s << line;
        }
        std::snprintf(line, sizeof line,
                      "\n  magnetic moment %.4g mm^2 per ampere (window margin %.0f mm, last change %.1f %%); %zu triangles (elements of %.2f mm "
                      "at the loop before refinement), %zu filaments, %.1f s\n",
                      m * 1e6, l.window_margin_mm, 100 * l.moment_change, l.triangles, l.element_size_mm, l.filaments, l.seconds);
        s << line;
        for (std::size_t j = 0; j < l.stitch_current.size(); ++j) {
            std::snprintf(line, sizeof line, "  stitch %zu carries %.3g of the loop current in the high-frequency limit\n", j + 1, l.stitch_current[j]);
            s << line;
        }
        for (const auto& z : l.impedance) {
            std::snprintf(line, sizeof line, "  at %-10.4g Hz  R %.4g mohm, L %.4g nH", z.frequency_hz, z.resistance_ohm * 1e3, z.inductance_nh);
            s << line;
            for (std::size_t j = 0; j < z.stitch_current.size(); ++j) {
                std::snprintf(line, sizeof line, "; stitch %zu %.3g at %.0f deg", j + 1, std::abs(z.stitch_current[j]),
                              std::arg(z.stitch_current[j]) * 180.0 / 3.14159265358979323846);
                s << line;
            }
            s << "\n";
        }
        if (!l.harmonics.empty()) {
            const char* site = l.field_site == "ground"       ? "over a reflecting ground, antenna height scanned"
                               : l.field_site == "small_loop" ? "small-loop far field, free space"
                                                              : "free space, worst direction";
            std::snprintf(line, sizeof line, "  radiated   largest %.1f dBuV/m at %.4g MHz (at %g m, %s)\n", dbuv(l.max_field_v_per_m),
                          l.max_field_frequency_hz * 1e-6, l.field_distance_m, site);
            s << line;
            for (const auto& band : std::vector<std::pair<double, double>>{{30e6, 230e6}, {230e6, 1e9}}) {
                double e = 0, f = 0;
                for (const auto& h : l.harmonics)
                    if (h.frequency_hz >= band.first && h.frequency_hz <= band.second && h.field_v_per_m > e) {
                        e = h.field_v_per_m;
                        f = h.frequency_hz;
                    }
                if (e > 0) {
                    std::snprintf(line, sizeof line, "             %4.0f-%4.0f MHz: largest %.1f dBuV/m at %.4g MHz\n", band.first * 1e-6,
                                  band.second * 1e-6, dbuv(e), f * 1e-6);
                    s << line;
                }
            }
        }
        for (const auto& w : l.warnings) s << "  warning: " << w << "\n";
    }
    if (!r.return_path.empty()) {
        std::size_t crossings = 0;
        for (const auto& f : r.return_path) crossings += f.kind != "partial";
        std::snprintf(line, sizeof line, "return path  %zu crossings or splits, %zu partial gaps\n", crossings, r.return_path.size() - crossings);
        s << line;
        std::size_t shown = 0;
        for (const auto& f : r.return_path) {
            if (shown++ >= 12) break;
            std::string refs;
            for (const auto& n : f.reference_nets) refs += (refs.empty() ? "" : "/") + n;
            std::snprintf(line, sizeof line, "  %-8s %-14s L%d over L%d %-10s at (%.3f, %.3f) mm, span %.3f mm, %.4f mm^2\n", f.kind.c_str(),
                          f.net.c_str(), f.layer, f.reference_layer, refs.c_str(), to_mm(f.x_nm), to_mm(f.y_nm), f.span_mm, f.area_mm2);
            s << line;
        }
    }
    for (const auto& c : r.coupling) {
        std::snprintf(line, sizeof line, "switch node %-10s to %-10s %8.3f pF (%.2f mm^2 facing it, %.2f mm^2 shielded)", c.net.c_str(),
                      c.plate.c_str(), c.capacitance_pf, c.unshielded_area_mm2, c.shielded_area_mm2);
        s << line;
        if (c.common_mode_current_ma > 0) {
            std::snprintf(line, sizeof line, "; %.3g mA at the given dV/dt", c.common_mode_current_ma);
            s << line;
        }
        s << "\n";
        if (c.panels > 0) {
            std::snprintf(line, sizeof line,
                          "  with fringing, %zu panels within %.1f mm: %.4g pF on the mesh, %.4g pF with elements twice as large; "
                          "parallel plates %.4g pF\n",
                          c.panels, c.window_margin_mm, c.mesh_pf, c.coarse_pf, c.parallel_plate_pf);
            s << line;
            for (std::size_t k = 0; k < c.nets_pf.size() && k < 5; ++k) {
                std::snprintf(line, sizeof line, "  to %-14s %.4g pF\n", c.nets_pf[k].first.c_str(), c.nets_pf[k].second);
                s << line;
            }
        }
        for (const auto& [net, layer, pf] : c.neighbours) {
            std::snprintf(line, sizeof line, "  overlap with %-14s on L%d: %.3f pF (parallel plates)\n", net.c_str(), layer, pf);
            s << line;
        }
    }
    for (const auto& w : r.warnings) s << "warning: " << w << "\n";
    (void)b;
    return s.str();
}

}  // namespace kairoq
