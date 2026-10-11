// Export of a loop's model to openEMS: the geometry the loop solver uses (the window, the
// copper islands in it, the plated holes and the elements), the dielectric of the stack
// under the window, and the settings of the simulation script. Description:
// docs/solvers/openems.md.
#include "kairoq/openems.hpp"
#include "kairoq/embedded.hpp"

#include "emi_loop.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

using nlohmann::json;
namespace fs = std::filesystem;

constexpr double kC0 = 299792458.0;

// Polygons of CSXCAD have no holes: a polygon with holes is cut into strips by vertical
// lines, one through each hole, so that every hole opens into a strip's boundary. A piece
// that keeps a hole (a cut touching it at a vertex only) is cut again; after four levels it
// is filled, counted in `filled`.
std::vector<Path64> hole_free(const Polygon& poly, std::size_t& filled, int depth = 0) {
    if (poly.holes.empty()) return {poly.outer};
    const Rect64 ob = bounds(poly);
    std::set<int64_t> cuts;
    for (const auto& h : poly.holes) {
        const Rect64 hb = bounds(Polygon{h, {}});
        cuts.insert((hb.left + hb.right) / 2);
    }
    std::vector<int64_t> xs{ob.left - 1};
    xs.insert(xs.end(), cuts.begin(), cuts.end());
    xs.push_back(ob.right + 1);
    std::vector<Path64> out;
    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        const Path64 strip{Point64(xs[i], ob.top - 1), Point64(xs[i + 1], ob.top - 1), Point64(xs[i + 1], ob.bottom + 1),
                           Point64(xs[i], ob.bottom + 1)};
        for (const auto& piece : intersection(PolygonSet{poly}, union_polygons({strip}))) {
            if (piece.holes.empty()) {
                out.push_back(piece.outer);
            } else if (depth > 4) {
                out.push_back(piece.outer);
                filled += piece.holes.size();
            } else {
                for (auto& p : hole_free(piece, filled, depth + 1)) out.push_back(std::move(p));
            }
        }
    }
    return out;
}

// Vertices in metres, without repeated or collinear ones: CSXCAD's polygon filling leaves
// parts of a polygon with collinear vertices empty.
json path_m(const Path64& p) {
    std::vector<Point64> v;
    for (const auto& q : p)
        if (v.empty() || q != v.back()) v.push_back(q);
    while (v.size() > 1 && v.front() == v.back()) v.pop_back();
    for (bool changed = true; changed && v.size() > 3;) {
        changed = false;
        for (std::size_t i = 0; i < v.size() && v.size() > 3; ++i) {
            const Point64 a = v[(i + v.size() - 1) % v.size()], b = v[i], c = v[(i + 1) % v.size()];
            const double cross = double(b.x - a.x) * double(c.y - b.y) - double(b.y - a.y) * double(c.x - b.x);
            const double len = std::hypot(double(b.x - a.x), double(b.y - a.y)) * std::hypot(double(c.x - b.x), double(c.y - b.y));
            if (std::abs(cross) <= 1e-9 * len) {
                v.erase(v.begin() + static_cast<std::ptrdiff_t>(i));
                changed = true;
            }
        }
    }
    json a = json::array();
    for (const auto& q : v) a.push_back({q.x * 1e-9, q.y * 1e-9});
    return a;
}

// The dielectric of the stack as slabs: the dielectric layers, and the copper layers with
// the permittivity of the dielectric next to them (towards the inside); the outer copper
// layers as air unless outer_filled.
json dielectric_slabs(const Board& board, double default_eps, bool outer_filled = true) {
    const auto zs = stackup_z(board);
    auto eps_of = [&](const StackupLayer& l) { return l.epsilon_r ? *l.epsilon_r : default_eps; };
    json slabs = json::array();
    const std::size_t nl = board.stackup.size();
    for (std::size_t i = 0; i < nl; ++i) {
        const auto& l = board.stackup[i];
        const StackupLayer* d = nullptr;
        auto search = [&](bool down) {
            for (std::size_t k = i; !d && (down ? k + 1 < nl : k > 0);) {
                k = down ? k + 1 : k - 1;
                if (board.stackup[k].type == LayerType::Dielectric) d = &board.stackup[k];
            }
        };
        if (l.type == LayerType::Dielectric) d = &l;
        else if (l.type == LayerType::Copper) {
            if (!outer_filled && (l.copper_index == 1 || l.copper_index == static_cast<int>(board.copper.size()))) continue;
            const bool lower = l.copper_index > static_cast<int>(board.copper.size()) / 2;
            search(!lower);
            search(lower);
        }
        if (!d || !(zs[i].top_nm > zs[i].bottom_nm)) continue;
        // Adjacent layers of one material are one slab.
        const double eps = eps_of(*d), tand = d->loss_tangent ? *d->loss_tangent : 0.0;
        if (!slabs.empty() && slabs.back()["epsilon_r"] == eps && slabs.back()["loss_tangent"] == tand &&
            std::abs(slabs.back()["z0"].get<double>() - zs[i].top_nm * 1e-9) < 1e-12) {
            slabs.back()["z0"] = zs[i].bottom_nm * 1e-9;
            slabs.back()["name"] = slabs.back()["name"].get<std::string>() + "+" + l.name;
            continue;
        }
        slabs.push_back({{"name", l.name},
                         {"z0", zs[i].bottom_nm * 1e-9},
                         {"z1", zs[i].top_nm * 1e-9},
                         {"epsilon_r", eps_of(*d)},
                         {"loss_tangent", d->loss_tangent ? *d->loss_tangent : 0.0}});
    }
    return slabs;
}

// Writes the embedded script `name` into the directory as `as`.
void write_script(const std::string& dir, const std::string& name, const std::string& as) {
    for (std::size_t i = 0; i < openems_files_count; ++i) {
        if (openems_files[i].path != name) continue;
        std::ofstream f(fs::path(dir) / as, std::ios::binary);
        f.write(reinterpret_cast<const char*>(openems_files[i].data), static_cast<std::streamsize>(openems_files[i].size));
        if (!f) throw std::runtime_error("openems: cannot write " + (fs::path(dir) / as).string());
        return;
    }
    throw std::runtime_error("openems: no embedded script " + name);
}

}  // namespace

json openems_model(const Board& board, const LoopSpec& loop, const EmiSettings& settings) {
    const OpenEmsSettings& s = settings.openems;
    if (!(s.f_max_hz > 0) || !(s.min_cell_nm > 0) || !(s.board_cell_nm >= s.min_cell_nm))
        throw std::runtime_error("openems: f_max_hz and min_cell_mm must be positive, board_cell_mm at least min_cell_mm");
    // The model holds the window the estimate settled on (the window grows with the moment).
    std::optional<LoopResult> est;
    LoopSettings settled = settings.loop;
    if (s.estimate) {
        est = solve_loop(board, loop, settings.loop);
        settled.window_margin_nm = est->window_margin_mm * 1e6;
    }
    const detail::LoopGeometry g = detail::loop_geometry(board, loop, settled);
    json m;
    m["format"] = "kairoq-openems 1";
    m["board"] = board.name;
    m["loop"] = loop.name;
    m["units"] = "m";
    json warnings = json::array();
    std::size_t filled = 0;

    // Copper: thin perfectly conducting sheets at the layers' mid-planes.
    json sheets = json::array();
    std::size_t npoly = 0;
    for (const auto& [layer, islands] : g.copper) {
        json polys = json::array();
        for (const auto& isl : islands)
            for (const auto& p : hole_free(isl, filled)) polys.push_back(path_m(p));
        npoly += polys.size();
        sheets.push_back({{"layer", layer}, {"z", g.z_m[layer]}, {"polygons", polys}});
    }
    m["sheets"] = sheets;

    // Dielectric under the window: the dielectric layers, and the copper layers with the
    // permittivity of the dielectric next to them (towards the inside).
    json window = json::array();
    for (const auto& w : g.window)
        for (const auto& p : hole_free(w, filled)) window.push_back(path_m(p));
    m["window"] = window;
    const json slabs = dielectric_slabs(board, settings.dielectric_epsilon_r);
    m["dielectric"] = slabs;

    json vias = json::array();
    for (const auto& v : g.vias)
        vias.push_back({{"x", v.at.x * 1e-9}, {"y", v.at.y * 1e-9}, {"radius", v.radius_nm * 1e-9}, {"z0", g.z_m[v.from_layer]},
                        {"z1", g.z_m[v.to_layer]}});
    m["vias"] = vias;

    // Elements: straight conductors between the pad centres; the first is the port.
    json elements = json::array();
    for (std::size_t k = 0; k < g.elements.size(); ++k) {
        const auto& e = g.elements[k];
        const bool along_x = std::abs(e.to.y - e.from.y) <= 1000, along_y = std::abs(e.to.x - e.from.x) <= 1000;
        const std::string name = e.from_pad + " - " + e.to_pad;
        if (k == 0 && along_x == along_y)
            throw std::runtime_error("openems: the port element " + name + " must lie along x or along y");
        if (k > 0 && along_x == along_y)
            warnings.push_back("element " + name + " lies along neither x nor y; it is a thin wire on the mesh lines");
        elements.push_back({{"name", name},
                            {"from", {e.from.x * 1e-9, e.from.y * 1e-9}},
                            {"to", {e.to.x * 1e-9, e.to.y * 1e-9}},
                            {"z_sheet", g.z_m[e.layer]},
                            {"z_axis", e.z_m},
                            {"radius", e.radius_m},
                            {"port", k == 0},
                            {"capacitance", e.capacitance_f}});
        if (e.stitch) {
            const auto& st = loop.stitches[k - loop.elements.size()];
            if (st.capacitance_f > 0 && along_x == along_y)
                throw std::runtime_error("openems: the stitch " + name + " with a capacitance must lie along x or along y");
            if (st.extra_inductance_nh > 0 || st.resistance_ohm > 0)
                warnings.push_back("the series inductance and resistance of the stitch " + name + " are not part of the full-wave model");
            continue;
        }
        if (loop.elements[k].extra_inductance_nh > 0)
            warnings.push_back("the package inductance of " + name + " (" + std::to_string(loop.elements[k].extra_inductance_nh) +
                               " nH) is not part of the full-wave model; the estimate leaves it out");
    }
    m["elements"] = elements;

    std::vector<double> ff = s.far_field_hz;
    if (ff.empty()) ff = {s.f_max_hz / 10, s.f_max_hz / 4, s.f_max_hz / 2};
    for (double f : ff)
        if (!(f > 0) || f > s.f_max_hz) throw std::runtime_error("openems: far-field frequencies must lie in (0, f_max]");
    const double lambda = kC0 / s.f_max_hz;
    m["simulation"] = {{"f_max_hz", s.f_max_hz},
                       {"min_cell", s.min_cell_nm * 1e-9},
                       {"board_cell", s.board_cell_nm * 1e-9},
                       {"max_cell", lambda / 20},
                       {"air_margin", s.air_margin_nm > 0 ? s.air_margin_nm * 1e-9 : lambda / 4},
                       {"port_resistance_ohm", s.port_resistance_ohm},
                       {"far_field_hz", ff},
                       {"distance_m", settings.loop.distance_m},
                       {"end_criterion_db", s.end_criterion_db},
                       {"max_timesteps", s.max_timesteps}};

    if (est) {
        const LoopResult& r = *est;
        // The far field of the solved currents in free space, as the near-to-far-field
        // transformation of the simulation gives it.
        LoopSettings free = settings.loop;
        free.field_site = "free";
        json far = json::array();
        for (double f : ff) far.push_back({{"frequency_hz", f}, {"v_per_m_per_a", loop_field_per_ampere(r.currents, r.moment_m2, f, free, true)}});
        m["estimate"] = {{"inductance_nh", r.inductance_nh - r.element_inductance_nh},
                         {"inductance_coarse_nh", r.inductance_coarse_nh > 0 ? r.inductance_coarse_nh - r.element_inductance_nh : 0.0},
                         {"moment_m2", r.moment_m2},
                         {"far_field", far},
                         {"triangles", r.triangles}};
        for (const auto& w : r.warnings) warnings.push_back(w);
    }
    m["counts"] = {{"sheet_polygons", npoly}, {"vias", vias.size()}, {"elements", elements.size()}};
    if (filled) warnings.push_back(std::to_string(filled) + " holes of polygons could not be opened for CSXCAD and are filled");
    m["warnings"] = warnings;
    return m;
}

json export_openems(const Board& board, const LoopSpec& loop, const EmiSettings& settings, const std::string& dir) {
    const json m = openems_model(board, loop, settings);
    fs::create_directories(dir);
    {
        std::ofstream f(fs::path(dir) / "model.json");
        f << m.dump(1) << "\n";
        if (!f) throw std::runtime_error("openems: cannot write " + (fs::path(dir) / "model.json").string());
    }
    write_script(dir, "simulate.py", "simulate.py");
    return m;
}

std::string openems_summary(const json& m) {
    std::ostringstream o;
    const auto& s = m["simulation"];
    o << "loop         " << m["loop"].get<std::string>() << " on " << m["board"].get<std::string>() << "\n";
    o << "geometry     " << m["counts"]["sheet_polygons"].get<std::size_t>() << " copper polygons on " << m["sheets"].size()
      << " layers, " << m["vias"].size() << " plated holes, " << m["elements"].size() << " elements (the first is the port)\n";
    char line[256];
    std::snprintf(line, sizeof line, "simulation   0 to %.4g GHz, cells %.3g to %.3g mm, air %.3g mm, port %.4g ohm\n",
                  s["f_max_hz"].get<double>() * 1e-9, s["min_cell"].get<double>() * 1e3, s["max_cell"].get<double>() * 1e3,
                  s["air_margin"].get<double>() * 1e3, s["port_resistance_ohm"].get<double>());
    o << line;
    if (m.contains("estimate")) {
        const auto& e = m["estimate"];
        const auto mo = e["moment_m2"].get<std::array<double, 3>>();
        std::snprintf(line, sizeof line, "estimate     %.4g nH in the copper and elements, moment %.4g mm^2 per ampere\n",
                      e["inductance_nh"].get<double>(), std::sqrt(mo[0] * mo[0] + mo[1] * mo[1] + mo[2] * mo[2]) * 1e6);
        o << line;
    }
    for (const auto& w : m["warnings"]) o << "warning      " << w.get<std::string>() << "\n";
    return o.str();
}

// ------------------------------------------------------------------ channels

json openems_channel_model(const Board& board, const SiSettings& si, const ChannelSpec& ch, const OpenEmsChannelSettings& s) {
    if (!(s.f_max_hz > 0) || s.points < 2 || !(s.margin_nm >= 0))
        throw std::runtime_error("openems: f_max_hz positive, points at least 2, margin not negative");
    if (ch.ports.empty()) throw std::runtime_error("openems: channel " + ch.name + " has no ports");
    json m;
    m["format"] = "kairoq-openems-ports 1";
    m["board"] = board.name;
    m["channel"] = ch.name;
    m["units"] = "m";
    json warnings = json::array();
    std::size_t filled = 0;

    // The dielectric with the SI analysis's materials, at the middle of the band.
    Board stack = board;
    double thinnest = INFINITY, densest = 1.0;
    for (auto& l : stack.stackup) {
        if (l.type != LayerType::Dielectric) continue;
        const auto it = si.layer_materials.find(l.name);
        DielectricModel d = it != si.layer_materials.end() ? it->second : si.default_dielectric;
        if (it == si.layer_materials.end() && l.epsilon_r) {
            d.epsilon_r = *l.epsilon_r;
            if (l.loss_tangent) d.loss_tangent = *l.loss_tangent;
        }
        const cplx e = d.wideband ? d.at(s.f_max_hz / 2) : cplx(d.epsilon_r, -d.epsilon_r * d.loss_tangent);
        l.epsilon_r = e.real();
        l.loss_tangent = std::abs(e.imag()) / e.real();
        if (l.thickness_nm > 0) thinnest = std::min(thinnest, l.thickness_nm);
        densest = std::max(densest, e.real());
    }
    const double eps_default = si.default_dielectric.epsilon_r;
    const json slabs = dielectric_slabs(stack, eps_default, false);

    // Sheets: outer layers at their face towards the board, inner layers at mid-plane.
    const auto zs = stackup_z(board);
    std::map<int, double> z_of;
    for (std::size_t i = 0; i < board.stackup.size(); ++i) {
        const auto& l = board.stackup[i];
        if (l.type != LayerType::Copper) continue;
        const int n = static_cast<int>(board.copper.size());
        z_of[l.copper_index] = 1e-9 * (l.copper_index == 1 ? zs[i].bottom_nm : l.copper_index == n ? zs[i].top_nm : 0.5 * (zs[i].bottom_nm + zs[i].top_nm));
    }

    // Ports' pads and the box around their nets' copper.
    auto find_pad = [&](const std::string& ref) -> const Pad& {
        const auto dot = ref.rfind('.');
        if (dot == std::string::npos) throw std::runtime_error("openems: port " + ref + " is not REF.PIN");
        const std::string comp = ref.substr(0, dot), pin = ref.substr(dot + 1);
        for (const auto& p : board.pads)
            if (p.component == comp && p.pin == pin) return p;
        throw std::runtime_error("openems: no pad " + ref);
    };
    std::set<int> nets;
    for (const auto& p : ch.ports) nets.insert(find_pad(p).net);
    if (!ch.parts.empty()) warnings.push_back("the channel's parts are not in the model: their pads are left open");
    if (!ch.loads.empty()) warnings.push_back("the channel's loads are not in the model");
    Rect64 box{INT64_MAX, INT64_MAX, INT64_MIN, INT64_MIN};
    for (const auto& layer : board.copper)
        for (const auto& nc : layer.nets)
            if (nets.count(nc.net) && nc.net != 0) {
                const Rect64 b = bounds(nc.polygons);
                box.left = std::min(box.left, b.left);
                box.top = std::min(box.top, b.top);
                box.right = std::max(box.right, b.right);
                box.bottom = std::max(box.bottom, b.bottom);
            }
    if (box.left > box.right) throw std::runtime_error("openems: the ports' nets have no copper");
    const auto mg = static_cast<std::int64_t>(s.margin_nm);
    box = {box.left - mg, box.top - mg, box.right + mg, box.bottom + mg};
    PolygonSet window = {Polygon{Path64{{box.left, box.top}, {box.right, box.top}, {box.right, box.bottom}, {box.left, box.bottom}}, {}}};
    if (!board.outline.empty()) window = intersection(window, board.outline);
    json win = json::array();
    for (const auto& w : window)
        for (const auto& p : hole_free(w, filled)) win.push_back(path_m(p));
    m["window"] = win;
    m["dielectric"] = slabs;

    json sheets = json::array();
    std::size_t npoly = 0;
    for (const auto& layer : board.copper) {
        json polys = json::array(), signal = json::array();
        for (const auto& nc : layer.nets)
            for (const auto& poly : intersection(nc.polygons, window))
                for (const auto& p : hole_free(poly, filled)) {
                    polys.push_back(path_m(p));
                    signal.push_back(nets.count(nc.net) > 0 && nc.net != 0);
                }
        npoly += polys.size();
        if (polys.empty()) continue;
        json sh = {{"layer", layer.index}, {"z", z_of.at(layer.index)}, {"polygons", polys}, {"signal", signal}};
        // The copper's conductivity and thickness, for lossy sheets.
        for (const auto& l : board.stackup)
            if (l.type == LayerType::Copper && l.copper_index == layer.index) {
                sh["conductivity"] = l.conductivity_s_per_m ? *l.conductivity_s_per_m : si.copper.conductivity_s_per_m;
                sh["thickness_m"] = l.thickness_nm > 0 ? l.thickness_nm * 1e-9 : 35e-6;
            }
        sheets.push_back(sh);
    }
    m["sheets"] = sheets;

    // Vias; around those of the channel's nets, the circles their pads and the holes of the
    // other copper (antipads) draw on each layer, which the mesh resolves.
    auto ring_distance = [](const Path64& ring, Point64 at) {
        double best = INFINITY;
        for (std::size_t k = 0; k < ring.size(); ++k) {
            const Point64 p = ring[k], q = ring[(k + 1) % ring.size()];
            const double ex = double(q.x - p.x), ey = double(q.y - p.y), len2 = ex * ex + ey * ey;
            const double t = len2 > 0 ? std::clamp((double(at.x - p.x) * ex + double(at.y - p.y) * ey) / len2, 0.0, 1.0) : 0.0;
            best = std::min(best, std::hypot(double(p.x) + t * ex - double(at.x), double(p.y) + t * ey - double(at.y)));
        }
        return best;
    };
    json vias = json::array();
    for (const auto& d : board.drills) {
        if (!d.plated || d.from_layer == d.to_layer || !contains(window, d.at)) continue;
        json v = {{"x", d.at.x * 1e-9}, {"y", d.at.y * 1e-9}, {"radius", 0.5 * d.diameter_nm * 1e-9},
                  {"z0", z_of.at(d.to_layer)}, {"z1", z_of.at(d.from_layer)}};
        if (nets.count(d.net) && d.net != 0) {
            std::set<long> rings;  // um
            for (const auto& layer : board.copper)
                for (const auto& nc : layer.nets)
                    for (const auto& poly : nc.polygons) {
                        const Rect64 pb = bounds(poly);
                        if (d.at.x < pb.left || d.at.x > pb.right || d.at.y < pb.top || d.at.y > pb.bottom) continue;
                        if (nc.net == d.net) {
                            if (!contains(poly, d.at)) continue;
                            double r = ring_distance(poly.outer, d.at);
                            for (const auto& h : poly.holes) r = std::min(r, ring_distance(h, d.at));
                            if (r > 0.5 * d.diameter_nm * 1.02 && r < 3e6) rings.insert(std::lround(r * 1e-3));
                        } else {
                            for (const auto& h : poly.holes)
                                if (contains(Polygon{h, {}}, d.at)) {
                                    const double r = ring_distance(h, d.at);
                                    if (r < 3e6) rings.insert(std::lround(r * 1e-3));
                                }
                        }
                    }
            json rs = json::array();
            for (long r : rings) rs.push_back(r * 1e-6);
            v["rings"] = rs;
        }
        vias.push_back(v);
    }
    m["vias"] = vias;

    // Ports: a vertical sheet from the pad's layer to the nearest copper of another net
    // beneath (or above) the pad's centre, across the track that ends at the pad.
    json ports = json::array();
    for (std::size_t k = 0; k < ch.ports.size(); ++k) {
        const Pad& pad = find_pad(ch.ports[k]);
        const double zp = z_of.at(pad.layer);
        double zr = NAN;
        for (const auto& layer : board.copper) {
            if (layer.index == pad.layer) continue;
            bool covered = false;
            for (const auto& nc : layer.nets)
                if (nc.net != pad.net && contains(nc.polygons, pad.at)) covered = true;
            const double z = z_of.at(layer.index);
            if (covered && (std::isnan(zr) || std::abs(z - zp) < std::abs(zr - zp))) zr = z;
        }
        if (std::isnan(zr)) throw std::runtime_error("openems: port " + ch.ports[k] + " has no copper of another net beneath it");
        // The track ending nearest the pad's centre gives the direction; the sheet spans
        // across it, as wide as the track.
        const Track* best = nullptr;
        double dbest = INFINITY;
        bool at_start = true;
        for (const auto& t : board.tracks) {
            if (t.layer != pad.layer || t.net != pad.net) continue;
            for (const bool st : {true, false}) {
                const Point64 q = st ? t.start : t.end;
                const double d = std::hypot(double(q.x - pad.at.x), double(q.y - pad.at.y));
                if (d < dbest) {
                    dbest = d;
                    best = &t;
                    at_start = st;
                }
            }
        }
        const Rect64 pb = bounds(pad.shape);
        int span = 0;  // axis the sheet spans: 0 x, 1 y
        double width = 1e-9 * double(std::min(pb.right - pb.left, pb.bottom - pb.top));
        json line;  // a straight track from the pad along x or y: a microstrip port can feed it
        if (best && dbest < 0.5 * double(std::max(pb.right - pb.left, pb.bottom - pb.top)) + best->width_nm) {
            const Point64 a = at_start ? best->start : best->end, b = at_start ? best->end : best->start;
            const double dx = double(b.x - a.x), dy = double(b.y - a.y);
            span = std::abs(dx) >= std::abs(dy) ? 1 : 0;  // the track runs along x: the sheet spans y
            width = best->width_nm * 1e-9;
            if (std::min(std::abs(dx), std::abs(dy)) > 0.05 * std::hypot(dx, dy))
                warnings.push_back("port " + ch.ports[k] + ": the track is not along x or y; the port spans the nearer axis");
            else if (!best->centre)
                line = {{"axis", span == 1 ? 0 : 1}, {"dir", (span == 1 ? dx : dy) > 0 ? 1 : -1}, {"length", 1e-9 * std::hypot(dx, dy)}};
        } else {
            warnings.push_back("port " + ch.ports[k] + ": no track ends at the pad; the port spans x, as wide as the pad");
        }
        json pj = {{"name", ch.ports[k]}, {"x", pad.at.x * 1e-9}, {"y", pad.at.y * 1e-9}, {"z_pad", zp}, {"z_ref", zr},
                   {"span", span}, {"width", width}, {"resistance_ohm", ch.impedance_ohm}};
        if (!line.is_null()) pj["line"] = line;
        ports.push_back(pj);
    }
    m["ports"] = ports;

    const double lambda = kC0 / s.f_max_hz;
    const double min_cell = s.min_cell_nm > 0 ? s.min_cell_nm * 1e-9 : std::min(1e-4, std::isfinite(thinnest) ? 0.25e-9 * thinnest : 1e-4);
    const double board_cell = s.board_cell_nm > 0 ? s.board_cell_nm * 1e-9 : std::min(lambda / std::sqrt(densest) / 20, 4 * min_cell);
    if (board_cell < min_cell) throw std::runtime_error("openems: board_cell below min_cell");
    std::vector<double> ff = s.far_field_hz;
    for (double f : ff)
        if (!(f > 0) || f > s.f_max_hz) throw std::runtime_error("openems: far-field frequencies must lie in (0, f_max]");
    m["simulation"] = {{"f_max_hz", s.f_max_hz},
                       {"points", s.points},
                       {"min_cell", min_cell},
                       {"board_cell", board_cell},
                       {"max_cell", lambda / 20},
                       {"air_margin", s.air_margin_nm > 0 ? s.air_margin_nm * 1e-9 : lambda / 4},
                       {"far_field_hz", ff},
                       {"end_criterion_db", s.end_criterion_db},
                       {"max_timesteps", s.max_timesteps}};
    const double area_mm2 = 1e-12 * double(box.right - box.left) * double(box.bottom - box.top);
    if (area_mm2 > 2500) warnings.push_back("the model spans " + std::to_string(int(area_mm2)) + " mm2; full-wave runs grow with its volume");
    m["counts"] = {{"sheet_polygons", npoly}, {"vias", vias.size()}, {"ports", ports.size()}};
    if (filled) warnings.push_back(std::to_string(filled) + " holes of polygons could not be opened for CSXCAD and are filled");
    m["warnings"] = warnings;
    return m;
}

json export_openems_channel(const Board& board, const SiSettings& si, const ChannelSpec& ch, const OpenEmsChannelSettings& s,
                            const std::string& dir) {
    const json m = openems_channel_model(board, si, ch, s);
    fs::create_directories(dir);
    std::ofstream f(fs::path(dir) / "model.json");
    f << m.dump(1) << "\n";
    if (!f) throw std::runtime_error("openems: cannot write " + (fs::path(dir) / "model.json").string());
    write_script(dir, "simulate_ports.py", "simulate.py");
    return m;
}

std::string openems_channel_summary(const json& m) {
    std::ostringstream o;
    const auto& s = m["simulation"];
    o << "channel      " << m["channel"].get<std::string>() << " on " << m["board"].get<std::string>() << "\n";
    o << "geometry     " << m["counts"]["sheet_polygons"].get<std::size_t>() << " copper polygons on " << m["sheets"].size() << " layers, "
      << m["vias"].size() << " plated holes, " << m["ports"].size() << " ports\n";
    char line[256];
    std::snprintf(line, sizeof line, "simulation   0 to %.4g GHz, cells %.3g to %.3g mm, air %.3g mm\n", s["f_max_hz"].get<double>() * 1e-9,
                  s["min_cell"].get<double>() * 1e3, s["max_cell"].get<double>() * 1e3, s["air_margin"].get<double>() * 1e3);
    o << line;
    for (const auto& p : m["ports"]) {
        std::snprintf(line, sizeof line, "port         %s: %.4g ohm, %.3g mm high, %.3g mm wide\n", p["name"].get<std::string>().c_str(),
                      p["resistance_ohm"].get<double>(), std::abs(p["z_pad"].get<double>() - p["z_ref"].get<double>()) * 1e3,
                      p["width"].get<double>() * 1e3);
        o << line;
    }
    for (const auto& w : m["warnings"]) o << "warning      " << w.get<std::string>() << "\n";
    return o.str();
}

}  // namespace kairoq
