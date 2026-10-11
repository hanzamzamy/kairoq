// PDN impedance (pdn.hpp, docs/solvers/pdn.md).
//
// Each pair of a power plane and its nearest ground plane is a cavity over the region both
// cover (cavity.hpp), with a port at every via that joins a load, a capacitor or the
// regulator to the power plane. The pairs' admittance matrices add at the ports, the
// branches are shunts there, and the loads' impedances follow from the inverse.
#include "kairoq/pdn.hpp"

#include "cavity.hpp"
#include "fem_util.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

using namespace detail;
using cd = std::complex<double>;

constexpr double kMu0 = 1.25663706212e-6;
constexpr double kC0 = 299792458.0;

// Inverse of a dense complex matrix (row major) by Gauss-Jordan elimination with partial
// pivoting.
std::vector<cd> invert(std::vector<cd> a, std::size_t n) {
    std::vector<cd> r(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) r[i * n + i] = 1.0;
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t i = c + 1; i < n; ++i)
            if (std::abs(a[i * n + c]) > std::abs(a[p * n + c])) p = i;
        if (std::abs(a[p * n + c]) == 0.0) throw std::runtime_error("pdn: singular admittance matrix");
        for (std::size_t k = 0; k < n; ++k) {
            std::swap(a[c * n + k], a[p * n + k]);
            std::swap(r[c * n + k], r[p * n + k]);
        }
        const cd d = 1.0 / a[c * n + c];
        for (std::size_t k = 0; k < n; ++k) {
            a[c * n + k] *= d;
            r[c * n + k] *= d;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == c || a[i * n + c] == 0.0) continue;
            const cd f = a[i * n + c];
            for (std::size_t k = 0; k < n; ++k) {
                a[i * n + k] -= f * a[c * n + k];
                r[i * n + k] -= f * r[c * n + k];
            }
        }
    }
    return r;
}

PolygonSet net_copper(const Board& b, int layer, int net) {
    for (const auto& cl : b.copper)
        if (cl.index == layer)
            for (const auto& nc : cl.nets)
                if (nc.net == net) return nc.polygons;
    return {};
}

// Surface impedance of a copper plate of conductivity sigma and thickness t, current on one
// face: (1 + j) / (sigma delta) coth((1 + j) t / delta), the sheet resistance 1 / (sigma t)
// at low frequency.
cd surface_impedance(double f, double sigma, double t) {
    const double delta = 1.0 / std::sqrt(kPi * f * kMu0 * sigma);
    const cd x = cd(1, 1) * t / delta;
    if (std::abs(x) < 1e-3) return 1.0 / (sigma * t);
    const cd e = std::exp(-2.0 * x);
    return cd(1, 1) / (sigma * delta) * (1.0 + e) / (1.0 - e);
}

// Inductance of a path of round wires that leaves a perfectly conducting plane at height
// z_plane and returns to it (the ends on the plane): half that of the closed loop of the
// path and its image (Neumann's formula, Gauss-Legendre on each pair of straight pieces;
// collinear consecutive pieces of one radius merged, partial self inductances in closed
// form).
double image_loop_inductance(const std::vector<std::array<double, 3>>& path, const std::vector<double>& radii, double z_plane) {
    std::vector<std::array<double, 3>> pts = path;
    std::vector<double> rad = radii;
    for (std::size_t k = path.size(); k-- > 0;) pts.push_back({path[k][0], path[k][1], 2 * z_plane - path[k][2]});
    for (std::size_t k = radii.size(); k-- > 0;) rad.push_back(radii[k]);
    struct Piece {
        std::array<double, 3> a, b;
        double r;
    };
    std::vector<Piece> pieces;
    for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
        const std::array<double, 3> a = pts[k], b = pts[k + 1];
        const double len = std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        if (len < 1e-12) continue;
        const double r = rad[std::min(k, rad.size() - 1)];
        if (!pieces.empty()) {
            auto& q = pieces.back();
            const std::array<double, 3> u{q.b[0] - q.a[0], q.b[1] - q.a[1], q.b[2] - q.a[2]}, v{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const double lu = std::hypot(u[0], u[1], u[2]);
            const double cosang = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (lu * len);
            if (cosang > 1 - 1e-12 && std::abs(q.r - r) < 1e-12 * r) {
                q.b = b;
                continue;
            }
        }
        pieces.push_back({a, b, r});
    }
    static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
    static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    double L = 0;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        const auto& p = pieces[i];
        const std::array<double, 3> u{p.b[0] - p.a[0], p.b[1] - p.a[1], p.b[2] - p.a[2]};
        const double lp = std::hypot(u[0], u[1], u[2]);
        L += kMu0 / (2 * kPi) * (lp * std::asinh(lp / p.r) - std::hypot(lp, p.r) + p.r);
        for (std::size_t j = 0; j < pieces.size(); ++j) {
            if (i == j) continue;
            const auto& q = pieces[j];
            const std::array<double, 3> v{q.b[0] - q.a[0], q.b[1] - q.a[1], q.b[2] - q.a[2]};
            const double dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            if (std::abs(dot) < 1e-30) continue;
            double sum = 0;
            // Two levels of 8 points each piece: 16 per piece.
            for (int hi = 0; hi < 2; ++hi)
                for (int ki = 0; ki < 8; ++ki) {
                    const double ti = 0.25 * (gx[ki] + 1) + 0.5 * hi;
                    const double xi[3] = {p.a[0] + ti * u[0], p.a[1] + ti * u[1], p.a[2] + ti * u[2]};
                    for (int hj = 0; hj < 2; ++hj)
                        for (int kj = 0; kj < 8; ++kj) {
                            const double tj = 0.25 * (gx[kj] + 1) + 0.5 * hj;
                            const double R = std::hypot(xi[0] - q.a[0] - tj * v[0], xi[1] - q.a[1] - tj * v[1], xi[2] - q.a[2] - tj * v[2]);
                            sum += 0.25 * gw[ki] * 0.25 * gw[kj] / std::max(R, 0.5 * (p.r + q.r));
                        }
                }
            L += kMu0 / (4 * kPi) * dot * sum;
        }
    }
    return 0.5 * L;
}

struct Pair {
    int power = 0, ground = 0;
    double upper_face = 0, lower_face = 0;  // nm, facing each other
    double outer_upper = 0, outer_lower = 0;  // nm, the planes' faces away from the pair
    PolygonSet region;
    std::vector<double> thickness_m;
    std::vector<DielectricModel> materials;
    double sigma[2] = {5.8e7, 5.8e7}, t_m[2] = {35e-6, 35e-6};
    std::vector<std::size_t> ports;  // global port indices inside the region
    std::shared_ptr<const Cavity> cavity;
};

struct Port {
    Point64 at;
    double radius_nm = 0;
    int from = 0, to = 0;  // layers the joining barrel spans
};

}  // namespace

PdnResult analyse_pdn(const Board& b, const PdnSettings& s) {
    const auto started = Clock::now();
    PdnResult res;
    const int pn = b.find_net(s.power_net), gn = b.find_net(s.ground_net);
    if (pn <= 0) throw std::runtime_error("pdn: no net " + s.power_net);
    if (gn <= 0) throw std::runtime_error("pdn: no net " + s.ground_net);
    if (!(s.f_min_hz > 0) || !(s.f_max_hz > s.f_min_hz) || s.points < 2) throw std::runtime_error("pdn: bad frequency range");
    const auto zs = stackup_z(b);
    std::map<int, std::size_t> stack_of;
    for (std::size_t i = 0; i < b.stackup.size(); ++i)
        if (b.stackup[i].type == LayerType::Copper) stack_of[b.stackup[i].copper_index] = i;
    auto zmid = [&](int layer) { return 0.5 * (zs[stack_of.at(layer)].top_nm + zs[stack_of.at(layer)].bottom_nm); };
    auto material = [&](std::size_t i) {
        DielectricModel m = s.dielectric;
        if (b.stackup[i].epsilon_r) m.epsilon_r = *b.stackup[i].epsilon_r;
        if (b.stackup[i].loss_tangent) m.loss_tangent = *b.stackup[i].loss_tangent;
        return m;
    };
    // Planes of the two nets.
    double board_area = area_nm2(b.outline);
    if (!(board_area > 0))
        for (const auto& cl : b.copper)
            for (const auto& nc : cl.nets) board_area = std::max(board_area, area_nm2(nc.polygons));
    std::vector<int> power_layers, ground_layers;
    for (const auto& cl : b.copper) {
        if (area_nm2(net_copper(b, cl.index, pn)) >= s.plane_fraction * board_area) power_layers.push_back(cl.index);
        if (area_nm2(net_copper(b, cl.index, gn)) >= s.plane_fraction * board_area) ground_layers.push_back(cl.index);
    }
    if (power_layers.empty()) throw std::runtime_error("pdn: " + s.power_net + " has no plane (copper over " + std::to_string(s.plane_fraction) + " of the board)");
    if (ground_layers.empty()) throw std::runtime_error("pdn: " + s.ground_net + " has no plane");
    // Pairs: each power plane with the nearest ground plane above and below.
    std::vector<Pair> pairs;
    for (int p : power_layers) {
        for (int side : {-1, 1}) {
            int best = 0;
            for (int g : ground_layers) {
                if ((zmid(g) - zmid(p)) * side <= 0) continue;
                if (best == 0 || std::abs(zmid(g) - zmid(p)) < std::abs(zmid(best) - zmid(p))) best = g;
            }
            if (best == 0) continue;
            Pair pr;
            pr.power = p;
            pr.ground = best;
            const int up = zmid(p) > zmid(best) ? p : best, lo = up == p ? best : p;
            pr.upper_face = zs[stack_of.at(up)].bottom_nm;
            pr.lower_face = zs[stack_of.at(lo)].top_nm;
            pr.outer_upper = zs[stack_of.at(up)].top_nm;
            pr.outer_lower = zs[stack_of.at(lo)].bottom_nm;
            pr.region = intersection(net_copper(b, p, pn), net_copper(b, best, gn));
            for (auto& poly : pr.region) {
                Paths64 kept;
                for (const auto& h : poly.holes)
                    if (std::abs(Clipper2Lib::Area(h)) > kPi * 1.25e6 * 1.25e6) kept.push_back(h);
                poly.holes = kept;
            }
            if (area_nm2(pr.region) < 0.01 * board_area) continue;
            for (std::size_t i = 0; i < b.stackup.size(); ++i)
                if (b.stackup[i].type == LayerType::Dielectric && zs[i].top_nm <= pr.upper_face + 1 && zs[i].bottom_nm >= pr.lower_face - 1) {
                    pr.thickness_m.push_back((zs[i].top_nm - zs[i].bottom_nm) * 1e-9);
                    pr.materials.push_back(material(i));
                }
            int k = 0;
            for (int l : {p, best}) {
                const auto& sl = b.stackup[stack_of.at(l)];
                if (sl.conductivity_s_per_m) pr.sigma[k] = *sl.conductivity_s_per_m;
                if (sl.thickness_nm > 0) pr.t_m[k] = sl.thickness_nm * 1e-9;
                ++k;
            }
            pairs.push_back(std::move(pr));
        }
    }
    if (pairs.empty()) throw std::runtime_error("pdn: no pair of " + s.power_net + " and " + s.ground_net + " planes");
    // Branches and the points where they join the power plane.
    std::vector<Port> ports;
    auto join = [&](const std::string& ref, int net, Port& out, double& z_pad) {
        double best = INFINITY;
        for (const auto& pad : b.pads) {
            if (pad.component != ref || pad.net != net) continue;
            // A pin on the plane layer, or the nearest via of the net within reach.
            for (int p : power_layers)
                if (net == pn && pad.layer == p && best > 0) {
                    best = 0;
                    out = {pad.at, std::sqrt(area_nm2(pad.shape) / kPi), p, p};
                    z_pad = zmid(pad.layer);
                }
            for (const auto& d : b.drills) {
                if (!d.plated || d.net != net) continue;
                const double dist = std::hypot(double(d.at.x - pad.at.x), double(d.at.y - pad.at.y));
                if (dist > s.via_reach_nm || dist >= best) continue;
                best = dist;
                out = {d.at, 0.5 * d.diameter_nm, d.from_layer, d.to_layer};
                z_pad = zmid(pad.layer);
            }
        }
        return std::isfinite(best);
    };
    struct Branch {
        PdnBranch info;
        std::size_t port = 0;
        double r = 0, l = 0, c = 0;  // series R, L (with mounting), C (0: none)
    };
    std::vector<Branch> branches;
    auto add_branch = [&](const std::string& ref, const std::string& kind, double r, double l, double c, std::optional<double> mount) {
        Port pt;
        double z_pad = 0;
        if (!join(ref, pn, pt, z_pad)) {
            res.warnings.push_back(ref + ": no pad of " + s.power_net + " joined to its plane; left out");
            return;
        }
        std::size_t idx = ports.size();
        for (std::size_t i = 0; i < ports.size(); ++i)
            if (ports[i].at == pt.at) idx = i;
        if (idx == ports.size()) ports.push_back(pt);
        Branch br;
        br.info.component = ref;
        br.info.kind = kind;
        br.info.x_mm = pt.at.x * 1e-6;
        br.info.y_mm = pt.at.y * 1e-6;
        br.port = idx;
        if (kind != "load") {
            // Mounting inductance: the loop of the pads, their tracks to the vias and the vias
            // down to the nearer face of the pair, over that face as a plane (images); the
            // part itself is its ESL.
            double m = 0;
            if (mount) {
                m = *mount;
            } else {
                Port gp;
                double zg = 0;
                Point64 power_pad = pt.at, ground_pad = pt.at;
                double pad_r = pt.radius_nm;
                for (const auto& pad : b.pads) {
                    if (pad.component != ref) continue;
                    const Rect64 bb = bounds(pad.shape);
                    const double half = 0.5 * std::min(double(bb.right - bb.left), double(bb.bottom - bb.top));
                    if (pad.net == pn) {
                        power_pad = Point64((bb.left + bb.right) / 2, (bb.top + bb.bottom) / 2);
                        pad_r = 0.5 * half;
                    }
                    if (pad.net == gn) ground_pad = Point64((bb.left + bb.right) / 2, (bb.top + bb.bottom) / 2);
                }
                const bool ground_via = join(ref, gn, gp, zg);
                if (!ground_via) gp = {Point64(power_pad.x + 1000000, power_pad.y), pt.radius_nm, 0, 0};
                double face = NAN, best = INFINITY;
                for (const auto& pr : pairs)
                    for (double f : {pr.outer_upper, pr.outer_lower})
                        if (std::abs(z_pad - f) < best) {
                            best = std::abs(z_pad - f);
                            face = f;
                        }
                const double s9 = 1e-9;
                const std::vector<std::array<double, 3>> path{{pt.at.x * s9, pt.at.y * s9, face * s9},
                                                              {pt.at.x * s9, pt.at.y * s9, z_pad * s9},
                                                              {power_pad.x * s9, power_pad.y * s9, z_pad * s9},
                                                              {ground_pad.x * s9, ground_pad.y * s9, z_pad * s9},
                                                              {gp.at.x * s9, gp.at.y * s9, z_pad * s9},
                                                              {gp.at.x * s9, gp.at.y * s9, face * s9}};
                const double rv = pt.radius_nm * s9, rg = gp.radius_nm * s9, rp = pad_r * s9;
                m = image_loop_inductance(path, {rv, rp, rp, rp, rg}, face * s9);
                if (!ground_via) res.warnings.push_back(ref + ": no via of " + s.ground_net + " near its pads; its ground via taken 1 mm away");
            }
            br.info.mount_h = m;
            br.r = r;
            br.l = l + m;
            br.c = c;
        }
        branches.push_back(br);
    };
    for (const auto& ld : s.loads) add_branch(ld, "load", 0, 0, 0, std::nullopt);
    for (const auto& c : s.capacitors) add_branch(c.component, "capacitor", c.esr_ohm, c.esl_h, c.c_f, c.mount_h);
    if (s.regulator) add_branch(s.regulator->component, "regulator", s.regulator->r_ohm, s.regulator->l_h, 0, std::nullopt);
    if (s.other_capacitors) {
        std::map<std::string, std::pair<bool, bool>> on;  // component -> pads on the power, ground net
        std::map<std::string, int> count;
        for (const auto& pad : b.pads) {
            if (pad.component.empty()) continue;
            ++count[pad.component];
            if (pad.net == pn) on[pad.component].first = true;
            if (pad.net == gn) on[pad.component].second = true;
        }
        for (const auto& [ref, both] : on) {
            if (!both.first || !both.second || count[ref] != 2 || ref.rfind(s.capacitor_prefix, 0) != 0) continue;
            bool listed = std::find(s.loads.begin(), s.loads.end(), ref) != s.loads.end() || (s.regulator && s.regulator->component == ref);
            for (const auto& c : s.capacitors) listed |= c.component == ref;
            if (listed) continue;
            add_branch(ref, "capacitor", s.other_capacitors->esr_ohm, s.other_capacitors->esl_h, s.other_capacitors->c_f, s.other_capacitors->mount_h);
        }
    }
    std::vector<std::size_t> load_port;
    for (const auto& br : branches)
        if (br.info.kind == "load") {
            res.loads.push_back(br.info.component);
            load_port.push_back(br.port);
        }
    if (load_port.empty()) throw std::runtime_error("pdn: no load joined to the planes");
    // Cavities: each pair with the ports inside its region that reach its power plane.
    double er_max = 1.0;
    for (const auto& pr : pairs)
        for (const auto& m : pr.materials) er_max = std::max(er_max, m.epsilon_r);
    for (auto& pr : pairs) {
        std::vector<CavityPort> cp;
        for (std::size_t i = 0; i < ports.size(); ++i)
            if (ports[i].from <= pr.power && ports[i].to >= pr.power && contains(pr.region, ports[i].at)) {
                pr.ports.push_back(i);
                cp.push_back({double(ports[i].at.x), double(ports[i].at.y), ports[i].radius_nm, false});
            }
        if (pr.ports.empty()) continue;
        const auto t_cavity = Clock::now();
        pr.cavity = std::make_shared<const Cavity>(pr.region, cp, 2 * kPi * s.f_max_hz * std::sqrt(er_max) / kC0);
        PdnPair info;
        info.power_layer = pr.power;
        info.ground_layer = pr.ground;
        info.height_m = (pr.upper_face - pr.lower_face) * 1e-9;
        info.area_m2 = pr.cavity->area_m2();
        info.modes = pr.cavity->modes();
        info.nodes = pr.cavity->nodes();
        info.lanczos_steps = pr.cavity->lanczos_steps();
        info.seconds = since(t_cavity);
        res.pairs.push_back(info);
    }
    for (std::size_t i = 0; i < ports.size(); ++i) {
        bool reached = false;
        for (const auto& pr : pairs) reached |= std::find(pr.ports.begin(), pr.ports.end(), i) != pr.ports.end();
        if (!reached) throw std::runtime_error("pdn: the via at (" + std::to_string(ports[i].at.x * 1e-6) + ", " + std::to_string(ports[i].at.y * 1e-6) +
                                               ") mm does not reach a pair of planes");
    }
    for (const auto& br : branches) res.branches.push_back(br.info);
    // Frequencies.
    const std::size_t np = ports.size();
    for (std::size_t k = 0; k < s.points; ++k)
        res.frequency_hz.push_back(s.f_min_hz * std::pow(s.f_max_hz / s.f_min_hz, double(k) / double(s.points - 1)));
    res.z_ohm.assign(load_port.size(), std::vector<cd>(s.points));
    std::vector<std::string> errors(s.points);
    parallel_for(s.points, [&](std::size_t k) {
        try {
            const double f = res.frequency_hz[k], w = 2 * kPi * f;
            std::vector<cd> Y(np * np, 0.0);
            for (const auto& pr : pairs) {
                if (!pr.cavity) continue;
                cd inv = 0.0;
                double t = 0;
                for (std::size_t i = 0; i < pr.materials.size(); ++i) {
                    inv += pr.thickness_m[i] / pr.materials[i].at(f);
                    t += pr.thickness_m[i];
                }
                const cd eps = t > 0 ? t / inv : cd(1.0);
                const cd zsurf = 0.5 * (surface_impedance(f, pr.sigma[0], pr.t_m[0]) + surface_impedance(f, pr.sigma[1], pr.t_m[1]));
                const std::size_t m = pr.ports.size();
                const auto yp = invert(pr.cavity->impedance(w, (pr.upper_face - pr.lower_face) * 1e-9, eps, zsurf), m);
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j) Y[pr.ports[i] * np + pr.ports[j]] += yp[i * m + j];
            }
            for (const auto& br : branches) {
                if (br.info.kind == "load") continue;
                const cd z = br.r + cd(0, w * br.l) + (br.c > 0 ? 1.0 / cd(0, w * br.c) : cd(0));
                Y[br.port * np + br.port] += 1.0 / z;
            }
            const auto Z = invert(Y, np);
            for (std::size_t l = 0; l < load_port.size(); ++l) res.z_ohm[l][k] = Z[load_port[l] * np + load_port[l]];
        } catch (const std::exception& e) {
            errors[k] = e.what();
        }
    });
    for (const auto& e : errors)
        if (!e.empty()) throw std::runtime_error(e);
    res.seconds = since(started);
    return res;
}

std::string pdn_summary(const PdnResult& r, const PdnSettings& s) {
    std::ostringstream o;
    char line[256];
    for (const auto& p : r.pairs) {
        std::snprintf(line, sizeof line, "pair         L%d (%s) and L%d (%s), %.3f mm apart, %.1f cm^2, %zu modes, %zu nodes\n", p.power_layer,
                      s.power_net.c_str(), p.ground_layer, s.ground_net.c_str(), p.height_m * 1e3, p.area_m2 * 1e4, p.modes, p.nodes);
        o << line;
    }
    for (const auto& b : r.branches) {
        if (b.kind == "load")
            std::snprintf(line, sizeof line, "%-12s %-10s at (%.3f, %.3f) mm\n", b.kind.c_str(), b.component.c_str(), b.x_mm, b.y_mm);
        else
            std::snprintf(line, sizeof line, "%-12s %-10s at (%.3f, %.3f) mm, mounting %.3f nH\n", b.kind.c_str(), b.component.c_str(), b.x_mm, b.y_mm,
                          b.mount_h * 1e9);
        o << line;
    }
    for (std::size_t l = 0; l < r.loads.size(); ++l) {
        const auto& z = r.z_ohm[l];
        std::size_t kmax = 0;
        for (std::size_t k = 0; k < z.size(); ++k)
            if (std::abs(z[k]) > std::abs(z[kmax])) kmax = k;
        std::snprintf(line, sizeof line, "load %s: |Z| %.4g ohm at %.4g Hz, %.4g ohm at %.4g Hz; largest %.4g ohm at %.4g Hz\n", r.loads[l].c_str(),
                      std::abs(z.front()), r.frequency_hz.front(), std::abs(z.back()), r.frequency_hz.back(), std::abs(z[kmax]), r.frequency_hz[kmax]);
        o << line;
        // Resonances: local maxima.
        std::string peaks;
        for (std::size_t k = 1; k + 1 < z.size(); ++k)
            if (std::abs(z[k]) > std::abs(z[k - 1]) && std::abs(z[k]) >= std::abs(z[k + 1])) {
                std::snprintf(line, sizeof line, "%s%.4g MHz (%.3g ohm)", peaks.empty() ? "" : ", ", r.frequency_hz[k] * 1e-6, std::abs(z[k]));
                peaks += line;
            }
        if (!peaks.empty()) o << "  peaks " << peaks << "\n";
        if (s.target_ohm) {
            std::string bands;
            for (std::size_t k = 0; k < z.size(); ++k) {
                if (!(std::abs(z[k]) > *s.target_ohm)) continue;
                const std::size_t k0 = k;
                while (k + 1 < z.size() && std::abs(z[k + 1]) > *s.target_ohm) ++k;
                std::snprintf(line, sizeof line, "%s%.4g to %.4g MHz", bands.empty() ? "" : ", ", r.frequency_hz[k0] * 1e-6, r.frequency_hz[k] * 1e-6);
                bands += line;
            }
            std::snprintf(line, sizeof line, "  target %.4g ohm: %s\n", *s.target_ohm, bands.empty() ? "met" : ("exceeded " + bands).c_str());
            o << line;
        }
    }
    std::snprintf(line, sizeof line, "time         %.2f s\n", r.seconds);
    o << line;
    for (const auto& w : r.warnings) o << "warning: " << w << "\n";
    return o.str();
}

std::string pdn_json(const PdnResult& r) {
    nlohmann::json j;
    j["frequency_hz"] = r.frequency_hz;
    j["loads"] = r.loads;
    nlohmann::json z = nlohmann::json::array();
    for (const auto& zl : r.z_ohm) {
        nlohmann::json v = nlohmann::json::array();
        for (const auto& c : zl) v.push_back({c.real(), c.imag()});
        z.push_back(v);
    }
    j["z_ohm"] = z;
    nlohmann::json pairs = nlohmann::json::array();
    for (const auto& p : r.pairs)
        pairs.push_back({{"power_layer", p.power_layer}, {"ground_layer", p.ground_layer}, {"height_m", p.height_m}, {"area_m2", p.area_m2},
                         {"modes", p.modes}, {"nodes", p.nodes}, {"lanczos_steps", p.lanczos_steps}, {"seconds", p.seconds}});
    j["pairs"] = pairs;
    nlohmann::json br = nlohmann::json::array();
    for (const auto& b : r.branches)
        br.push_back({{"component", b.component}, {"kind", b.kind}, {"x_mm", b.x_mm}, {"y_mm", b.y_mm}, {"mount_h", b.mount_h}});
    j["branches"] = br;
    j["seconds"] = r.seconds;
    j["warnings"] = r.warnings;
    return j.dump();
}

}  // namespace kairoq
