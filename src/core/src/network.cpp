// Networks of a board's nets as transmission lines (docs/solvers/networks.md).
//
// The tracks of the channel's nets are cut and solved as in the signal-integrity analysis;
// each sample stands for its length of track as a uniform line, single, or coupled with
// the matching sample of its partner where a pair runs together. A via is one barrel node
// with capacitances to the reference planes around it, whose nodes are joined by the
// impedance of each pair of planes (the physics-based via model); without planes, its
// barrel is a coaxial line between the copper layers it joins. Lines, parts, loads and vias
// are stamped into a nodal admittance matrix at each frequency, which gives the impedance
// matrix of the ports and from it the S-parameters.
#include "kairoq/network.hpp"

#include "fem_util.hpp"
#include "cavity.hpp"
#include "via_model.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace kairoq {

namespace {

using namespace detail;
namespace fs = std::filesystem;

constexpr double kC0 = 299792458.0;
constexpr std::size_t kGround = SIZE_MAX;  // the reference node of the matrix, not stored

// ------------------------------------------------------------------ dense complex algebra

struct Cm {
    std::size_t n = 0;
    std::vector<cplx> a;
    explicit Cm(std::size_t size = 0) : n(size), a(size * size, 0.0) {}
    cplx& operator()(std::size_t i, std::size_t j) { return a[i * n + j]; }
    const cplx& operator()(std::size_t i, std::size_t j) const { return a[i * n + j]; }
    static Cm identity(std::size_t size) {
        Cm m(size);
        for (std::size_t i = 0; i < size; ++i) m(i, i) = 1.0;
        return m;
    }
};

Cm mul(const Cm& x, const Cm& y) {
    Cm r(x.n);
    for (std::size_t i = 0; i < x.n; ++i)
        for (std::size_t k = 0; k < x.n; ++k) {
            const cplx v = x(i, k);
            if (v == 0.0) continue;
            for (std::size_t j = 0; j < x.n; ++j) r(i, j) += v * y(k, j);
        }
    return r;
}

Cm combine(const Cm& x, const Cm& y, cplx s) {  // x + s y
    Cm r = x;
    for (std::size_t i = 0; i < r.a.size(); ++i) r.a[i] += s * y.a[i];
    return r;
}

// Solves A X = B for X (A n x n, B n x k, row major), partial pivoting.
std::vector<cplx> solve(std::vector<cplx> A, std::size_t n, std::vector<cplx> B, std::size_t k) {
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::abs(A[r * n + c]) > std::abs(A[p * n + c])) p = r;
        if (std::abs(A[p * n + c]) == 0.0) throw std::runtime_error("network: singular admittance matrix");
        if (p != c) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A[c * n + j], A[p * n + j]);
            for (std::size_t j = 0; j < k; ++j) std::swap(B[c * k + j], B[p * k + j]);
        }
        const cplx d = A[c * n + c];
        for (std::size_t r = c + 1; r < n; ++r) {
            const cplx f = A[r * n + c] / d;
            if (f == 0.0) continue;
            for (std::size_t j = c; j < n; ++j) A[r * n + j] -= f * A[c * n + j];
            for (std::size_t j = 0; j < k; ++j) B[r * k + j] -= f * B[c * k + j];
        }
    }
    for (std::size_t c = n; c-- > 0;) {
        for (std::size_t j = 0; j < k; ++j) {
            cplx v = B[c * k + j];
            for (std::size_t m = c + 1; m < n; ++m) v -= A[c * n + m] * B[m * k + j];
            B[c * k + j] = v / A[c * n + c];
        }
    }
    return B;
}

Cm inverse(const Cm& m) {
    Cm id = Cm::identity(m.n);
    Cm r(m.n);
    r.a = solve(m.a, m.n, id.a, m.n);
    return r;
}

// Matrix exponential: scaling and squaring with the (6, 6) Pade approximant.
Cm expm(const Cm& m) {
    double norm = 0;
    for (std::size_t i = 0; i < m.n; ++i) {
        double row = 0;
        for (std::size_t j = 0; j < m.n; ++j) row += std::abs(m(i, j));
        norm = std::max(norm, row);
    }
    int s = norm > 0.5 ? static_cast<int>(std::ceil(std::log2(norm / 0.5))) : 0;
    Cm x = m;
    for (auto& v : x.a) v /= std::pow(2.0, s);
    static const double c[] = {1.0, 1.0 / 2, 5.0 / 44, 1.0 / 66, 1.0 / 792, 1.0 / 15840, 1.0 / 665280};
    Cm num = Cm::identity(m.n), den = Cm::identity(m.n), p = Cm::identity(m.n);
    for (int k = 1; k <= 6; ++k) {
        p = mul(p, x);
        num = combine(num, p, c[k]);
        den = combine(den, p, (k % 2 ? -1.0 : 1.0) * c[k]);
    }
    Cm r = mul(inverse(den), num);
    for (int k = 0; k < s; ++k) r = mul(r, r);
    return r;
}

// ------------------------------------------------------------------ elements

struct CavityBarrel {
    double x_m = 0, y_m = 0, radius_m = 0;
    bool shorted = false;                         // a barrel joining the two planes
};

// A barrel's segment in a pair of planes: between two of its nodes, its share of the
// pair's height.
struct CavityBranch {
    std::size_t upper = 0, lower = 0;
    std::size_t via = 0;  // index among the pair's barrels
    double weight = 0;
};

struct Element {
    enum Kind { Line, Barrel, Series, Shunt, Part, Cavity, Capacitances } kind = Line;
    std::vector<std::size_t> near, far;  // nodes (one per conductor)
    std::size_t cs = 0;                  // Line: cross-section
    std::vector<std::size_t> conductors; // Line: conductors of the cross-section
    double length = 0;                   // m
    double zv = 0, eps = 1;              // Barrel
    double r = 0, l = 0, c = 0;          // Series (in series), Shunt (r and c in parallel)
    const Touchstone* ts = nullptr;      // Part
    // Cavity: a pair of planes h apart, filled with the dielectrics of the given
    // thicknesses, and the barrels through it.
    double h = 0;
    std::vector<double> thickness_m;
    std::vector<DielectricModel> materials;
    std::vector<CavityBarrel> ports;
    std::shared_ptr<const kairoq::Cavity> cavity;  // the pair's outline; without it, infinite planes
    std::vector<CavityBranch> branches;  // Cavity: the barrels' segments in it
    // Capacitances: Maxwell matrix of the nodes `near` (F); Cavity: the branches' own
    // inductance from the pair's higher modes (H).
    std::vector<double> matrix;
};

// Admittance y between two nodes, either of which may be the reference.
void stamp_branch(std::vector<cplx>& Y, std::size_t nn, std::size_t a, std::size_t b, cplx y) {
    if (a != kGround) Y[a * nn + a] += y;
    if (b != kGround) Y[b * nn + b] += y;
    if (a != kGround && b != kGround) {
        Y[a * nn + b] -= y;
        Y[b * nn + a] -= y;
    }
}

// Admittance y from port q to port p, each a pair of nodes (signal, return).
void stamp_ports(std::vector<cplx>& Y, std::size_t nn, std::size_t sp, std::size_t rp, std::size_t sq, std::size_t rq, cplx y) {
    auto add = [&](std::size_t i, std::size_t j, cplx v) {
        if (i != kGround && j != kGround) Y[i * nn + j] += v;
    };
    add(sp, sq, y);
    add(sp, rq, -y);
    add(rp, sq, -y);
    add(rp, rq, y);
}

// Admittance blocks of a 2N-port line section from its ABCD matrix (currents into the
// section at both ends): Y11 = D B^-1, Y12 = C - D B^-1 A, Y21 = -B^-1, Y22 = B^-1 A.
void stamp_abcd(const Cm& T, const std::vector<std::size_t>& near, const std::vector<std::size_t>& far, std::vector<cplx>& Y,
                std::size_t nn) {
    const std::size_t n = near.size();
    Cm A(n), B(n), C(n), D(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            A(i, j) = T(i, j);
            B(i, j) = T(i, j + n);
            C(i, j) = T(i + n, j);
            D(i, j) = T(i + n, j + n);
        }
    const Cm Bi = inverse(B);
    const Cm DBi = mul(D, Bi);
    const Cm Y11 = DBi, Y12 = combine(C, mul(DBi, A), -1.0), Y22 = mul(Bi, A);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            Y[near[i] * nn + near[j]] += Y11(i, j);
            Y[near[i] * nn + far[j]] += Y12(i, j);
            Y[far[i] * nn + near[j]] += -Bi(i, j);
            Y[far[i] * nn + far[j]] += Y22(i, j);
        }
}

Cm section_abcd(const CsFrequency& f, std::size_t n, const std::vector<std::size_t>& conductors, double length) {
    const std::size_t k = conductors.size();
    const double w = 2 * kPi * f.frequency_hz;
    if (k == 1) {
        const std::size_t i = conductors[0] * n + conductors[0];
        const cplx z(f.r[i], w * f.l[i]), y(f.g[i], w * f.c[i]);
        const auto t = line_abcd(std::sqrt(z / y), std::sqrt(z * y), length);
        Cm T(2);
        T.a = {t[0], t[1], t[2], t[3]};
        return T;
    }
    Cm M(2 * k);
    for (std::size_t a = 0; a < k; ++a)
        for (std::size_t b = 0; b < k; ++b) {
            const std::size_t i = conductors[a] * n + conductors[b];
            M(a, b + k) = length * cplx(f.r[i], w * f.l[i]);
            M(a + k, b) = length * cplx(f.g[i], w * f.c[i]);
        }
    return expm(M);
}

// ------------------------------------------------------------------ channel model

struct Slot {
    std::size_t track, sample;  // sample: index into SiResult::samples
    std::size_t a, b;           // nodes at the start and end of the sample's length
    double weight;              // m
    double ux, uy;              // direction of the track at the sample
    bool matched = false;
};

class Channel {
public:
    Channel(const Board& b, const SiSettings& si, const ChannelSpec& spec, const NetworkSettings& ns, const std::string& base)
        : b_(b), si_(si), spec_(spec), ns_(ns), base_(base) {}
    ChannelResult run();

private:
    const Board& b_;
    const SiSettings& si_;
    const ChannelSpec& spec_;
    const NetworkSettings& ns_;
    const std::string& base_;
    std::vector<std::size_t> parent_;
    std::unordered_map<std::string, std::size_t> key_;
    std::vector<Element> elements_;
    std::vector<Touchstone> touchstones_;
    std::vector<std::string> warnings_;
    std::map<std::string, std::shared_ptr<const kairoq::Cavity>> cavity_cache_;

    std::size_t node(const std::string& key) {
        auto [it, fresh] = key_.emplace(key, parent_.size());
        if (fresh) parent_.push_back(parent_.size());
        return it->second;
    }
    std::size_t fresh_node() {
        parent_.push_back(parent_.size());
        return parent_.size() - 1;
    }
    std::size_t root(std::size_t v) {
        while (parent_[v] != v) v = parent_[v] = parent_[parent_[v]];
        return v;
    }
    void join(std::size_t a, std::size_t c) { parent_[root(a)] = root(c); }
    std::string end_key(int layer, Point64 p) const {
        return "t" + std::to_string(layer) + ":" + std::to_string(p.x / 1000) + ":" + std::to_string(p.y / 1000);
    }
    const Pad* find_pad(const std::string& ref, int layer = 0) const {
        const auto dot = ref.find('.');
        if (dot == std::string::npos) throw std::runtime_error("channel " + spec_.name + ": pad '" + ref + "' is not REF.PIN");
        const std::string comp = ref.substr(0, dot), pin = ref.substr(dot + 1);
        const Pad* best = nullptr;
        for (const auto& p : b_.pads)
            if (p.component == comp && p.pin == pin && (layer == 0 || p.layer == layer)) {
                if (!best || p.layer < best->layer) best = &p;
            }
        if (!best && layer == 0) throw std::runtime_error("channel " + spec_.name + ": no pad " + ref);
        return best;
    }
    std::string pad_key(const Pad& p) const { return "p" + p.component + "." + p.pin + ":" + std::to_string(p.layer); }
};

ChannelResult Channel::run() {
    const auto started = Clock::now();
    ChannelResult res;
    res.name = spec_.name;
    res.ports = spec_.ports;
    res.impedance_ohm = spec_.impedance_ohm;
    if (spec_.ports.empty()) throw std::runtime_error("channel " + spec_.name + " has no ports");
    // The channel's nets: those of its ports, and through its parts.
    std::set<int> nets;
    for (const auto& p : spec_.ports) nets.insert(find_pad(p)->net);
    for (bool grew = true; grew;) {
        grew = false;
        for (const auto& part : spec_.parts) {
            const int a = find_pad(part.from_pad)->net, c = find_pad(part.to_pad)->net;
            if (nets.count(a) && !nets.count(c)) grew = nets.insert(c).second;
            if (nets.count(c) && !nets.count(a)) grew = nets.insert(a).second;
        }
    }
    nets.erase(0);
    for (int id : nets) res.nets.push_back(b_.nets[static_cast<std::size_t>(id)]);
    // Their tracks as lines: the signal-integrity analysis of these nets.
    SiSettings si = si_;
    si.nets.clear();
    si.net_names = res.nets;
    si.targets.clear();
    const SiResult sr = analyse_si(b_, si);
    for (const auto& w : sr.warnings) warnings_.push_back(w);

    // Nodes: track ends, pads, via layers.
    std::map<std::size_t, std::vector<std::size_t>> samples_of;
    for (std::size_t k = 0; k < sr.samples.size(); ++k) samples_of[sr.samples[k].track].push_back(k);
    std::vector<std::size_t> tracks;
    for (std::size_t ti = 0; ti < b_.tracks.size(); ++ti)
        if (nets.count(b_.tracks[ti].net)) tracks.push_back(ti);
    for (auto ti : tracks) {
        node(end_key(b_.tracks[ti].layer, b_.tracks[ti].start));
        node(end_key(b_.tracks[ti].layer, b_.tracks[ti].end));
    }
    for (const auto& p : b_.pads) {
        if (!nets.count(p.net)) continue;
        const std::size_t pn = node(pad_key(p));
        const Rect64 pb = bounds(p.shape);
        for (auto ti : tracks) {
            const auto& tr = b_.tracks[ti];
            if (tr.layer != p.layer) continue;
            for (const Point64 q : {tr.start, tr.end})
                if (q.x >= pb.left && q.x <= pb.right && q.y >= pb.top && q.y <= pb.bottom && contains(p.shape, q))
                    join(node(end_key(tr.layer, q)), pn);
        }
    }
    // Vias: a node on each copper layer of the barrel, joined by the via model.
    const auto zs = stackup_z(b_);
    std::map<int, double> zmid;
    std::map<int, std::size_t> stack_of;  // copper layer -> stackup index
    for (std::size_t i = 0; i < b_.stackup.size(); ++i)
        if (b_.stackup[i].type == LayerType::Copper) {
            zmid[b_.stackup[i].copper_index] = 0.5 * (zs[i].top_nm + zs[i].bottom_nm) * 1e-9;
            stack_of[b_.stackup[i].copper_index] = i;
        }
    auto material = [&](std::size_t i) {
        const auto& l = b_.stackup[i];
        if (auto it = si_.layer_materials.find(l.name); it != si_.layer_materials.end()) return it->second;
        DielectricModel m = si_.default_dielectric;
        if (l.epsilon_r) m.epsilon_r = *l.epsilon_r;
        if (l.loss_tangent) m.loss_tangent = *l.loss_tangent;
        return m;
    };
    std::set<int> refs;
    for (const auto& name : si_.reference_nets) refs.insert(b_.find_net(name));
    if (refs.empty()) {
        const double area = area_nm2(b_.outline);
        for (const auto& cl : b_.copper)
            for (const auto& nc : cl.nets)
                if (nc.net != 0 && area_nm2(nc.polygons) >= si_.plane_fraction * area) refs.insert(nc.net);
    }
    const bool physical = ns_.via_model == "physical";
    // The reference copper around a point on a layer: none (0), a hole of the given
    // equal-area diameter (> 0), or solid copper (-1).
    auto plane_at = [&](int layer, Point64 at) {
        for (const auto& cl : b_.copper) {
            if (cl.index != layer) continue;
            for (const auto& nc : cl.nets) {
                if (!refs.count(nc.net)) continue;
                for (const auto& poly : nc.polygons) {
                    const Rect64 pb = bounds(poly);
                    if (at.x < pb.left || at.x > pb.right || at.y < pb.top || at.y > pb.bottom) continue;
                    if (!contains(Polygon{poly.outer, {}}, at)) continue;
                    for (const auto& h : poly.holes)
                        if (contains(Polygon{h, {}}, at)) return 2 * std::sqrt(std::abs(Clipper2Lib::Area(h)) / kPi);
                    return -1.0;
                }
            }
        }
        return 0.0;
    };
    // The pad around a barrel: the largest circle about its centre inside the net's copper.
    auto pad_radius = [&](int layer, int net, Point64 at) {
        double r = 0;
        for (const auto& cl : b_.copper) {
            if (cl.index != layer) continue;
            for (const auto& nc : cl.nets) {
                if (nc.net != net) continue;
                for (const auto& poly : nc.polygons) {
                    const Rect64 pb = bounds(poly);
                    if (at.x < pb.left || at.x > pb.right || at.y < pb.top || at.y > pb.bottom || !contains(poly, at)) continue;
                    double d = INFINITY;
                    std::vector<const Path64*> rings{&poly.outer};
                    for (const auto& h : poly.holes) rings.push_back(&h);
                    for (const auto* ring : rings)
                        for (std::size_t k = 0; k < ring->size(); ++k) {
                            const Point64 p = (*ring)[k], q = (*ring)[(k + 1) % ring->size()];
                            const double ex = double(q.x - p.x), ey = double(q.y - p.y), len2 = ex * ex + ey * ey;
                            const double t = len2 > 0 ? std::clamp((double(at.x - p.x) * ex + double(at.y - p.y) * ey) / len2, 0.0, 1.0) : 0.0;
                            d = std::min(d, std::hypot(double(p.x) + t * ex - double(at.x), double(p.y) + t * ey - double(at.y)));
                        }
                    r = std::max(r, d);
                }
            }
        }
        return r;
    };
    std::size_t fallback_antipads = 0;
    // Per pair of planes (upper, lower layer): the channel's barrels through it and their
    // segments in it.
    struct PairUse {
        std::vector<CavityBarrel> ports;
        std::vector<std::size_t> drill_of;
        std::vector<CavityBranch> branches;
        std::vector<std::pair<const ViaLadder*, std::size_t>> ladder_of;  // each branch's segment
    };
    std::map<std::pair<int, int>, PairUse> cavities;
    std::map<std::string, ViaLadder> ladder_cache;
    double er_max = 1.0;
    for (std::size_t i = 0; i < b_.stackup.size(); ++i)
        if (b_.stackup[i].type == LayerType::Dielectric) er_max = std::max(er_max, material(i).epsilon_r);
    // Ladder segments of a twentieth of the shortest wavelength.
    const double max_segment_nm = kC0 / (ns_.f_max_hz * std::sqrt(er_max)) / 20 * 1e9;
    const int ncopper = static_cast<int>(zmid.size());
    for (std::size_t d = 0; d < b_.drills.size(); ++d) {
        const auto& dr = b_.drills[d];
        if (!dr.plated || !nets.count(dr.net)) continue;
        const double r = 0.5 * dr.diameter_nm;
        std::vector<int> layers;
        for (int l = dr.from_layer; l <= dr.to_layer; ++l) {
            const std::size_t vn = node("v" + std::to_string(d) + ":" + std::to_string(l));
            layers.push_back(l);
            for (auto ti : tracks) {
                const auto& tr = b_.tracks[ti];
                if (tr.layer != l) continue;
                for (const Point64 q : {tr.start, tr.end})
                    if (std::hypot(double(q.x - dr.at.x), double(q.y - dr.at.y)) <= std::max(r, 0.5 * tr.width_nm) + 1000.0)
                        join(node(end_key(l, q)), vn);
            }
            for (const auto& p : b_.pads)
                if (p.net == dr.net && p.layer == l && contains(p.shape, dr.at)) join(node(pad_key(p)), vn);
        }
        auto vnode = [&](int l) { return node("v" + std::to_string(d) + ":" + std::to_string(l)); };
        if (!ns_.vias || layers.size() < 2) {
            for (std::size_t k = 1; k < layers.size(); ++k) join(vnode(layers[k - 1]), vnode(layers[k]));
            continue;
        }
        // Physics-based model: the reference planes around the barrel.
        if (physical) {
            ViaCut cut;
            cut.drill_nm = dr.diameter_nm;
            std::vector<int> planes;
            std::string key = std::to_string(std::llround(dr.diameter_nm)) + ":" + std::to_string(dr.from_layer) + "-" +
                              std::to_string(dr.to_layer);
            bool crossed_solid = false;
            for (int l = 1; l <= ncopper; ++l) {
                const std::size_t i = stack_of[l];
                ViaCut::Copper c;
                c.z_bottom_nm = zs[i].bottom_nm;
                c.z_top_nm = zs[i].top_nm;
                c.pad_nm = 2 * pad_radius(l, dr.net, dr.at);
                const double hole = plane_at(l, dr.at);
                const bool on_barrel = l >= dr.from_layer && l <= dr.to_layer;
                if (hole > 0 && hole > std::max(dr.diameter_nm, c.pad_nm) * 1.02) {
                    c.antipad_nm = hole;
                    planes.push_back(l);
                } else if (hole < 0 && !on_barrel) {
                    c.solid = true;
                    planes.push_back(l);
                } else if (hole != 0 && on_barrel) {
                    crossed_solid = true;
                }
                if (!on_barrel) c.pad_nm = 0;
                key += ":" + std::to_string(std::llround(c.pad_nm)) + "/" + std::to_string(std::llround(c.antipad_nm)) + (c.solid ? "s" : "");
                cut.copper.push_back(c);
            }
            if (!crossed_solid && !planes.empty()) {
                cut.from = static_cast<std::size_t>(dr.from_layer - 1);
                cut.to = static_cast<std::size_t>(dr.to_layer - 1);
                for (std::size_t i = 0; i < b_.stackup.size(); ++i)
                    if (b_.stackup[i].type == LayerType::Dielectric)
                        cut.dielectrics.push_back({zs[i].bottom_nm, zs[i].top_nm, material(i).epsilon_r});
                auto it = ladder_cache.find(key);
                if (it == ladder_cache.end()) it = ladder_cache.emplace(key, via_ladder(cut, max_segment_nm, 0.05)).first;
                const ViaLadder& lad = it->second;
                std::vector<std::size_t> ln;
                for (const auto& nd : lad.nodes) ln.push_back(nd.layer >= 0 ? vnode(nd.layer + 1) : fresh_node());
                {
                    Element e;
                    e.kind = Element::Capacitances;
                    e.near = ln;
                    e.matrix = lad.capacitance;
                    elements_.push_back(e);
                }
                const std::size_t ns = lad.segments.size();
                for (std::size_t k = 0; k < ns; ++k) {
                    const auto& sg = lad.segments[k];
                    if (sg.plane_upper < 0) {
                        Element e;
                        e.kind = Element::Series;
                        e.near = {ln[sg.upper]};
                        e.far = {ln[sg.lower]};
                        e.l = lad.inductance[k * ns + k];
                        elements_.push_back(e);
                        continue;
                    }
                    auto& g = cavities[{sg.plane_upper + 1, sg.plane_lower + 1}];
                    std::size_t vi = 0;
                    while (vi < g.ports.size() && g.drill_of[vi] != d) ++vi;
                    if (vi == g.ports.size()) {
                        g.ports.push_back({dr.at.x * 1e-9, dr.at.y * 1e-9, r * 1e-9, false});
                        g.drill_of.push_back(d);
                    }
                    g.branches.push_back({ln[sg.upper], ln[sg.lower], vi, sg.weight});
                    g.ladder_of.push_back({&lad, k});
                }
                ++res.vias;
                continue;
            }
        }
        // Coaxial model. The antipad: the hole in the reference copper around the barrel,
        // on the layers it passes without connecting; its diameter of equal area.
        double dsum = 0;
        int dn = 0;
        for (int l = 1; l <= ncopper; ++l) {
            const double hole = plane_at(l, dr.at);
            if (hole > 0) {
                dsum += hole;
                ++dn;
            }
        }
        double antipad = dn ? dsum / dn : 0.0;
        if (!(antipad > dr.diameter_nm * 1.05)) {
            antipad = 3 * dr.diameter_nm;
            ++fallback_antipads;
        }
        for (std::size_t k = 1; k < layers.size(); ++k) {
            // Permittivity: the dielectrics between the two layers.
            double eps_t = 0, t_sum = 0;
            for (std::size_t i = 0; i < b_.stackup.size(); ++i) {
                const auto& l = b_.stackup[i];
                const double zc = 0.5 * (zs[i].top_nm + zs[i].bottom_nm) * 1e-9;
                const double lo = std::min(zmid[layers[k - 1]], zmid[layers[k]]), hi = std::max(zmid[layers[k - 1]], zmid[layers[k]]);
                if (l.type != LayerType::Dielectric || zc < lo || zc > hi) continue;
                const double t = l.thickness_nm * 1e-9;
                eps_t += t * (l.epsilon_r ? *l.epsilon_r : si_.default_dielectric.epsilon_r);
                t_sum += t;
            }
            Element e;
            e.kind = Element::Barrel;
            e.near = {vnode(layers[k - 1])};
            e.far = {vnode(layers[k])};
            e.length = std::abs(zmid[layers[k - 1]] - zmid[layers[k]]);
            e.eps = t_sum > 0 ? eps_t / t_sum : si_.default_dielectric.epsilon_r;
            e.zv = 60.0 / std::sqrt(e.eps) * std::log(antipad / dr.diameter_nm);
            elements_.push_back(e);
        }
        ++res.vias;
    }
    if (fallback_antipads)
        warnings_.push_back(std::to_string(fallback_antipads) + " vias without an antipad in a reference plane; three drill diameters assumed");
    // The pairs of planes: their dielectrics, and the barrels of reference nets near the
    // channel's, which join the two planes.
    auto set_branches = [](Element& e, const PairUse& use) {
        e.branches = use.branches;
        const std::size_t nb = use.branches.size();
        e.matrix.assign(nb * nb, 0.0);
        for (std::size_t i = 0; i < nb; ++i)
            for (std::size_t j = 0; j < nb; ++j)
                if (use.ladder_of[i].first == use.ladder_of[j].first && use.branches[i].via == use.branches[j].via) {
                    const ViaLadder& lad = *use.ladder_of[i].first;
                    const std::size_t ns = lad.segments.size();
                    e.matrix[i * nb + j] = lad.inductance[use.ladder_of[i].second * ns + use.ladder_of[j].second];
                }
    };
    for (auto& [pair, use] : cavities) {
        const auto& ports = use.ports;
        Element e;
        e.kind = Element::Cavity;
        const double upper = zs[stack_of[pair.first]].bottom_nm, lower = zs[stack_of[pair.second]].top_nm;
        e.h = (upper - lower) * 1e-9;
        for (std::size_t i = 0; i < b_.stackup.size(); ++i)
            if (b_.stackup[i].type == LayerType::Dielectric && zs[i].top_nm <= upper + 1 && zs[i].bottom_nm >= lower - 1) {
                e.thickness_m.push_back((zs[i].top_nm - zs[i].bottom_nm) * 1e-9);
                e.materials.push_back(material(i));
            }
        // The finite pair: the overlap of the two planes, holes smaller than 2.5 mm across
        // (antipads, clearances) filled; the reference barrels joining them short it.
        {
            auto planes_of = [&](int layer) {
                PolygonSet u;
                for (const auto& cl : b_.copper)
                    if (cl.index == layer)
                        for (const auto& nc : cl.nets)
                            if (refs.count(nc.net)) u = unite(u, nc.polygons);
                return u;
            };
            PolygonSet region = intersection(planes_of(pair.first), planes_of(pair.second));
            for (auto& poly : region) {
                Paths64 kept;
                for (const auto& h : poly.holes)
                    if (std::abs(Clipper2Lib::Area(h)) > kPi * 1.25e6 * 1.25e6) kept.push_back(h);
                poly.holes = kept;
            }
            std::vector<CavityPort> cports;
            std::string key = std::to_string(std::llround(area_nm2(region))) + ":" + std::to_string(region.size());
            for (const auto& p : ports) {
                cports.push_back({p.x_m * 1e9, p.y_m * 1e9, p.radius_m * 1e9, false});
                key += ":" + std::to_string(std::llround(p.x_m * 1e6)) + "," + std::to_string(std::llround(p.y_m * 1e6));
            }
            for (const auto& g : b_.drills)
                if (g.plated && refs.count(g.net) && g.from_layer <= pair.first && g.to_layer >= pair.second && contains(region, g.at)) {
                    cports.push_back({double(g.at.x), double(g.at.y), 0.5 * g.diameter_nm, true});
                    key += ":g" + std::to_string(g.at.x / 1000) + "," + std::to_string(g.at.y / 1000);
                }
            // The finite pair where its modes to three times the sweep's wavenumber are few;
            // on large pairs, infinite planes with the nearby joining barrels shorted, which
            // then decide the impedance around the barrel.
            double er_pair = 1.0;
            for (const auto& m : e.materials) er_pair = std::max(er_pair, m.epsilon_r);
            const double k_eig = 3 * 2 * kPi * ns_.f_max_hz * std::sqrt(er_pair) / kC0;
            const double modes = area_nm2(region) * 1e-18 * k_eig * k_eig / (4 * kPi);
            bool inside = !region.empty() && modes <= 300;
            for (const auto& p : ports) inside = inside && contains(region, Point64(std::llround(p.x_m * 1e9), std::llround(p.y_m * 1e9)));
            if (inside) {
                double er_max = 1.0;
                for (const auto& m : e.materials) er_max = std::max(er_max, m.epsilon_r);
                auto it = cavity_cache_.find(key);
                if (it == cavity_cache_.end()) {
                    try {
                        it = cavity_cache_.emplace(key, std::make_shared<const kairoq::Cavity>(region, cports, 2 * kPi * ns_.f_max_hz * std::sqrt(er_max) / kC0)).first;
                    } catch (const std::exception& ex) {
                        warnings_.push_back(std::string("plane pair taken as infinite: ") + ex.what());
                    }
                }
                if (it != cavity_cache_.end()) {
                    e.cavity = it->second;
                    e.ports = ports;
                    set_branches(e, use);
                    elements_.push_back(e);
                    continue;
                }
            }
        }
        // The joining barrels nearest each of the channel's (eight within 5 mm or ten
        // heights): they decide the pair's impedance around it.
        const double reach = std::max(5e-3, 10 * e.h);
        std::vector<std::pair<double, CavityBarrel>> shorts;
        std::set<std::size_t> taken;
        for (const auto& p : ports) {
            std::vector<std::pair<double, std::size_t>> near;
            for (std::size_t gi = 0; gi < b_.drills.size(); ++gi) {
                const auto& g = b_.drills[gi];
                if (!g.plated || !refs.count(g.net) || g.from_layer > pair.first || g.to_layer < pair.second) continue;
                const double dist = std::hypot(g.at.x * 1e-9 - p.x_m, g.at.y * 1e-9 - p.y_m);
                if (dist <= reach) near.push_back({dist, gi});
            }
            std::sort(near.begin(), near.end());
            for (std::size_t k = 0; k < near.size() && k < 8; ++k)
                if (taken.insert(near[k].second).second) {
                    const auto& g = b_.drills[near[k].second];
                    shorts.push_back({near[k].first, {g.at.x * 1e-9, g.at.y * 1e-9, 0.5 * g.diameter_nm * 1e-9, true}});
                }
        }
        e.ports = ports;
        for (const auto& sp : shorts) e.ports.push_back(sp.second);
        set_branches(e, use);
        elements_.push_back(e);
    }

    // Parts and loads.
    for (const auto& part : spec_.parts) {
        const std::size_t a = node(pad_key(*find_pad(part.from_pad))), c = node(pad_key(*find_pad(part.to_pad)));
        Element e;
        e.near = {a};
        e.far = {c};
        if (!part.touchstone.empty()) {
            fs::path p(part.touchstone);
            if (p.is_relative()) p = fs::path(base_) / p;
            touchstones_.push_back(read_touchstone(p.string()));
            if (touchstones_.back().ports != 2) throw std::runtime_error("channel " + spec_.name + ": " + part.touchstone + " is not a two-port");
            e.kind = Element::Part;
        } else if (part.r_ohm == 0 && part.l_h == 0 && part.c_f == 0) {
            join(a, c);
            continue;
        } else {
            e.kind = Element::Series;
            e.r = part.r_ohm;
            e.l = part.l_h;
            e.c = part.c_f;
        }
        elements_.push_back(e);
    }
    for (std::size_t k = 0, t = 0; k < elements_.size(); ++k)
        if (elements_[k].kind == Element::Part) elements_[k].ts = &touchstones_[t++];
    for (const auto& load : spec_.loads) {
        Element e;
        e.kind = Element::Shunt;
        e.near = {node(pad_key(*find_pad(load.pad)))};
        e.r = load.r_ohm;
        e.c = load.c_f;
        elements_.push_back(e);
    }

    // Tracks: each sample a line of its length; pairs coupled where their samples match.
    std::vector<Slot> slots;
    std::map<int, std::vector<std::size_t>> slots_of_net;
    std::size_t borrowed = 0;
    for (auto ti : tracks) {
        const auto& tr = b_.tracks[ti];
        auto it = samples_of.find(ti);
        const std::size_t a = node(end_key(tr.layer, tr.start)), z = node(end_key(tr.layer, tr.end));
        std::vector<std::size_t> list;
        if (it != samples_of.end()) list = it->second;
        std::sort(list.begin(), list.end(), [&](std::size_t x, std::size_t y) { return sr.samples[x].along_nm < sr.samples[y].along_nm; });
        // Samples without a cross-section (no reference) borrow the nearest one of the track.
        std::vector<std::size_t> usable;
        for (auto k : list)
            if (sr.samples[k].z0 > 0) usable.push_back(k);
        if (usable.empty()) {
            // No line at all: a wire of 1 nH/mm.
            Element e;
            e.kind = Element::Series;
            e.near = {a};
            e.far = {z};
            e.l = 1e-6 * std::hypot(double(tr.end.x - tr.start.x), double(tr.end.y - tr.start.y)) * 1e-9;
            if (e.l > 0) elements_.push_back(e);
            else join(a, z);
            ++borrowed;
            continue;
        }
        std::size_t prev = a;
        for (std::size_t q = 0; q < list.size(); ++q) {
            const auto& s = sr.samples[list[q]];
            const std::size_t next = q + 1 == list.size() ? z : fresh_node();
            std::size_t src = list[q];
            if (!(s.z0 > 0)) {
                std::size_t best = usable.front();
                for (auto u : usable)
                    if (std::abs(sr.samples[u].along_nm - s.along_nm) < std::abs(sr.samples[best].along_nm - s.along_nm)) best = u;
                src = best;
                ++borrowed;
            }
            // The track direction at the sample.
            double ux = double(tr.end.x - tr.start.x), uy = double(tr.end.y - tr.start.y);
            if (tr.centre) {
                const double rx = s.x_nm - tr.centre->x, ry = s.y_nm - tr.centre->y;
                ux = tr.ccw ? -ry : ry;
                uy = tr.ccw ? rx : -rx;
            }
            const double un = std::hypot(ux, uy);
            slots.push_back(Slot{ti, src, prev, next, 0.0, ux / un, uy / un});
            slots_of_net[tr.net].push_back(slots.size() - 1);
            prev = next;
        }
    }
    // Weights: each sample stands for an equal share of its track (as in the analysis).
    {
        std::map<std::size_t, std::vector<std::size_t>> by_track;
        for (std::size_t k = 0; k < slots.size(); ++k) by_track[slots[k].track].push_back(k);
        for (const auto& [ti, list] : by_track) {
            const auto& tr = b_.tracks[ti];
            double len = std::hypot(double(tr.end.x - tr.start.x), double(tr.end.y - tr.start.y));
            if (tr.centre) {
                const double r = std::hypot(double(tr.start.x - tr.centre->x), double(tr.start.y - tr.centre->y));
                const double a0 = std::atan2(double(tr.start.y - tr.centre->y), double(tr.start.x - tr.centre->x));
                double sweep = std::atan2(double(tr.end.y - tr.centre->y), double(tr.end.x - tr.centre->x)) - a0;
                if (tr.ccw && sweep <= 0) sweep += 2 * kPi;
                if (!tr.ccw && sweep >= 0) sweep -= 2 * kPi;
                len = r * std::abs(sweep);
            }
            for (auto k : list) slots[k].weight = len * 1e-9 / double(list.size());
        }
    }
    if (borrowed) warnings_.push_back(std::to_string(borrowed) + " track sections without a reference take the line of their neighbours (or 1 nH/mm)");
    // Coupled sections: a sample with copper of other channel nets in its cut (a pair's
    // partner, an aggressor beside a victim) and those nets' samples nearest to that copper
    // become one section of several conductors, with the coupled parameters of the cut.
    for (auto& s : slots) {
        if (s.matched) continue;
        const auto& smp = sr.samples[s.sample];
        struct Member {
            std::size_t slot, conductor;
        };
        std::vector<Member> group;
        for (const auto& o : smp.others) {
            if (!nets.count(o.net) || o.net == smp.net) continue;
            auto it = slots_of_net.find(o.net);
            if (it == slots_of_net.end()) continue;
            std::size_t best = SIZE_MAX;
            double dbest = INFINITY;
            for (auto c : it->second) {
                const auto& os = sr.samples[slots[c].sample];
                if (slots[c].matched || os.layer != smp.layer) continue;
                bool taken = false;
                for (const auto& m : group) taken |= m.slot == c;
                if (taken) continue;
                const double d = std::hypot(os.x_nm - o.x_nm, os.y_nm - o.y_nm);
                if (d < dbest) {
                    dbest = d;
                    best = c;
                }
            }
            if (best == SIZE_MAX || dbest > 0.6e9 * std::max(s.weight, slots[best].weight) + 1e5) continue;
            group.push_back({best, o.conductor});
        }
        if (group.empty()) continue;
        s.matched = true;
        double lc = s.weight, lmax = s.weight;
        for (const auto& m : group) {
            slots[m.slot].matched = true;
            lc = std::min(lc, slots[m.slot].weight);
            lmax = std::max(lmax, slots[m.slot].weight);
        }
        const bool close = lmax - lc < 1e-5;
        if (close) lc = 0.5 * (lc + lmax);
        Element e;
        e.kind = Element::Line;
        e.cs = smp.cross_section;
        e.length = lc;
        // Each member from its start (in the direction of the sample's track) to a middle
        // node, the rest of its length a single line.
        auto member = [&](const Slot& m, std::size_t conductor, bool reversed) {
            const std::size_t mn = reversed ? m.b : m.a, mf = reversed ? m.a : m.b;
            std::size_t mid = mf;
            if (!close && m.weight > lc + 1e-9) {
                mid = fresh_node();
                Element r;
                r.kind = Element::Line;
                r.near = {mid};
                r.far = {mf};
                r.cs = sr.samples[m.sample].cross_section;
                r.conductors = {0};
                r.length = m.weight - lc;
                elements_.push_back(r);
            }
            e.near.push_back(mn);
            e.far.push_back(mid);
            e.conductors.push_back(conductor);
        };
        member(s, 0, false);
        for (const auto& m : group) {
            const Slot& o = slots[m.slot];
            member(o, m.conductor, s.ux * o.ux + s.uy * o.uy < 0);
        }
        elements_.push_back(e);
        ++res.coupled_sections;
    }
    for (const auto& s : slots) {
        if (s.matched) continue;
        Element e;
        e.kind = Element::Line;
        e.near = {s.a};
        e.far = {s.b};
        e.cs = sr.samples[s.sample].cross_section;
        e.conductors = {0};
        e.length = s.weight;
        elements_.push_back(e);
    }
    for (const auto& e : elements_) res.sections += e.kind == Element::Line;
    {
        double fh = 0;
        for (const auto& e : elements_)
            if (e.kind == Element::Line && sr.cross_sections[e.cs].result.model.microstrip.enabled)
                fh = std::max(fh, ns_.f_max_hz * sr.cross_sections[e.cs].result.model.microstrip.height_m * 1e-6);
        if (fh > kMicrostripFnMax) {
            char text[160];
            std::snprintf(text, sizeof text, "microstrip dispersion held at its formulas' range above %.1f GHz (f h above %.0f GHz mm)",
                          kMicrostripFnMax / fh * ns_.f_max_hz * 1e-9, kMicrostripFnMax);
            warnings_.push_back(text);
        }
    }
    // Compact node numbering.
    std::vector<std::size_t> port_nodes;
    for (const auto& p : spec_.ports) port_nodes.push_back(node(pad_key(*find_pad(p))));
    std::unordered_map<std::size_t, std::size_t> index;
    auto idx = [&](std::size_t v) {
        if (v == kGround) return kGround;
        const std::size_t r = root(v);
        auto [it, fresh] = index.emplace(r, index.size());
        return it->second;
    };
    for (auto& e : elements_) {
        for (auto& v : e.near) v = idx(v);
        for (auto& v : e.far) v = idx(v);
        for (auto& b : e.branches) {
            b.upper = idx(b.upper);
            b.lower = idx(b.lower);
        }
    }
    for (auto& v : port_nodes) v = idx(v);
    for (std::size_t i = 0; i < port_nodes.size(); ++i)
        for (std::size_t j = 0; j < i; ++j)
            if (port_nodes[i] == port_nodes[j])
                throw std::runtime_error("channel " + spec_.name + ": ports " + spec_.ports[j] + " and " + spec_.ports[i] + " are one node");
    const std::size_t nn = index.size();
    res.nodes = nn;

    // Frequencies, and the solution at each.
    const std::size_t np = port_nodes.size();
    const double z0 = spec_.impedance_ohm;
    for (std::size_t k = 1; k <= ns_.points; ++k) res.frequency_hz.push_back(ns_.f_max_hz * double(k) / double(ns_.points));
    res.s.resize(res.frequency_hz.size());
    std::vector<std::string> errors(res.frequency_hz.size());
    parallel_for(res.frequency_hz.size(), [&](std::size_t fi) {
        try {
            const double f = res.frequency_hz[fi], w = 2 * kPi * f;
            std::unordered_map<std::size_t, CsFrequency> lines;
            std::vector<cplx> Y(nn * nn, 0.0);
            for (const auto& e : elements_) {
                switch (e.kind) {
                case Element::Line: {
                    auto it = lines.find(e.cs);
                    if (it == lines.end()) it = lines.emplace(e.cs, sr.cross_sections[e.cs].result.model.at(f)).first;
                    const Cm T = section_abcd(it->second, sr.cross_sections[e.cs].result.model.n, e.conductors, e.length);
                    stamp_abcd(T, e.near, e.far, Y, nn);
                    break;
                }
                case Element::Barrel: {
                    const auto t = line_abcd(e.zv, cplx(0, w * std::sqrt(e.eps) / kC0), e.length);
                    Cm T(2);
                    T.a = {t[0], t[1], t[2], t[3]};
                    stamp_abcd(T, e.near, e.far, Y, nn);
                    break;
                }
                case Element::Series: {
                    const cplx z = e.r + cplx(0, w * e.l) + (e.c > 0 ? 1.0 / cplx(0, w * e.c) : cplx(0));
                    stamp_branch(Y, nn, e.near[0], e.far[0], 1.0 / z);
                    break;
                }
                case Element::Shunt:
                    Y[e.near[0] * nn + e.near[0]] += (e.r > 0 ? 1.0 / e.r : 0.0) + cplx(0, w * e.c);
                    break;
                case Element::Cavity: {
                    // The pair's fundamental mode between the barrels (finite outline, or
                    // infinite planes with the nearby joining barrels short-circuited),
                    // shared among each barrel's segments by their heights; the higher
                    // modes' inductance of each barrel's own segments.
                    cplx inv = 0.0;
                    double t = 0;
                    for (std::size_t k = 0; k < e.materials.size(); ++k) {
                        inv += e.thickness_m[k] / e.materials[k].at(f);
                        t += e.thickness_m[k];
                    }
                    const cplx eps = t > 0 ? t / inv : cplx(1.0);
                    const cplx kk = w / kC0 * std::sqrt(eps);
                    std::size_t nv = 0;
                    for (const auto& p : e.ports) nv += !p.shorted;
                    Cm Zv(nv);
                    if (e.cavity) {
                        Zv.a = e.cavity->impedance(w, e.h, eps);
                    } else {
                        const std::size_t m = e.ports.size();
                        Cm Z(m);
                        for (std::size_t i = 0; i < m; ++i)
                            for (std::size_t j = 0; j < m; ++j)
                                Z(i, j) = i == j ? radial_self(w, e.h, kk, e.ports[i].radius_m)
                                                 : radial_transfer(w, e.h, kk,
                                                                   std::max(std::hypot(e.ports[i].x_m - e.ports[j].x_m, e.ports[i].y_m - e.ports[j].y_m),
                                                                            e.ports[i].radius_m + e.ports[j].radius_m));
                        // Shorted barrels: zero voltage, so the open ports' impedance is the
                        // inverse of the open block of the admittance.
                        const Cm Yall = inverse(Z);
                        Cm Yo(nv);
                        for (std::size_t i = 0; i < nv; ++i)
                            for (std::size_t j = 0; j < nv; ++j) Yo(i, j) = Yall(i, j);
                        Zv = inverse(Yo);
                    }
                    const std::size_t nb = e.branches.size();
                    Cm Zb(nb);
                    for (std::size_t i = 0; i < nb; ++i)
                        for (std::size_t j = 0; j < nb; ++j)
                            Zb(i, j) = Zv(e.branches[i].via, e.branches[j].via) * e.branches[i].weight * e.branches[j].weight +
                                       cplx(0, w * e.matrix[i * nb + j]);
                    const Cm Yb = inverse(Zb);
                    for (std::size_t i = 0; i < nb; ++i)
                        for (std::size_t j = 0; j < nb; ++j)
                            stamp_ports(Y, nn, e.branches[i].upper, e.branches[i].lower, e.branches[j].upper, e.branches[j].lower, Yb(i, j));
                    break;
                }
                case Element::Capacitances: {
                    const std::size_t m = e.near.size();
                    for (std::size_t i = 0; i < m; ++i)
                        for (std::size_t j = 0; j < m; ++j) Y[e.near[i] * nn + e.near[j]] += cplx(0, w * e.matrix[i * m + j]);
                    break;
                }
                case Element::Part: {
                    // Y of the two-port from its S: Y = (I - S)(I + S)^-1 / R.
                    const auto s = e.ts->at(f);
                    Cm S(2), I = Cm::identity(2);
                    S.a = s;
                    const Cm yp = mul(combine(I, S, -1.0), inverse(combine(I, S, 1.0)));
                    const std::size_t nodes[2] = {e.near[0], e.far[0]};
                    for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 2; ++j) Y[nodes[i] * nn + nodes[j]] += yp(i, j) / e.ts->reference_ohm;
                    break;
                }
                }
            }
            // Port impedances: unit currents into each port node.
            std::vector<cplx> rhs(nn * np, 0.0);
            for (std::size_t p = 0; p < np; ++p) rhs[port_nodes[p] * np + p] = 1.0;
            const auto x = solve(Y, nn, rhs, np);
            Cm Zn(np);
            for (std::size_t i = 0; i < np; ++i)
                for (std::size_t j = 0; j < np; ++j) Zn(i, j) = x[port_nodes[i] * np + j] / z0;
            const Cm I = Cm::identity(np);
            const Cm S = mul(combine(Zn, I, -1.0), inverse(combine(Zn, I, 1.0)));
            res.s[fi] = S.a;
        } catch (const std::exception& ex) {
            errors[fi] = ex.what();
        }
    });
    for (const auto& e : errors)
        if (!e.empty()) throw std::runtime_error("channel " + spec_.name + ": " + e);

    // Mixed-mode parameters: differential and common modes of the port pairs.
    const std::size_t nd = spec_.differential.size();
    if (nd) {
        Cm M(np);  // rows: d1..dK, c1..cK, then the single ports left
        std::set<std::size_t> used;
        for (std::size_t k = 0; k < nd; ++k) {
            const auto [p, q] = spec_.differential[k];
            if (p >= np || q >= np || p == q) throw std::runtime_error("channel " + spec_.name + ": differential ports out of range");
            M(k, p) = 1 / std::sqrt(2.0);
            M(k, q) = -1 / std::sqrt(2.0);
            M(nd + k, p) = 1 / std::sqrt(2.0);
            M(nd + k, q) = 1 / std::sqrt(2.0);
            used.insert(p);
            used.insert(q);
        }
        std::size_t row = 2 * nd;
        for (std::size_t p = 0; p < np && row < np; ++p)
            if (!used.count(p)) M(row++, p) = 1.0;
        Cm Mt(np);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j) Mt(i, j) = M(j, i);
        for (const auto& s : res.s) {
            Cm S(np);
            S.a = s;
            const Cm mm = mul(mul(M, S), Mt);
            std::vector<cplx> out(4 * nd * nd);
            for (std::size_t i = 0; i < 2 * nd; ++i)
                for (std::size_t j = 0; j < 2 * nd; ++j) out[i * 2 * nd + j] = mm(i, j);
            res.mixed.push_back(std::move(out));
        }
    }

    // Time-domain reflectometry of each port (each differential port), from the reflection
    // on the linear frequency grid: a Gaussian step of the given rise time, delayed by four
    // standard deviations so that it starts from rest.
    double total = 0;
    for (const auto& n : sr.nets) total += n.delay_ps * 1e-12;
    const double duration = ns_.tdr_duration_s > 0 ? ns_.tdr_duration_s : std::max(1e-9, 2.0 * total);
    const double df = ns_.f_max_hz / double(ns_.points);
    // The step must be band-limited by the sweep: a rise time of at least 0.9 / f_max,
    // where its spectrum has fallen to 9 %; otherwise the truncation rings.
    const double rise = std::max(ns_.rise_time_s, 0.9 / ns_.f_max_hz);
    if (rise > ns_.rise_time_s * 1.0001) {
        char text[160];
        std::snprintf(text, sizeof text, "TDR rise time %.3g ps, the shortest the sweep to %.3g GHz resolves", rise * 1e12, ns_.f_max_hz * 1e-9);
        res.warnings.push_back(text);
    }
    const double sigma = rise / 2.563, t0 = 4 * sigma;
    // The sampled spectrum makes the response periodic in 1 / df: a later time folds back.
    if (duration + t0 > 1.0 / df) {
        char text[160];
        std::snprintf(text, sizeof text, "TDR window %.3g ns beyond the sweep's period %.3g ns: the response wraps; more points lengthen the period",
                      (duration + t0) * 1e9, 1e9 / df);
        res.warnings.push_back(text);
    }
    const double dt = std::min(rise / 10, 1.0 / (4 * ns_.f_max_hz));
    // Response to the step of one parameter (a reflection or a transmission) over time.
    const std::size_t nt = static_cast<std::size_t>(std::ceil((duration + t0) / dt));
    std::vector<double> times;
    for (std::size_t it = 0; it <= nt; ++it)
        if (double(it) * dt >= t0) times.push_back(double(it) * dt - t0);
    auto step_response = [&](auto param) {
        const std::size_t nf = res.frequency_hz.size();
        const double g0 = param(0).real();  // the lowest frequency's value, extended to DC
        std::vector<double> out;
        double step = 0;
        for (std::size_t it = 0; it <= nt; ++it) {
            const double t = double(it) * dt;
            double h = g0;
            for (std::size_t k = 0; k < nf; ++k) {
                const double f = res.frequency_hz[k];
                const double H = std::exp(-0.5 * std::pow(2 * kPi * f * sigma, 2));
                h += 2 * (param(k) * std::exp(cplx(0, 2 * kPi * f * (t - t0)))).real() * H;
            }
            step += h * df * dt;
            if (t >= t0) out.push_back(step);
        }
        return out;
    };
    auto tdr = [&](const std::string& name, auto gamma, double zref) {
        TdrTrace tr;
        tr.port = name;
        tr.time_s = times;
        for (double rho : step_response(gamma)) {
            rho = std::clamp(rho, -0.999999, 0.999999);
            tr.impedance_ohm.push_back(zref * (1 + rho) / (1 - rho));
        }
        res.tdr.push_back(std::move(tr));
    };
    if (nd) {
        for (std::size_t k = 0; k < nd; ++k)
            tdr("d" + std::to_string(k + 1), [&](std::size_t fi) { return res.mixed[fi][k * 2 * nd + k]; }, 2 * z0);
    } else {
        for (std::size_t p = 0; p < np; ++p) tdr(spec_.ports[p], [&](std::size_t fi) { return res.s[fi][p * np + p]; }, z0);
    }
    // Waveforms at every port for a step into each driven port, as fractions of the step:
    // the near- and far-end crosstalk of an aggressor, the received edge of a victim.
    for (std::size_t d : spec_.drive) {
        if (d >= np) throw std::runtime_error("channel " + spec_.name + ": driven port out of range");
        StepResponse sr2;
        sr2.driven = spec_.ports[d];
        sr2.time_s = times;
        for (std::size_t j = 0; j < np; ++j) {
            sr2.ports.push_back(spec_.ports[j]);
            sr2.values.push_back(step_response([&](std::size_t fi) { return res.s[fi][j * np + d]; }));
        }
        res.steps.push_back(std::move(sr2));
    }
    res.warnings.insert(res.warnings.begin(), warnings_.begin(), warnings_.end());
    res.seconds = since(started);
    return res;
}

}  // namespace

std::array<cplx, 4> line_abcd(cplx zc, cplx gamma, double length_m) {
    const cplx x = gamma * length_m;
    const cplx ch = std::cosh(x), sh = std::sinh(x);
    return {ch, zc * sh, sh / zc, ch};
}

std::array<cplx, 4> s_from_abcd(const std::array<cplx, 4>& t, double z0) {
    const cplx A = t[0], B = t[1], C = t[2], D = t[3];
    const cplx den = A + B / z0 + C * z0 + D;
    return {(A + B / z0 - C * z0 - D) / den, 2.0 * (A * D - B * C) / den, 2.0 / den, (-A + B / z0 - C * z0 + D) / den};
}

std::vector<ChannelResult> analyse_channels(const Board& board, const SiSettings& si, const std::vector<ChannelSpec>& channels,
                                            const NetworkSettings& settings, const std::string& base_dir) {
    if (!(settings.f_max_hz > 0) || settings.points < 2) throw std::runtime_error("network: f_max_hz must be positive, points at least 2");
    // A track that ends on the inside of another track of its net and layer (a stub routed
    // onto a line) meets it at a junction: the other track is split there, in a copy.
    Board split = board;
    bool any_split = false;
    // Ends of tracks on one layer within 1 um of each other are one point: arcs and lines
    // drawn to the same place may differ by rounding.
    {
        std::map<std::tuple<int, std::int64_t, std::int64_t>, std::vector<Point64>> cells;
        auto snap = [&](int layer, Point64& p) {
            const std::int64_t cx = p.x / 1000, cy = p.y / 1000;
            for (std::int64_t dx = -1; dx <= 1; ++dx)
                for (std::int64_t dy = -1; dy <= 1; ++dy) {
                    auto it = cells.find({layer, cx + dx, cy + dy});
                    if (it == cells.end()) continue;
                    for (const Point64& q : it->second)
                        if (std::abs(q.x - p.x) <= 1000 && std::abs(q.y - p.y) <= 1000) {
                            if (q != p) {
                                p = q;
                                any_split = true;
                            }
                            return;
                        }
                }
            cells[{layer, cx, cy}].push_back(p);
        };
        for (auto& t : split.tracks) {
            snap(t.layer, t.start);
            snap(t.layer, t.end);
        }
    }
    for (bool again = true; again;) {
        again = false;
        for (std::size_t i = 0; i < split.tracks.size() && !again; ++i) {
            for (const Point64 e : {split.tracks[i].start, split.tracks[i].end}) {
                for (std::size_t j = 0; j < split.tracks.size(); ++j) {
                    if (j == i) continue;
                    const Track& t = split.tracks[j];
                    if (t.centre || t.layer != split.tracks[i].layer || t.net != split.tracks[i].net) continue;
                    if (e == t.start || e == t.end) continue;
                    const double dx = double(t.end.x - t.start.x), dy = double(t.end.y - t.start.y), len2 = dx * dx + dy * dy;
                    if (len2 <= 0) continue;
                    const double u = (double(e.x - t.start.x) * dx + double(e.y - t.start.y) * dy) / len2;
                    if (u <= 1e-6 || u >= 1 - 1e-6) continue;
                    const double px = t.start.x + u * dx - e.x, py = t.start.y + u * dy - e.y;
                    if (px * px + py * py > 1e6) continue;  // within 1 um of the centre line
                    Track second = t;
                    second.start = e;
                    split.tracks[j].end = e;
                    split.tracks.push_back(second);
                    any_split = again = true;
                    break;
                }
                if (again) break;
            }
        }
    }
    const Board& base = any_split ? split : board;
    std::vector<ChannelResult> out;
    for (const auto& c : channels) {
        bool points = false;
        for (const auto& p : c.ports) points |= p.find('@') != std::string::npos;
        if (!points) {
            out.push_back(Channel(base, si, c, settings, base_dir).run());
            continue;
        }
        // Ports at points (NAME@X,Y[,LAYER], mm): a small pad named NAME.1 on the end of the
        // track nearest the point, on that layer (1 when not given), in a copy of the board.
        Board b = base;
        ChannelSpec spec = c;
        for (auto& p : spec.ports) {
            const auto at = p.find('@');
            if (at == std::string::npos) continue;
            const std::string name = p.substr(0, at);
            double x = 0, y = 0;
            int layer = 1;
            const std::string rest = p.substr(at + 1);
            char extra[32] = {0};
            const int n = std::sscanf(rest.c_str(), "%lf,%lf,%31s", &x, &y, extra);
            if (n < 2 || name.empty()) throw std::runtime_error("channel " + c.name + ": port '" + p + "' is not NAME@X,Y[,LAYER]");
            if (n == 3) layer = std::atoi(extra[0] == 'L' || extra[0] == 'l' ? extra + 1 : extra);
            const Point64 q(std::llround(x * 1e6), std::llround(y * 1e6));
            const Track* best = nullptr;
            Point64 end;
            double dmin = 2e6;  // within 2 mm
            for (const auto& t : b.tracks)
                if (t.layer == layer)
                    for (const Point64 e : {t.start, t.end}) {
                        const double d = std::hypot(double(e.x - q.x), double(e.y - q.y));
                        if (d < dmin) {
                            dmin = d;
                            best = &t;
                            end = e;
                        }
                    }
            if (!best) throw std::runtime_error("channel " + c.name + ": no track end on layer " + std::to_string(layer) + " within 2 mm of port " + name);
            Pad pad;
            pad.component = name;
            pad.pin = "1";
            pad.function = "SMDPad";
            pad.net = best->net;
            pad.layer = layer;
            pad.at = end;
            const std::int64_t h = 25000;  // half side, nm
            pad.shape = {Polygon{{Point64(end.x - h, end.y - h), Point64(end.x + h, end.y - h), Point64(end.x + h, end.y + h), Point64(end.x - h, end.y + h)}, {}}};
            b.pads.push_back(pad);
            p = name + ".1";
        }
        out.push_back(Channel(b, si, spec, settings, base_dir).run());
    }
    return out;
}

// ------------------------------------------------------------------ Touchstone

std::vector<cplx> Touchstone::at(double f) const {
    if (frequency_hz.empty()) throw std::runtime_error("Touchstone: no data");
    if (f <= frequency_hz.front()) return s.front();
    if (f >= frequency_hz.back()) return s.back();
    const auto it = std::upper_bound(frequency_hz.begin(), frequency_hz.end(), f);
    const std::size_t k = static_cast<std::size_t>(it - frequency_hz.begin());
    const double w = (f - frequency_hz[k - 1]) / (frequency_hz[k] - frequency_hz[k - 1]);
    std::vector<cplx> out(s[k].size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        const cplx a = s[k - 1][i], b = s[k][i];
        double pa = std::arg(a), pb = std::arg(b);
        while (pb - pa > kPi) pb -= 2 * kPi;
        while (pb - pa < -kPi) pb += 2 * kPi;
        out[i] = std::polar((1 - w) * std::abs(a) + w * std::abs(b), (1 - w) * pa + w * pb);
    }
    return out;
}

Touchstone parse_touchstone(const std::string& text, std::size_t ports) {
    Touchstone t;
    t.ports = ports;
    double unit = 1e9;  // GHz default
    std::string format = "MA";
    std::vector<double> numbers;
    std::istringstream in(text);
    std::string line;
    bool option = false;
    while (std::getline(in, line)) {
        if (const auto bang = line.find('!'); bang != std::string::npos) line = line.substr(0, bang);
        std::istringstream ls(line);
        std::string first;
        if (!(ls >> first)) continue;
        if (first == "#") {
            if (option) continue;  // only the first option line counts
            option = true;
            std::string w;
            while (ls >> w) {
                for (auto& ch : w) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (w == "HZ") unit = 1;
                else if (w == "KHZ") unit = 1e3;
                else if (w == "MHZ") unit = 1e6;
                else if (w == "GHZ") unit = 1e9;
                else if (w == "MA" || w == "DB" || w == "RI") format = w;
                else if (w == "R") ls >> t.reference_ohm;
                else if (w != "S") throw std::runtime_error("Touchstone: only S-parameters are read");
            }
            continue;
        }
        if (first[0] == '[') throw std::runtime_error("Touchstone version 2 keywords are not read");
        numbers.push_back(std::stod(first));
        double v;
        while (ls >> v) numbers.push_back(v);
    }
    const std::size_t per = 1 + 2 * ports * ports;
    if (numbers.size() % per != 0) throw std::runtime_error("Touchstone: data do not fill whole frequency rows");
    for (std::size_t r = 0; r < numbers.size(); r += per) {
        t.frequency_hz.push_back(numbers[r] * unit);
        std::vector<cplx> m(ports * ports);
        for (std::size_t k = 0; k < ports * ports; ++k) {
            const double a = numbers[r + 1 + 2 * k], b = numbers[r + 2 + 2 * k];
            cplx v;
            if (format == "RI") v = cplx(a, b);
            else if (format == "MA") v = std::polar(a, b * kPi / 180);
            else v = std::polar(std::pow(10.0, a / 20), b * kPi / 180);
            // Two-ports list S11 S21 S12 S22; others row by row.
            std::size_t i = k / ports, j = k % ports;
            if (ports == 2) std::swap(i, j);
            m[i * ports + j] = v;
        }
        t.s.push_back(std::move(m));
    }
    return t;
}

Touchstone read_touchstone(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot read " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string ext = fs::path(path).extension().string();
    std::size_t ports = 0;
    if (ext.size() >= 4 && (ext[1] == 's' || ext[1] == 'S') && (ext.back() == 'p' || ext.back() == 'P'))
        ports = static_cast<std::size_t>(std::stoul(ext.substr(2, ext.size() - 3)));
    if (!ports) throw std::runtime_error(path + ": the extension (.sNp) gives no port count");
    return parse_touchstone(ss.str(), ports);
}

std::string touchstone(const ChannelResult& r, bool mixed) {
    const auto& data = mixed ? r.mixed : r.s;
    if (data.empty()) throw std::runtime_error("channel " + r.name + ": no " + std::string(mixed ? "mixed-mode " : "") + "parameters");
    const std::size_t n = static_cast<std::size_t>(std::llround(std::sqrt(double(data.front().size()))));
    std::ostringstream o;
    o << "! kairoq channel " << r.name << (mixed ? ", mixed-mode: differential ports then common" : "") << "\n";
    o << "! ports:";
    for (const auto& p : r.ports) o << " " << p;
    o << "\n# Hz S RI R " << (mixed ? "100" : "50") << "\n";
    char buf[64];
    for (std::size_t fi = 0; fi < r.frequency_hz.size(); ++fi) {
        std::snprintf(buf, sizeof buf, "%.9g", r.frequency_hz[fi]);
        o << buf;
        const auto& m = data[fi];
        std::size_t count = 0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                const cplx v = n == 2 ? m[j * n + i] : m[i * n + j];  // two-ports: S11 S21 S12 S22
                std::snprintf(buf, sizeof buf, " %.12g %.12g", v.real(), v.imag());
                o << buf;
                if (n > 2 && ++count % 4 == 0 && !(i == n - 1 && j == n - 1)) o << "\n";
            }
        o << "\n";
    }
    return o.str();
}

std::string channels_json(const std::vector<ChannelResult>& rs) {
    using nlohmann::json;
    json out = json::array();
    auto mat = [](const std::vector<std::vector<cplx>>& v) {
        json a = json::array();
        for (const auto& m : v) {
            json row = json::array();
            for (const auto& c : m) row.push_back({c.real(), c.imag()});
            a.push_back(row);
        }
        return a;
    };
    for (const auto& r : rs) {
        json t = json::array();
        for (const auto& tr : r.tdr) t.push_back({{"port", tr.port}, {"time_s", tr.time_s}, {"impedance_ohm", tr.impedance_ohm}});
        json st = json::array();
        for (const auto& x : r.steps) st.push_back({{"driven", x.driven}, {"ports", x.ports}, {"time_s", x.time_s}, {"values", x.values}});
        out.push_back({{"name", r.name}, {"ports", r.ports}, {"nets", r.nets}, {"frequency_hz", r.frequency_hz}, {"s", mat(r.s)},
                       {"mixed", mat(r.mixed)}, {"tdr", t}, {"steps", st}, {"nodes", r.nodes}, {"sections", r.sections},
                       {"coupled_sections", r.coupled_sections}, {"vias", r.vias}, {"seconds", r.seconds}, {"warnings", r.warnings}});
    }
    return out.dump();
}

std::string channels_summary(const std::vector<ChannelResult>& rs) {
    std::ostringstream o;
    char line[256];
    for (const auto& r : rs) {
        std::snprintf(line, sizeof line, "channel %s: %zu ports, %zu nets, %zu lines (%zu coupled), %zu vias, %zu nodes, %.1f s\n",
                      r.name.c_str(), r.ports.size(), r.nets.size(), r.sections, r.coupled_sections, r.vias, r.nodes, r.seconds);
        o << line;
        const bool mm = !r.mixed.empty();
        const auto& data = mm ? r.mixed : r.s;
        const std::size_t n = static_cast<std::size_t>(std::llround(std::sqrt(double(data.front().size()))));
        const std::size_t half = mm ? n / 2 : 0;
        o << (mm ? "  f GHz     SDD11 dB  SDD21 dB  SCD21 dB\n" : "  f GHz     S11 dB    S21 dB\n");
        for (double f : {1e9, 2.5e9, 5e9, 10e9, 20e9}) {
            if (f > r.frequency_hz.back() * 1.0001) continue;
            std::size_t k = 0;
            for (std::size_t i = 0; i < r.frequency_hz.size(); ++i)
                if (std::abs(r.frequency_hz[i] - f) < std::abs(r.frequency_hz[k] - f)) k = i;
            const auto& m = data[k];
            auto db = [](cplx v) { return 20 * std::log10(std::max(std::abs(v), 1e-15)); };
            if (mm) {
                const double scd = half >= 2 ? db(m[(half + 1) * n + 0]) : -999.0;
                std::snprintf(line, sizeof line, "  %6.2f  %9.2f %9.2f %9.2f\n", r.frequency_hz[k] * 1e-9, db(m[0]), half >= 2 ? db(m[1 * n + 0]) : -999.0, scd);
            } else {
                std::snprintf(line, sizeof line, "  %6.2f  %9.2f %9.2f\n", r.frequency_hz[k] * 1e-9, db(m[0]), n >= 2 ? db(m[1 * n + 0]) : -999.0);
            }
            o << line;
        }
        for (const auto& t : r.tdr) {
            double lo = INFINITY, hi = -INFINITY;
            for (double z : t.impedance_ohm) {
                lo = std::min(lo, z);
                hi = std::max(hi, z);
            }
            std::snprintf(line, sizeof line, "  TDR %s: %.1f to %.1f ohm over %.2f ns\n", t.port.c_str(), lo, hi, t.time_s.empty() ? 0.0 : t.time_s.back() * 1e9);
            o << line;
        }
        for (const auto& x : r.steps)
            for (std::size_t j = 0; j < x.ports.size(); ++j) {
                if (x.ports[j] == x.driven) continue;
                double peak = 0, at = 0;
                for (std::size_t k = 0; k < x.values[j].size(); ++k)
                    if (std::abs(x.values[j][k]) > std::abs(peak)) {
                        peak = x.values[j][k];
                        at = x.time_s[k];
                    }
                std::snprintf(line, sizeof line, "  step into %s: at %s peak %+.2f %% at %.0f ps, final %+.2f %%\n", x.driven.c_str(),
                              x.ports[j].c_str(), 100 * peak, at * 1e12, 100 * x.values[j].back());
                o << line;
            }
        for (const auto& w : r.warnings) o << "  warning: " << w << "\n";
    }
    return o.str();
}

}  // namespace kairoq
