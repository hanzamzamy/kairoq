// Signal integrity of a board's tracks: cuts across each track, their cross-sections
// solved (cross_section.cpp), the profile, delay and loss of each net and the findings.
// Method: docs/solvers/si.md.
#include "kairoq/si.hpp"

#include "fem_util.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace kairoq {

namespace {

using namespace detail;

bool glob(const std::string& pattern, const std::string& text) {
    // '*' any sequence, '?' one character; case-sensitive, as net names are.
    std::size_t p = 0, t = 0, star = std::string::npos, mark = 0;
    while (t < text.size()) {
        if (p < pattern.size() && (pattern[p] == '?' || pattern[p] == text[t])) {
            ++p;
            ++t;
        } else if (p < pattern.size() && pattern[p] == '*') {
            star = p++;
            mark = t;
        } else if (star != std::string::npos) {
            p = star + 1;
            t = ++mark;
        } else {
            return false;
        }
    }
    while (p < pattern.size() && pattern[p] == '*') ++p;
    return p == pattern.size();
}

bool matches(const std::vector<std::string>& patterns, const std::string& name) {
    for (const auto& p : patterns)
        if (glob(p, name)) return true;
    return false;
}

// Two names differing only in their last P/N, p/n or +/- form a pair (P or + first).
std::optional<std::pair<std::string, std::string>> name_pair(const std::string& a, const std::string& b) {
    if (a.size() != b.size() || a.empty()) return std::nullopt;
    std::size_t diff = std::string::npos;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i]) {
            if (diff != std::string::npos) return std::nullopt;
            diff = i;
        }
    if (diff == std::string::npos) return std::nullopt;
    // The differing character must be the last letter or sign of the name, or the last
    // one before a trailing number.
    std::size_t last = a.size();
    while (last > 0 && std::isdigit(static_cast<unsigned char>(a[last - 1]))) --last;
    if (diff + 1 != last && diff + 1 != a.size()) return std::nullopt;
    const char x = a[diff], y = b[diff];
    for (const auto& [pc, nc] : {std::pair{'P', 'N'}, std::pair{'p', 'n'}, std::pair{'+', '-'}}) {
        if (x == pc && y == nc) return std::pair{a, b};
        if (x == nc && y == pc) return std::pair{b, a};
    }
    return std::nullopt;
}

// Arc sweep (signed, radians) and length (nm) of a track.
double track_sweep(const Track& tr) {
    if (!tr.centre) return 0.0;
    const double a0 = std::atan2(double(tr.start.y - tr.centre->y), double(tr.start.x - tr.centre->x));
    double sweep = std::atan2(double(tr.end.y - tr.centre->y), double(tr.end.x - tr.centre->x)) - a0;
    if (tr.ccw && sweep <= 0) sweep += 2 * kPi;
    if (!tr.ccw && sweep >= 0) sweep -= 2 * kPi;
    return sweep;
}
double track_length(const Track& tr) {
    if (tr.centre) return std::hypot(double(tr.start.x - tr.centre->x), double(tr.start.y - tr.centre->y)) * std::abs(track_sweep(tr));
    return std::hypot(double(tr.end.x - tr.start.x), double(tr.end.y - tr.start.y));
}

struct Interval {
    double s0, s1;  // nm along the cut
    int net;
};

// Copper of one layer for cuts: per net, its paths and bounding box.
struct LayerCopper {
    struct Net {
        int net;
        Paths64 paths;
        Rect64 box;
    };
    std::vector<Net> nets;
};

struct Task {
    int net = 0, layer = 0;
    std::size_t track = 0;
    double px = 0, py = 0;    // nm
    double ux = 0, uy = 0;    // unit tangent
    double width = 0;         // nm
    double weight = 0;        // track length this sample stands for (nm)
    double along = 0;         // position along the track from its start (nm)
};

struct Cut {
    std::string key;
    std::vector<CsConductor> conductors;
    std::vector<CsDielectric> dielectrics;
    std::vector<std::string> roles;
    std::vector<int> reference_layers;
    bool no_reference = false;
    bool partner = false;
    double partner_s = 0;  // middle of the partner's copper along the cut (nm)
    struct Other {
        int net;
        std::size_t index;  // among the signal conductors
        double s;           // middle along the cut (nm)
    };
    std::vector<Other> others;
    std::optional<CsModel::Microstrip> microstrip;  // signal on an outer layer over one plane
};

class Analysis {
public:
    Analysis(const Board& b, const SiSettings& s) : b_(b), s_(s) {}
    SiResult run();

private:
    const Board& b_;
    const SiSettings& s_;
    std::set<int> refs_;
    std::map<int, int> partner_;
    std::vector<LayerCopper> copper_;   // index layer - 1
    std::vector<LayerZ> z_;             // per stackup layer
    std::vector<std::size_t> stack_of_copper_;  // copper index -> stackup index
    std::vector<std::string> warnings_;

    std::vector<Interval> intersect(int layer, double px, double py, double nx, double ny, double h) const;
    std::optional<Cut> make_cut(const Task& t) const;
    std::vector<CsDielectric> dielectrics(int top_layer, int bottom_layer, const std::vector<std::array<double, 4>>& top_rects,
                                          const std::vector<std::array<double, 4>>& bottom_rects,
                                          const std::vector<std::pair<double, double>>& board) const;
    Paths64 outline_;
    DielectricModel material_of(std::size_t stack_index) const;
};

DielectricModel Analysis::material_of(std::size_t i) const {
    const auto& l = b_.stackup[i];
    if (auto it = s_.layer_materials.find(l.name); it != s_.layer_materials.end()) return it->second;
    DielectricModel m = s_.default_dielectric;
    if (l.epsilon_r) m.epsilon_r = *l.epsilon_r;
    if (l.loss_tangent) m.loss_tangent = *l.loss_tangent;
    return m;
}

std::vector<Interval> Analysis::intersect(int layer, double px, double py, double nx, double ny, double h) const {
    std::vector<Interval> out;
    const Point64 a(std::llround(px - h * nx), std::llround(py - h * ny)), b(std::llround(px + h * nx), std::llround(py + h * ny));
    const Rect64 cb(std::min(a.x, b.x), std::min(a.y, b.y), std::max(a.x, b.x), std::max(a.y, b.y));
    for (const auto& n : copper_[static_cast<std::size_t>(layer - 1)].nets) {
        if (n.box.right < cb.left || n.box.left > cb.right || n.box.bottom < cb.top || n.box.top > cb.bottom) continue;
        Clipper2Lib::Clipper64 c;
        c.AddOpenSubject({Path64{a, b}});
        c.AddClip(n.paths);
        Paths64 closed, open;
        c.Execute(Clipper2Lib::ClipType::Intersection, Clipper2Lib::FillRule::NonZero, closed, open);
        for (const auto& piece : open) {
            if (piece.size() < 2) continue;
            double s0 = INFINITY, s1 = -INFINITY;
            for (const auto& q : piece) {
                const double s = (q.x - px) * nx + (q.y - py) * ny;
                s0 = std::min(s0, s);
                s1 = std::max(s1, s);
            }
            if (s1 - s0 > 1.0) out.push_back({s0, s1, n.net});
        }
    }
    std::sort(out.begin(), out.end(), [](const Interval& x, const Interval& y) { return x.s0 < y.s0; });
    return out;
}

std::vector<CsDielectric> Analysis::dielectrics(int top_layer, int bottom_layer, const std::vector<std::array<double, 4>>& top_rects,
                                                const std::vector<std::array<double, 4>>& bottom_rects,
                                                const std::vector<std::pair<double, double>>& board) const {
    // Stack layers between the copper layers of the cut (with the masks of the outer
    // layers): the dielectric layers with their material, inner copper layers filled
    // with the material next to them, conformal solder mask on the outer layers.
    std::vector<CsDielectric> out;
    const int n = static_cast<int>(b_.copper.size());
    const std::size_t i0 = top_layer == 1 ? 0 : stack_of_copper_[static_cast<std::size_t>(top_layer)];
    const std::size_t i1 = bottom_layer == n ? b_.stackup.size() - 1 : stack_of_copper_[static_cast<std::size_t>(bottom_layer)];
    // Across the board only: the cut's spans inside the outline.
    auto rect = [&](double z0, double z1) {
        Paths64 parts;
        for (const auto& [a, c] : board)
            parts.push_back({Point64(std::llround(a), std::llround(z0)), Point64(std::llround(c), std::llround(z0)),
                             Point64(std::llround(c), std::llround(z1)), Point64(std::llround(a), std::llround(z1))});
        return union_polygons(parts);
    };
    auto neighbour_material = [&](std::size_t i, bool down) -> std::optional<DielectricModel> {
        for (std::size_t k = i;;) {
            if (down ? k + 1 >= b_.stackup.size() : k == 0) return std::nullopt;
            k = down ? k + 1 : k - 1;
            if (b_.stackup[k].type == LayerType::Dielectric) return material_of(k);
        }
    };
    std::vector<CsDielectric> masks;
    for (std::size_t i = i0; i <= i1; ++i) {
        const auto& l = b_.stackup[i];
        const double z0 = z_[i].bottom_nm, z1 = z_[i].top_nm;
        if (!(z1 > z0)) continue;
        if (l.type == LayerType::Dielectric) {
            out.push_back({rect(z0, z1), material_of(i)});
        } else if (l.type == LayerType::Copper && l.copper_index > 1 && l.copper_index < n) {
            const bool lower = l.copper_index > n / 2;
            auto m = neighbour_material(i, !lower);
            if (!m) m = neighbour_material(i, lower);
            if (m) out.push_back({rect(z0, z1), *m});
        } else if (l.type == LayerType::SolderMask) {
            // Conformal: over the laminate and over the outer conductors.
            const bool top = i < stack_of_copper_[1];
            const std::size_t cu = stack_of_copper_[static_cast<std::size_t>(top ? 1 : n)];
            const double t = z1 - z0;
            const double base = top ? z_[cu].bottom_nm : z_[cu].top_nm;
            Paths64 parts = to_paths(rect(top ? base : base - t, top ? base + t : base));
            for (const auto& r : top ? top_rects : bottom_rects) {
                const double x0 = r[0] - t, x1 = r[1] + t, y0 = top ? r[2] : r[2] - t, y1 = top ? r[3] + t : r[3];
                parts.push_back({Point64(std::llround(x0), std::llround(y0)), Point64(std::llround(x1), std::llround(y0)),
                                 Point64(std::llround(x1), std::llround(y1)), Point64(std::llround(x0), std::llround(y1))});
            }
            DielectricModel m;
            m.epsilon_r = l.epsilon_r ? *l.epsilon_r : s_.mask_epsilon_r;
            m.loss_tangent = l.loss_tangent ? *l.loss_tangent : s_.mask_loss_tangent;
            m.frequency_hz = s_.default_dielectric.frequency_hz;
            // Clipped to the board, over the mask's full height.
            const PolygonSet mask = union_polygons(parts);
            const Rect64 mb = bounds(mask);
            masks.push_back({intersection(mask, rect(double(mb.top), double(mb.bottom))), m});
        }
    }
    out.insert(out.end(), masks.begin(), masks.end());
    return out;
}

std::optional<Cut> Analysis::make_cut(const Task& t) const {
    const int n = static_cast<int>(b_.copper.size());
    const double nx = -t.uy, ny = t.ux, h = s_.lateral_nm;
    std::vector<std::vector<Interval>> iv(static_cast<std::size_t>(n + 1));
    for (int l = 1; l <= n; ++l) iv[static_cast<std::size_t>(l)] = intersect(l, t.px, t.py, nx, ny, h);
    // The track's own copper.
    std::optional<Interval> self;
    for (const auto& i : iv[static_cast<std::size_t>(t.layer)])
        if (i.net == t.net && i.s0 <= 1000.0 && i.s1 >= -1000.0) self = i;
    if (!self) return std::nullopt;
    const double w = self->s1 - self->s0;
    auto is_ref = [&](const Interval& i) {
        return refs_.count(i.net) || (i.net != t.net && i.s1 - i.s0 >= std::max(s_.large_copper_nm, 5 * w));
    };
    // Shields: the nearest layers above and below with reference copper across the track
    // and a track width to either side.
    auto covers = [&](int l) {
        double reach = self->s0 - w;
        for (const auto& i : iv[static_cast<std::size_t>(l)])
            if (is_ref(i) && i.s0 <= reach + 1.0 && i.s1 > reach) reach = std::max(reach, i.s1);
        return reach >= self->s1 + w;
    };
    int top = 1, bottom = n;
    for (int l = t.layer - 1; l >= 1; --l)
        if (covers(l)) {
            top = l;
            break;
        }
    for (int l = t.layer + 1; l <= n; ++l)
        if (covers(l)) {
            bottom = l;
            break;
        }
    const auto pit = partner_.find(t.net);
    std::optional<Interval> partner;
    double partner_gap = INFINITY;
    auto gap = [&](const Interval& i) { return std::max({0.0, i.s0 - self->s1, self->s0 - i.s1}); };
    if (pit != partner_.end())
        for (int l = top; l <= bottom; ++l)
            for (const auto& i : iv[static_cast<std::size_t>(l)])
                if (i.net == pit->second && gap(i) <= s_.neighbour_nm && gap(i) < partner_gap && l == t.layer) {
                    partner = i;
                    partner_gap = gap(i);
                }
    // The cut reaches twenty times the distance to the farther shield: the field of a
    // microstrip falls as the inverse square of the distance, and the impedance of the
    // truncated cut is then within about 0.1 % of an unbounded one. Other nets' copper on
    // the track's layer within `neighbour_heights` times that distance (at most
    // `neighbour_nm`) is a quiet line of its own; other copper is held at the reference
    // potential, which leaves the track's capacitance unchanged and its inductance nearly so.
    const std::size_t st = stack_of_copper_[static_cast<std::size_t>(t.layer)];
    double dv = 0;
    for (int l : {top, bottom}) {
        if (l == t.layer) continue;
        const std::size_t sl = stack_of_copper_[static_cast<std::size_t>(l)];
        dv = std::max(dv, std::max(z_[st].bottom_nm - z_[sl].top_nm, z_[sl].bottom_nm - z_[st].top_nm));
    }
    if (!(dv > 0)) dv = std::max(w, 1.0e5);
    // Beyond the track and its partner, whole.
    double extent = std::max(std::abs(self->s0), std::abs(self->s1));
    if (partner) extent = std::max({extent, std::abs(partner->s0), std::abs(partner->s1)});
    const double reach = std::min(h, std::max(20 * dv, 3 * w) + extent);
    const double near = std::min(s_.neighbour_nm, s_.neighbour_heights * dv);
    Cut c;
    std::vector<std::array<double, 4>> top_rects, bottom_rects;
    std::ostringstream key;
    key << top << '/' << bottom << '/' << t.layer;
    std::size_t signals = 0;
    auto add = [&](const Interval& i, int l, const std::string& role, bool reference) {
        const std::size_t si = stack_of_copper_[static_cast<std::size_t>(l)];
        const double z0 = z_[si].bottom_nm, z1 = z_[si].top_nm;
        const double s0 = std::max(i.s0, -reach), s1 = std::min(i.s1, reach);
        if (!(s1 - s0 > 1000.0)) return;  // clipped to nothing
        // Etched conductors are narrower on the face away from the laminate they were
        // etched on: outer layers outward, inner layers their top unless the project names
        // the face.
        const double e = std::min(s_.etch_nm, 0.45 * (s1 - s0));
        bool narrow_bottom = l == n && n > 1;
        for (const std::string& name : {b_.copper[static_cast<std::size_t>(l - 1)].name, b_.stackup[si].name})
            if (const auto it = s_.etch_bottom.find(name); it != s_.etch_bottom.end()) narrow_bottom = it->second;
        Path64 p;
        if (narrow_bottom)
            p = {Point64(std::llround(s0 + e), std::llround(z0)), Point64(std::llround(s1 - e), std::llround(z0)),
                 Point64(std::llround(s1), std::llround(z1)), Point64(std::llround(s0), std::llround(z1))};
        else
            p = {Point64(std::llround(s0), std::llround(z0)), Point64(std::llround(s1), std::llround(z0)),
                 Point64(std::llround(s1 - e), std::llround(z1)), Point64(std::llround(s0 + e), std::llround(z1))};
        CsConductor k;
        k.name = role;
        k.shape = union_polygons({p});
        k.reference = reference;
        k.metal = s_.copper;
        if (b_.stackup[si].conductivity_s_per_m) k.metal.conductivity_s_per_m = *b_.stackup[si].conductivity_s_per_m;
        c.conductors.push_back(std::move(k));
        c.roles.push_back(role);
        if (!reference) {
            if (role != "self") c.others.push_back({i.net, signals, 0.5 * (s0 + s1)});
            ++signals;
        }
        if (l == 1) top_rects.push_back({s0, s1, z0, z1});
        if (l == n) bottom_rects.push_back({s0, s1, z0, z1});
        // Track and partner to 1 um, other copper to 5 um.
        const double q = role == "self" || role == "partner" ? 1000.0 : 5000.0;
        key << '|' << role.substr(0, 1) << l << ':' << std::llround(s0 / q) << ',' << std::llround(s1 / q);
    };
    add(*self, t.layer, "self", false);
    if (partner) {
        add(*partner, t.layer, "partner", false);
        c.partner = true;
        c.partner_s = 0.5 * (partner->s0 + partner->s1);
    }
    bool any_reference = false;
    // Reference copper of a layer, its pieces joined across gaps that are at most a tenth
    // of their distance to the track (a gap g at distance D changes the track's fields by
    // about (g/D)^2): fan-outs and pierced planes far from the track would otherwise give
    // the cut a hundred conductors to mesh.
    auto distance_to_track = [&](double sm, int l) {
        double lateral = std::max({0.0, self->s0 - sm, sm - self->s1});
        if (partner) lateral = std::min(lateral, std::max({0.0, partner->s0 - sm, sm - partner->s1}));
        const std::size_t sl = stack_of_copper_[static_cast<std::size_t>(l)];
        const double vertical = l == t.layer ? 0.0 : std::max(z_[st].bottom_nm - z_[sl].top_nm, z_[sl].bottom_nm - z_[st].top_nm);
        return std::hypot(lateral, std::max(vertical, 0.0));
    };
    for (int l = top; l <= bottom; ++l) {
        std::vector<Interval> refs;
        for (const auto& i : iv[static_cast<std::size_t>(l)]) {
            if (l == t.layer && i.s0 == self->s0 && i.s1 == self->s1 && i.net == self->net) continue;
            if (partner && l == t.layer && i.s0 == partner->s0 && i.s1 == partner->s1) continue;
            if (i.s1 < -reach || i.s0 > reach) continue;
            if (is_ref(i)) {
                refs.push_back(i);
                any_reference = true;
            } else if (l == t.layer && gap(i) <= near) {
                add(i, l, "neighbour " + b_.nets[static_cast<std::size_t>(i.net)], false);
            } else {
                refs.push_back(i);
            }
        }
        std::sort(refs.begin(), refs.end(), [](const Interval& a, const Interval& b) { return a.s0 < b.s0; });
        std::vector<Interval> merged;
        for (const auto& i : refs) {
            if (!merged.empty()) {
                Interval& m = merged.back();
                const double g = i.s0 - m.s1;
                if (g <= 0 || g <= s_.reference_merge * distance_to_track(0.5 * (m.s1 + i.s0), l)) {
                    m.s1 = std::max(m.s1, i.s1);
                    continue;
                }
            }
            merged.push_back(i);
        }
        for (const auto& i : merged) add(i, l, "reference", true);
    }
    // The reference layers are the shields: planes across the track above and below.
    for (int l : {top, bottom})
        if (l != t.layer && covers(l)) c.reference_layers.push_back(l);
    std::sort(c.reference_layers.begin(), c.reference_layers.end());
    c.reference_layers.erase(std::unique(c.reference_layers.begin(), c.reference_layers.end()), c.reference_layers.end());
    c.no_reference = !any_reference;
    // A microstrip: the track on an outer layer with one reference layer, on the side of
    // the board; its substrate the dielectrics between, in series.
    if ((t.layer == 1 || t.layer == n) && c.reference_layers.size() == 1) {
        const int rl = c.reference_layers.front();
        const std::size_t sr = stack_of_copper_[static_cast<std::size_t>(rl)];
        double thick = 0, inv = 0;
        CsModel::Microstrip ms;
        for (std::size_t i = std::min(st, sr) + 1; i < std::max(st, sr); ++i) {
            if (b_.stackup[i].type != LayerType::Dielectric) continue;
            const DielectricModel m = material_of(i);
            const double ti = z_[i].top_nm - z_[i].bottom_nm;
            if (ti > thick) ms.substrate = m;
            thick = std::max(thick, ti);
            inv += ti / m.epsilon_r;
        }
        double total = 0;
        for (std::size_t i = std::min(st, sr) + 1; i < std::max(st, sr); ++i)
            if (b_.stackup[i].type == LayerType::Dielectric) total += z_[i].top_nm - z_[i].bottom_nm;
        if (total > 0 && inv > 0) {
            ms.substrate.epsilon_r = total / inv;
            ms.enabled = true;
            ms.width_m = (self->s1 - self->s0) * 1e-9;
            ms.height_m = dv * 1e-9;
            c.microstrip = ms;
        }
    }
    // The board along the cut, within its reach.
    std::vector<std::pair<double, double>> board;
    {
        const Point64 a(std::llround(t.px - h * nx), std::llround(t.py - h * ny)), b2(std::llround(t.px + h * nx), std::llround(t.py + h * ny));
        Clipper2Lib::Clipper64 cl;
        cl.AddOpenSubject({Path64{a, b2}});
        cl.AddClip(outline_);
        Paths64 closed, open;
        cl.Execute(Clipper2Lib::ClipType::Intersection, Clipper2Lib::FillRule::NonZero, closed, open);
        for (const auto& piece : open) {
            double s0 = INFINITY, s1 = -INFINITY;
            for (const auto& q : piece) {
                const double sv = (q.x - t.px) * nx + (q.y - t.py) * ny;
                s0 = std::min(s0, sv);
                s1 = std::max(s1, sv);
            }
            s0 = std::max(s0, -reach);
            s1 = std::min(s1, reach);
            if (s1 - s0 > 1000.0) board.emplace_back(s0, s1);
        }
        if (board.empty()) board.emplace_back(-reach, reach);
    }
    c.dielectrics = dielectrics(top, bottom, top_rects, bottom_rects, board);
    for (const auto& [a, b2] : board) key << "|b" << std::llround(a / 5000.0) << ',' << std::llround(b2 / 5000.0);
    c.key = key.str();
    return c;
}

SiResult Analysis::run() {
    const auto t0 = Clock::now();
    SiResult res;
    const int n = static_cast<int>(b_.copper.size());
    if (n == 0) throw std::runtime_error("si: the board has no copper layers");
    // Reference nets.
    for (const auto& name : s_.reference_nets) {
        const int id = b_.find_net(name);
        if (id < 0) throw std::runtime_error("si: no net '" + name + "'");
        refs_.insert(id);
    }
    if (refs_.empty()) {
        const double board_area = area_nm2(b_.outline);
        for (const auto& cl : b_.copper)
            for (const auto& nc : cl.nets)
                if (nc.net != 0 && area_nm2(nc.polygons) >= s_.plane_fraction * board_area) refs_.insert(nc.net);
    }
    if (refs_.empty()) warnings_.push_back("no reference net found; only large copper of any net is a reference");
    // Pairs.
    auto add_pair = [&](const std::string& p, const std::string& q) {
        const int a = b_.find_net(p), c = b_.find_net(q);
        if (a < 0 || c < 0) throw std::runtime_error("si: pair " + p + " / " + q + ": no such net");
        partner_[a] = c;
        partner_[c] = a;
    };
    for (const auto& p : s_.pairs) add_pair(p.p, p.n);
    std::set<int> with_tracks;
    for (const auto& t : b_.tracks)
        if (t.net > 0) with_tracks.insert(t.net);
    if (s_.pair_by_name) {
        std::map<std::string, int> by_name;
        for (int id : with_tracks) by_name[b_.nets[static_cast<std::size_t>(id)]] = id;
        for (const auto& [name, id] : by_name) {
            if (partner_.count(id)) continue;
            // Candidates: the same name with the last P/N or +/- swapped.
            for (std::size_t i = name.size(); i-- > 0;) {
                const char ch = name[i];
                const char sw = ch == 'P' ? 'N' : ch == 'N' ? 'P' : ch == 'p' ? 'n' : ch == 'n' ? 'p' : ch == '+' ? '-' : ch == '-' ? '+' : 0;
                if (!sw) {
                    if (std::isdigit(static_cast<unsigned char>(ch))) continue;
                    break;
                }
                std::string other = name;
                other[i] = sw;
                auto it = by_name.find(other);
                if (it != by_name.end() && !partner_.count(it->second) && name_pair(name, other)) {
                    partner_[id] = it->second;
                    partner_[it->second] = id;
                }
                break;
            }
        }
    }
    // Nets to analyse.
    std::set<int> nets;
    for (int id : with_tracks) {
        if (refs_.count(id)) continue;
        const auto& name = b_.nets[static_cast<std::size_t>(id)];
        const bool named = std::find(s_.net_names.begin(), s_.net_names.end(), name) != s_.net_names.end();
        if ((s_.nets.empty() && s_.net_names.empty()) || named || matches(s_.nets, name)) nets.insert(id);
    }
    // Stack positions.
    z_ = stackup_z(b_);
    stack_of_copper_.assign(static_cast<std::size_t>(n + 1), 0);
    for (std::size_t i = 0; i < b_.stackup.size(); ++i)
        if (b_.stackup[i].type == LayerType::Copper) stack_of_copper_[static_cast<std::size_t>(b_.stackup[i].copper_index)] = i;
    for (std::size_t i = 0; i < b_.stackup.size(); ++i) {
        const auto& l = b_.stackup[i];
        if (l.type != LayerType::Dielectric || s_.layer_materials.count(l.name)) continue;
        char text[160];
        if (!l.epsilon_r) {
            std::snprintf(text, sizeof text, "no permittivity given for %s; %.3g assumed", l.name.c_str(), s_.default_dielectric.epsilon_r);
            warnings_.push_back(text);
        }
        if (!l.loss_tangent) {
            std::snprintf(text, sizeof text, "no loss tangent given for %s; %.3g assumed", l.name.c_str(), s_.default_dielectric.loss_tangent);
            warnings_.push_back(text);
        }
    }
    outline_ = to_paths(b_.outline);
    copper_.resize(static_cast<std::size_t>(n));
    for (const auto& cl : b_.copper)
        for (const auto& nc : cl.nets) copper_[static_cast<std::size_t>(cl.index - 1)].nets.push_back({nc.net, to_paths(nc.polygons), bounds(nc.polygons)});

    // Samples along the tracks of the nets.
    std::vector<Task> tasks;
    for (std::size_t ti = 0; ti < b_.tracks.size(); ++ti) {
        const auto& tr = b_.tracks[ti];
        if (!nets.count(tr.net) || !(tr.width_nm > 0)) continue;
        // Points along the track and their tangents.
        std::function<void(double, double&, double&, double&, double&)> at;
        const double len = track_length(tr);
        if (tr.centre) {
            const double cx = double(tr.centre->x), cy = double(tr.centre->y);
            const double r = std::hypot(tr.start.x - cx, tr.start.y - cy);
            const double a0 = std::atan2(tr.start.y - cy, tr.start.x - cx);
            const double sweep = track_sweep(tr);
            at = [=](double s, double& x, double& y, double& ux, double& uy) {
                const double a = a0 + sweep * (len > 0 ? s / len : 0);
                x = cx + r * std::cos(a);
                y = cy + r * std::sin(a);
                const double d = sweep > 0 ? 1.0 : -1.0;
                ux = -std::sin(a) * d;
                uy = std::cos(a) * d;
            };
        } else {
            const double dx = double(tr.end.x - tr.start.x), dy = double(tr.end.y - tr.start.y);
            at = [=](double s, double& x, double& y, double& ux, double& uy) {
                x = tr.start.x + (len > 0 ? dx * s / len : 0);
                y = tr.start.y + (len > 0 ? dy * s / len : 0);
                ux = len > 0 ? dx / len : 1;
                uy = len > 0 ? dy / len : 0;
            };
        }
        if (!(len > 0)) continue;
        // Cuts clear of the ends by a track width, where bends and pads are.
        const double clear = std::min(tr.width_nm, 0.5 * len);
        const double span = len - 2 * clear;
        const int k = std::max(1, static_cast<int>(std::ceil(span / s_.step_nm)));
        for (int i = 0; i < k; ++i) {
            Task t;
            t.net = tr.net;
            t.layer = tr.layer;
            t.track = ti;
            t.width = tr.width_nm;
            t.along = clear + span * (i + 0.5) / k;
            t.weight = len / k;
            at(t.along, t.px, t.py, t.ux, t.uy);
            tasks.push_back(t);
        }
    }
    // Cuts, and the distinct cross-sections among them.
    std::vector<std::optional<Cut>> cuts(tasks.size());
    parallel_for(tasks.size(), [&](std::size_t i) { cuts[i] = make_cut(tasks[i]); });
    std::unordered_map<std::string, std::size_t> by_key;
    std::vector<std::size_t> cs_of(tasks.size(), SIZE_MAX);
    std::vector<const Cut*> distinct;
    for (std::size_t i = 0; i < tasks.size(); ++i) {
        if (!cuts[i] || cuts[i]->no_reference) continue;
        auto [it, fresh] = by_key.emplace(cuts[i]->key, distinct.size());
        if (fresh) distinct.push_back(&*cuts[i]);
        cs_of[i] = it->second;
    }
    res.cuts = tasks.size();
    std::vector<CsResult> solved(distinct.size());
    CsSettings cs;
    cs.target_error = s_.target_error;
    cs.margin = 0.5;
    cs.frequencies_hz = {s_.reference_hz};
    std::vector<std::string> errors(distinct.size());
    parallel_for(distinct.size(), [&](std::size_t k) {
        try {
            CsSettings ck = cs;
            for (const auto& r : distinct[k]->roles)
                if (r != "reference") ck.weights.push_back(r == "self" || r == "partner" ? 1.0 : 0.1);
            solved[k] = solve_cross_section(distinct[k]->conductors, distinct[k]->dielectrics, ck);
            if (s_.dispersion && distinct[k]->microstrip) {
                solved[k].model.microstrip = *distinct[k]->microstrip;
                for (auto& fq : solved[k].frequencies) fq = solved[k].model.at(fq.frequency_hz);
            }
        } catch (const std::exception& e) {
            errors[k] = e.what();
        }
    });
    for (std::size_t k = 0; k < distinct.size(); ++k) {
        if (!errors[k].empty()) warnings_.push_back("cross-section " + distinct[k]->key + ": " + errors[k]);
        res.cross_sections.push_back({distinct[k]->roles, std::move(solved[k])});
    }
    // Samples.
    for (std::size_t i = 0; i < tasks.size(); ++i) {
        const auto& t = tasks[i];
        SiSample s;
        s.net = t.net;
        s.layer = t.layer;
        s.x_nm = t.px;
        s.y_nm = t.py;
        if (!cuts[i]) continue;
        s.track = t.track;
        s.along_nm = t.along;
        for (const auto& o : cuts[i]->others)
            s.others.push_back({o.net, o.index, t.px - t.uy * o.s, t.py + t.ux * o.s});
        s.reference_layers = cuts[i]->reference_layers;
        s.no_reference = cuts[i]->no_reference;
        if (cs_of[i] != SIZE_MAX && errors[cs_of[i]].empty()) {
            s.cross_section = cs_of[i];
            const auto& r = res.cross_sections[cs_of[i]].result;
            const auto& f = r.frequencies.front();
            const std::size_t nc = r.conductors.size();
            const auto q = line_quantities(f, nc, 0);
            s.z0 = q.z0.real();
            s.epsilon_eff = q.epsilon_eff;
            s.delay_ps_per_mm = q.delay_s_per_m * 1e9;
            s.loss_db_per_m = q.loss_db_per_m;
            if (cuts[i]->partner && nc >= 2) {
                const auto p = pair_quantities(f, nc, 0, 1);
                s.z_differential = p.z_differential;
                s.z_common = p.z_common;
                s.coupled = true;
                s.partner_x_nm = t.px - t.uy * cuts[i]->partner_s;
                s.partner_y_nm = t.py + t.ux * cuts[i]->partner_s;
            }
        }
        res.samples.push_back(s);
    }
    std::vector<std::size_t> task_of;  // sample -> task
    for (std::size_t i = 0; i < tasks.size(); ++i)
        if (cuts[i]) task_of.push_back(i);

    // Distances along each net from its first pad: a graph of track ends, joined where
    // they meet, in a pad of the net or in a plated hole of the net.
    std::map<int, std::vector<std::size_t>> tracks_of;
    for (std::size_t ti = 0; ti < b_.tracks.size(); ++ti)
        if (nets.count(b_.tracks[ti].net)) tracks_of[b_.tracks[ti].net].push_back(ti);
    std::map<std::size_t, std::array<double, 2>> end_distance;  // track -> distance of start, end
    for (const auto& [net, list] : tracks_of) {
        std::map<std::tuple<int, std::int64_t, std::int64_t>, std::size_t> node_of;
        std::vector<std::vector<std::pair<std::size_t, double>>> adj;
        auto node = [&](int layer, Point64 p) {
            auto key = std::make_tuple(layer, p.x / 1000, p.y / 1000);
            auto [it, fresh] = node_of.emplace(key, adj.size());
            if (fresh) adj.emplace_back();
            return it->second;
        };
        auto link = [&](std::size_t a, std::size_t c, double w) {
            adj[a].push_back({c, w});
            adj[c].push_back({a, w});
        };
        std::vector<std::array<std::size_t, 2>> ends;
        for (auto ti : list) {
            const auto& tr = b_.tracks[ti];
            const double len = track_length(tr);
            const std::size_t a = node(tr.layer, tr.start), e = node(tr.layer, tr.end);
            link(a, e, len);
            ends.push_back({a, e});
        }
        // Ends inside a pad or a plated hole of the net are joined to it.
        std::vector<std::size_t> pad_nodes;
        std::size_t start = SIZE_MAX;
        std::vector<const Pad*> pads;
        for (const auto& p : b_.pads)
            if (p.net == net) pads.push_back(&p);
        std::sort(pads.begin(), pads.end(), [](const Pad* x, const Pad* y) {
            return std::tie(x->component, x->pin, x->layer) < std::tie(y->component, y->pin, y->layer);
        });
        std::map<std::pair<std::string, std::string>, std::size_t> pin_node;
        for (const auto* p : pads) {
            auto [it, fresh] = pin_node.emplace(std::pair{p->component, p->pin}, adj.size());
            if (fresh) adj.emplace_back();
            const std::size_t pn = it->second;
            if (start == SIZE_MAX) start = pn;
            const Rect64 pb = bounds(p->shape);
            for (std::size_t k = 0; k < list.size(); ++k) {
                const auto& tr = b_.tracks[list[k]];
                if (tr.layer != p->layer) continue;
                for (int e = 0; e < 2; ++e) {
                    const Point64 q = e ? tr.end : tr.start;
                    if (q.x >= pb.left && q.x <= pb.right && q.y >= pb.top && q.y <= pb.bottom && contains(p->shape, q))
                        link(ends[k][static_cast<std::size_t>(e)], pn, 0.0);
                }
            }
        }
        for (const auto& d : b_.drills) {
            if (!d.plated || d.net != net) continue;
            const std::size_t dn = adj.size();
            adj.emplace_back();
            const double r = 0.5 * d.diameter_nm + 1000.0;
            for (std::size_t k = 0; k < list.size(); ++k) {
                const auto& tr = b_.tracks[list[k]];
                for (int e = 0; e < 2; ++e) {
                    const Point64 q = e ? tr.end : tr.start;
                    if (std::hypot(double(q.x - d.at.x), double(q.y - d.at.y)) <= std::max(r, 0.5 * tr.width_nm))
                        link(ends[k][static_cast<std::size_t>(e)], dn, 0.0);
                }
            }
        }
        if (start == SIZE_MAX && !ends.empty()) start = ends.front()[0];
        std::vector<double> dist(adj.size(), INFINITY);
        if (start != SIZE_MAX) {
            using Item = std::pair<double, std::size_t>;
            std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
            dist[start] = 0;
            pq.push({0, start});
            while (!pq.empty()) {
                const auto [d, v] = pq.top();
                pq.pop();
                if (d > dist[v]) continue;
                for (const auto& [w, l] : adj[v])
                    if (d + l < dist[w]) {
                        dist[w] = d + l;
                        pq.push({dist[w], w});
                    }
            }
        }
        for (std::size_t k = 0; k < list.size(); ++k) end_distance[list[k]] = {dist[ends[k][0]], dist[ends[k][1]]};
    }
    for (std::size_t k = 0; k < res.samples.size(); ++k) {
        const auto& t = tasks[task_of[k]];
        const auto it = end_distance.find(t.track);
        double d = -1;
        if (it != end_distance.end()) {
            const double total = track_length(b_.tracks[t.track]);
            const double a = it->second[0] + t.along, c = it->second[1] + (total - t.along);
            d = std::min(a, c);
            if (!std::isfinite(d)) d = -1;
        }
        res.samples[k].distance_nm = d;
    }

    // Nets: totals, ranges and targets.
    auto target_of = [&](const std::string& name, bool differential) -> std::optional<SiTarget> {
        for (const auto& t : s_.targets)
            if (t.differential == differential && (t.nets.empty() || matches(t.nets, name))) return t;
        return std::nullopt;
    };
    std::map<int, std::vector<std::size_t>> samples_of;
    for (std::size_t k = 0; k < res.samples.size(); ++k) samples_of[res.samples[k].net].push_back(k);
    for (int id : nets) {
        SiNet sn;
        sn.name = b_.nets[static_cast<std::size_t>(id)];
        if (auto it = partner_.find(id); it != partner_.end()) sn.partner = b_.nets[static_cast<std::size_t>(it->second)];
        std::set<int> layers;
        double wz = 0, wzd = 0, sz = 0, szd = 0;
        sn.z_min = sn.zd_min = INFINITY;
        sn.z_max = sn.zd_max = -INFINITY;
        for (auto k : samples_of[id]) {
            const auto& s = res.samples[k];
            const auto& t = tasks[task_of[k]];
            sn.length_nm += t.weight;
            layers.insert(s.layer);
            if (!(s.z0 > 0)) continue;
            sn.delay_ps += s.delay_ps_per_mm * t.weight * 1e-6;
            sn.loss_db += s.loss_db_per_m * t.weight * 1e-9;
            sn.z_min = std::min(sn.z_min, s.z0);
            sn.z_max = std::max(sn.z_max, s.z0);
            wz += t.weight;
            sz += s.z0 * t.weight;
            if (s.z_differential > 0) {
                sn.zd_min = std::min(sn.zd_min, s.z_differential);
                sn.zd_max = std::max(sn.zd_max, s.z_differential);
                wzd += t.weight;
                szd += s.z_differential * t.weight;
                sn.coupled_nm += t.weight;
            }
        }
        sn.z_mean = wz > 0 ? sz / wz : 0;
        sn.zd_mean = wzd > 0 ? szd / wzd : 0;
        if (!(wz > 0)) sn.z_min = sn.z_max = 0;
        if (!(wzd > 0)) sn.zd_min = sn.zd_max = 0;
        sn.layers.assign(layers.begin(), layers.end());
        sn.target = target_of(sn.name, !sn.partner.empty() && sn.coupled_nm > 0);
        if (!sn.target) sn.target = target_of(sn.name, false);
        res.nets.push_back(sn);
    }
    // Crosstalk: per victim and aggressor, the coupling coefficients of the victim's cuts
    // where the aggressor's copper lies in them (pair partners excluded).
    {
        struct Acc {
            double len = 0, delay = 0, kb = 0, kbmax = 0, kf = 0, x = 0, y = 0;
            std::map<int, double> layers;
        };
        std::map<std::pair<int, int>, Acc> acc;  // victim, aggressor
        for (std::size_t k = 0; k < res.samples.size(); ++k) {
            const auto& s = res.samples[k];
            if (!(s.z0 > 0) || s.others.empty()) continue;
            const auto& r = res.cross_sections[s.cross_section].result;
            const auto& f = r.frequencies.front();
            const double ls = f.l[0], cs = f.c[0];
            const double w = tasks[task_of[k]].weight;  // nm
            const auto pit = partner_.find(s.net);
            for (const auto& o : s.others) {
                if (pit != partner_.end() && pit->second == o.net) continue;
                if (o.net == s.net || o.net == 0) continue;
                const double lm = f.l[o.conductor], cm = -f.c[o.conductor];
                const double kb = 0.25 * (lm / ls + cm / cs);
                const double kf = 0.5 * (cm / cs - lm / ls) * std::sqrt(ls * cs);  // s/m
                auto& a = acc[{s.net, o.net}];
                a.len += w;
                a.delay += w * s.delay_ps_per_mm * 1e-6;
                a.kb += kb * w;
                if (kb > a.kbmax) {
                    a.kbmax = kb;
                    a.x = s.x_nm;
                    a.y = s.y_nm;
                }
                a.kf += kf * w * 1e-9;  // s
                a.layers[s.layer] += w;
            }
        }
        const double tr = s_.rise_time_s;
        std::map<std::pair<std::string, std::string>, std::pair<double, double>> where;
        for (const auto& [key, a] : acc) {
            SiCoupling c;
            c.victim = b_.nets[static_cast<std::size_t>(key.first)];
            c.aggressor = b_.nets[static_cast<std::size_t>(key.second)];
            c.length_nm = a.len;
            c.delay_ps = a.delay;
            c.kb_mean = a.kb / a.len;
            c.kb_max = a.kbmax;
            c.kf_ps_per_mm = a.kf / (a.len * 1e-9) * 1e9;  // s/m to ps/mm
            c.layer = std::max_element(a.layers.begin(), a.layers.end(), [](auto& x, auto& y) { return x.second < y.second; })->first;
            // Near end: saturates once the coupled run's round trip outlasts the rise.
            c.next = c.kb_mean * std::min(1.0, 2 * c.delay_ps * 1e-12 / tr);
            // Far end: proportional to the coupled length and the slope of the step.
            c.fext = std::min(1.0, std::abs(a.kf) / tr);
            res.coupling.push_back(c);
            where[{c.victim, c.aggressor}] = {a.x, a.y};
        }
        std::sort(res.coupling.begin(), res.coupling.end(),
                  [](const SiCoupling& x, const SiCoupling& y) { return std::max(x.next, x.fext) > std::max(y.next, y.fext); });
        for (const auto& c : res.coupling)
            if (std::max(c.next, c.fext) > s_.crosstalk_threshold) {
                char text[200];
                std::snprintf(text, sizeof text, "from %s: near end %.1f %%, far end %.1f %% over %.2f mm", c.aggressor.c_str(), 100 * c.next,
                              100 * c.fext, c.length_nm * 1e-6);
                const auto [x, y] = where[{c.victim, c.aggressor}];
                res.findings.push_back({c.victim, "crosstalk", c.layer, x, y, c.length_nm, std::max(c.next, c.fext), s_.crosstalk_threshold, text});
            }
    }
    // Findings, per net along its samples in order of distance.
    for (const auto& sn : res.nets) {
        const int id = b_.find_net(sn.name);
        auto list = samples_of[id];
        std::sort(list.begin(), list.end(), [&](std::size_t a, std::size_t c) {
            return std::tie(res.samples[a].distance_nm, tasks[task_of[a]].track, tasks[task_of[a]].along) <
                   std::tie(res.samples[c].distance_nm, tasks[task_of[c]].track, tasks[task_of[c]].along);
        });
        // Runs of consecutive samples with the same finding are one finding; runs shorter
        // than a fifth (finding_rise_fraction) of the rise in delay are electrically short.
        auto push = [&](SiFinding f) {
            if (f.length_nm < s_.min_finding_nm) return;
            f.electrically_short = f.delay_ps < s_.finding_rise_fraction * s_.finding_rise_s * 1e12;
            res.findings.push_back(std::move(f));
        };
        auto run_findings = [&](const std::string& kind, auto test, auto value, double target) {
            SiFinding cur;
            bool open = false;
            for (auto k : list) {
                const auto& s = res.samples[k];
                const bool hit = test(s);
                if (hit) {
                    const double v = value(s);
                    if (!open) {
                        cur = SiFinding{sn.name, kind, s.layer, s.x_nm, s.y_nm, 0, v, target, ""};
                        open = true;
                    } else if (std::abs(v - target) > std::abs(cur.value - target)) {
                        cur.value = v;
                        cur.x_nm = s.x_nm;
                        cur.y_nm = s.y_nm;
                        cur.layer = s.layer;
                    }
                    cur.length_nm += tasks[task_of[k]].weight;
                    cur.delay_ps += tasks[task_of[k]].weight * 1e-6 * s.delay_ps_per_mm;
                } else if (open) {
                    push(cur);
                    open = false;
                }
            }
            if (open) push(cur);
        };
        if (sn.target && !sn.target->differential && sn.target->ohm > 0) {
            const double t = sn.target->ohm, tol = sn.target->tolerance;
            run_findings("impedance", [&](const SiSample& s) { return s.z0 > 0 && std::abs(s.z0 - t) > tol * t; },
                         [](const SiSample& s) { return s.z0; }, t);
        }
        if (sn.target && sn.target->differential && sn.target->ohm > 0) {
            const double t = sn.target->ohm, tol = sn.target->tolerance;
            run_findings("differential impedance",
                         [&](const SiSample& s) { return s.z_differential > 0 && std::abs(s.z_differential - t) > tol * t; },
                         [](const SiSample& s) { return s.z_differential; }, t);
        }
        run_findings("no reference", [](const SiSample& s) { return s.no_reference; }, [](const SiSample&) { return 0.0; }, 0.0);
        if (!sn.partner.empty())
            run_findings("uncoupled", [](const SiSample& s) { return s.z0 > 0 && !(s.z_differential > 0); },
                         [](const SiSample&) { return 0.0; }, 0.0);
        // Changes of the reference layers along the net, between runs of samples with the
        // same reference; a run without reference layers that is electrically short (a via's
        // or pad's surroundings) joins the runs on either side.
        struct RefRun {
            const std::vector<int>* refs;
            double delay_ps;
            const SiSample* first;
        };
        std::vector<RefRun> runs;
        for (auto k : list) {
            const auto& s = res.samples[k];
            if (s.no_reference || s.distance_nm < 0) continue;
            if (runs.empty() || *runs.back().refs != s.reference_layers) runs.push_back({&s.reference_layers, 0.0, &s});
            runs.back().delay_ps += tasks[task_of[k]].weight * 1e-6 * s.delay_ps_per_mm;
        }
        const double short_ps = s_.finding_rise_fraction * s_.finding_rise_s * 1e12;
        std::vector<RefRun> kept;
        for (const auto& r : runs) {
            if (r.refs->empty() && r.delay_ps < short_ps) continue;
            if (!kept.empty() && *kept.back().refs == *r.refs) continue;
            kept.push_back(r);
        }
        auto names = [](const std::vector<int>& v) {
            std::string t;
            for (int l : v) t += (t.empty() ? "L" : ",L") + std::to_string(l);
            return t.empty() ? std::string("none") : t;
        };
        for (std::size_t i = 1; i < kept.size(); ++i) {
            const SiSample& s = *kept[i].first;
            std::string text = names(*kept[i - 1].refs) + " to " + names(*kept[i].refs);
            if (kept[i - 1].first->layer != s.layer)
                text += " (signal L" + std::to_string(kept[i - 1].first->layer) + " to L" + std::to_string(s.layer) + ")";
            res.findings.push_back({sn.name, "reference change", s.layer, s.x_nm, s.y_nm, 0, 0, 0, text});
        }
    }
    res.warnings = warnings_;
    std::size_t skipped = 0;
    for (const auto& c : cuts) skipped += !c;
    if (skipped) res.warnings.push_back(std::to_string(skipped) + " cuts found no copper of their track and were skipped");
    res.finding_rise_ps = s_.finding_rise_s * 1e12;
    res.finding_rise_fraction = s_.finding_rise_fraction;
    res.seconds = since(t0);
    return res;
}

}  // namespace

SiResult analyse_si(const Board& board, const SiSettings& settings) { return Analysis(board, settings).run(); }

std::string si_json(const Board& b, const SiResult& r) {
    using nlohmann::json;
    const auto zs = stackup_z(b);
    std::map<int, double> ztop;
    for (std::size_t i = 0; i < b.stackup.size(); ++i)
        if (b.stackup[i].type == LayerType::Copper) ztop[b.stackup[i].copper_index] = zs[i].top_nm * 1e-6;
    json j;
    j["nets"] = json::array();
    for (const auto& n : r.nets) {
        json e{{"name", n.name}, {"partner", n.partner}, {"length_mm", n.length_nm * 1e-6}, {"delay_ps", n.delay_ps},
               {"z_min", n.z_min}, {"z_mean", n.z_mean}, {"z_max", n.z_max}, {"zd_min", n.zd_min}, {"zd_mean", n.zd_mean},
               {"zd_max", n.zd_max}, {"coupled_mm", n.coupled_nm * 1e-6}, {"loss_db", n.loss_db}, {"layers", n.layers}};
        if (n.target) e["target"] = {{"ohm", n.target->ohm}, {"differential", n.target->differential}, {"tolerance", n.target->tolerance}};
        j["nets"].push_back(e);
    }
    j["samples"] = json::array();
    std::map<std::size_t, std::vector<const SiSample*>> by_track;
    for (const auto& s : r.samples) {
        by_track[s.track].push_back(&s);
        j["samples"].push_back({{"net", b.nets[static_cast<std::size_t>(s.net)]}, {"layer", s.layer},
                                {"position_mm", {s.x_nm * 1e-6, s.y_nm * 1e-6}}, {"distance_mm", s.distance_nm < 0 ? -1.0 : s.distance_nm * 1e-6},
                                {"z0", s.z0}, {"z_differential", s.z_differential}, {"z_common", s.z_common},
                                {"epsilon_eff", s.epsilon_eff}, {"delay_ps_per_mm", s.delay_ps_per_mm}, {"loss_db_per_m", s.loss_db_per_m},
                                {"reference_layers", s.reference_layers}, {"cross_section", s.cross_section},
                                {"no_reference", s.no_reference}, {"track", s.track}, {"along_mm", s.along_nm * 1e-6}});
    }
    // Tracks with their samples, as polylines (arcs in steps of at most 0.2 mm).
    j["tracks"] = json::array();
    for (const auto& [ti, list] : by_track) {
        const auto& tr = b.tracks[ti];
        json pts = json::array();
        if (tr.centre) {
            const double cx = double(tr.centre->x), cy = double(tr.centre->y);
            const double rad = std::hypot(double(tr.start.x) - cx, double(tr.start.y) - cy);
            const double a0 = std::atan2(double(tr.start.y) - cy, double(tr.start.x) - cx), sweep = track_sweep(tr);
            const int n = std::max(2, static_cast<int>(std::ceil(rad * std::abs(sweep) / 2.0e5)) + 1);
            for (int k = 0; k < n; ++k) {
                const double a = a0 + sweep * k / (n - 1);
                pts.push_back({(cx + rad * std::cos(a)) * 1e-6, (cy + rad * std::sin(a)) * 1e-6});
            }
        } else {
            pts.push_back({tr.start.x * 1e-6, tr.start.y * 1e-6});
            pts.push_back({tr.end.x * 1e-6, tr.end.y * 1e-6});
        }
        json samples = json::array();
        for (const auto* s : list) samples.push_back({s->along_nm * 1e-6, s->z0, s->z_differential});
        j["tracks"].push_back({{"index", ti}, {"layer", tr.layer}, {"net", b.nets[static_cast<std::size_t>(tr.net)]},
                               {"width_mm", tr.width_nm * 1e-6}, {"length_mm", track_length(tr) * 1e-6}, {"z_mm", ztop[tr.layer]},
                               {"points_mm", pts}, {"samples", samples}});
    }
    j["coupling"] = json::array();
    for (const auto& c : r.coupling)
        j["coupling"].push_back({{"aggressor", c.aggressor}, {"victim", c.victim}, {"layer", c.layer}, {"length_mm", c.length_nm * 1e-6},
                                 {"delay_ps", c.delay_ps}, {"kb_max", c.kb_max}, {"kb_mean", c.kb_mean}, {"kf_ps_per_mm", c.kf_ps_per_mm},
                                 {"next", c.next}, {"fext", c.fext}});
    j["findings"] = json::array();
    for (const auto& f : r.findings)
        j["findings"].push_back({{"net", f.net}, {"kind", f.kind}, {"layer", f.layer}, {"position_mm", {f.x_nm * 1e-6, f.y_nm * 1e-6}},
                                 {"length_mm", f.length_nm * 1e-6}, {"value", f.value}, {"target", f.target}, {"text", f.text},
                                 {"delay_ps", f.delay_ps}, {"electrically_short", f.electrically_short}});
    j["cross_sections"] = json::array();
    for (const auto& c : r.cross_sections)
        j["cross_sections"].push_back({{"conductors", c.conductors}, {"triangles", c.result.triangles},
                                       {"estimated_error", c.result.estimated_error}, {"seconds", c.result.seconds},
                                       {"passes", c.result.passes}});
    j["cuts"] = r.cuts;
    j["seconds"] = r.seconds;
    j["warnings"] = r.warnings;
    return j.dump();
}

std::string si_summary(const Board& b, const SiResult& r) {
    std::ostringstream o;
    char line[512];
    std::snprintf(line, sizeof line, "%zu nets, %zu cuts, %zu cross-sections solved in %.1f s\n", r.nets.size(), r.cuts,
                  r.cross_sections.size(), r.seconds);
    o << line;
    o << "net                                layers   length mm  delay ps   Z ohm (min-mean-max)    Zdiff ohm (min-mean-max)   loss dB\n";
    for (const auto& n : r.nets) {
        std::string layers;
        for (int l : n.layers) layers += (layers.empty() ? "" : ",") + std::to_string(l);
        char zd[64] = "";
        if (n.zd_mean > 0) std::snprintf(zd, sizeof zd, "%6.1f %6.1f %6.1f", n.zd_min, n.zd_mean, n.zd_max);
        std::snprintf(line, sizeof line, "%-34s %-8s %9.2f %9.1f  %6.1f %6.1f %6.1f     %-25s %7.3f\n", n.name.substr(0, 34).c_str(),
                      layers.c_str(), n.length_nm * 1e-6, n.delay_ps, n.z_min, n.z_mean, n.z_max, zd, n.loss_db);
        o << line;
    }
    if (!r.coupling.empty()) {
        o << "crosstalk (largest; near and far end of the victim, for the step's rise time)\n";
        for (std::size_t k = 0; k < r.coupling.size() && k < 10; ++k) {
            const auto& c = r.coupling[k];
            std::snprintf(line, sizeof line, "  %-28s from %-28s %7.2f mm  Kb %.4f  near %5.2f %%  far %5.2f %%\n", c.victim.substr(0, 28).c_str(),
                          c.aggressor.substr(0, 28).c_str(), c.length_nm * 1e-6, c.kb_mean, 100 * c.next, 100 * c.fext);
            o << line;
        }
    }
    // Findings per net and kind: the electrically significant ones each, the short ones
    // counted with their worst value.
    std::size_t nshort = 0;
    for (const auto& f : r.findings) nshort += f.electrically_short;
    if (!r.findings.empty()) {
        std::snprintf(line, sizeof line, "findings: %zu, of which %zu electrically short (delay below %.2g of the %.0f ps rise), counted per net\n",
                      r.findings.size(), nshort, r.finding_rise_fraction, r.finding_rise_ps);
        o << line;
    }
    std::map<std::pair<std::string, std::string>, std::pair<std::size_t, const SiFinding*>> shorts;
    std::vector<std::pair<std::string, std::string>> order;
    for (const auto& f : r.findings) {
        if (!f.electrically_short) continue;
        auto& e = shorts[{f.net, f.kind}];
        if (!e.first) order.push_back({f.net, f.kind});
        ++e.first;
        if (!e.second || std::abs(f.value - f.target) > std::abs(e.second->value - e.second->target)) e.second = &f;
    }
    for (const auto& f : r.findings) {
        if (f.electrically_short) continue;
        if (f.kind == "reference change" || f.kind == "crosstalk")
            std::snprintf(line, sizeof line, "  %-22s %-30s L%d at (%.3f, %.3f) mm: %s\n", f.kind.c_str(), f.net.substr(0, 30).c_str(), f.layer,
                          f.x_nm * 1e-6, f.y_nm * 1e-6, f.text.c_str());
        else if (f.target > 0)
            std::snprintf(line, sizeof line, "  %-22s %-30s L%d at (%.3f, %.3f) mm: %.1f ohm against %.1f, over %.2f mm\n", f.kind.c_str(),
                          f.net.substr(0, 30).c_str(), f.layer, f.x_nm * 1e-6, f.y_nm * 1e-6, f.value, f.target, f.length_nm * 1e-6);
        else
            std::snprintf(line, sizeof line, "  %-22s %-30s L%d at (%.3f, %.3f) mm, over %.2f mm\n", f.kind.c_str(), f.net.substr(0, 30).c_str(),
                          f.layer, f.x_nm * 1e-6, f.y_nm * 1e-6, f.length_nm * 1e-6);
        o << line;
    }
    for (const auto& k : order) {
        const auto& [count, worst] = shorts[k];
        if (worst->target > 0)
            std::snprintf(line, sizeof line, "  %-22s %-30s %zu short, worst %.1f ohm against %.1f, L%d at (%.3f, %.3f) mm\n", k.second.c_str(),
                          k.first.substr(0, 30).c_str(), count, worst->value, worst->target, worst->layer, worst->x_nm * 1e-6, worst->y_nm * 1e-6);
        else
            std::snprintf(line, sizeof line, "  %-22s %-30s %zu short, L%d at (%.3f, %.3f) mm\n", k.second.c_str(), k.first.substr(0, 30).c_str(),
                          count, worst->layer, worst->x_nm * 1e-6, worst->y_nm * 1e-6);
        o << line;
    }
    for (const auto& w : r.warnings) o << "warning: " << w << "\n";
    (void)b;
    return o.str();
}

}  // namespace kairoq
