#include "build.hpp"

#include <algorithm>
#include <functional>
#include <cstdio>
#include <map>
#include <cmath>
#include <cctype>
#include <memory>
#include <numeric>
#include <set>

namespace kairoq::detail {

namespace {

struct UnionFind {
    std::vector<std::size_t> parent;
    explicit UnionFind(std::size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    std::size_t find(std::size_t a) {
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        return a;
    }
    void unite(std::size_t a, std::size_t b) { parent[find(a)] = find(b); }
};

// Polygons of one layer with a grid index, for point lookups.
struct LayerLookup {
    std::vector<const Polygon*> polygons;
    std::vector<std::size_t> ids;  // caller-defined id per polygon
    std::unique_ptr<GridIndex> index;

    void build() {
        std::vector<Rect64> boxes;
        boxes.reserve(polygons.size());
        for (const auto* p : polygons) boxes.push_back(bounds(*p));
        index = std::make_unique<GridIndex>(boxes, 1.0e6);
    }
    // Ids of the polygons containing the point.
    std::vector<std::size_t> at(Point64 p) const {
        std::vector<std::size_t> out;
        if (!index) return out;
        for (std::size_t i : index->query(p))
            if (contains(*polygons[i], p)) out.push_back(ids[i]);
        return out;
    }
};

int clamp_layer(int l, int n) { return std::clamp(l, 1, std::max(n, 1)); }

}  // namespace

std::string normalise_net_name(const std::string& name) {
    std::string out = name;
    for (auto& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return out;
}

std::vector<StackupLayer> default_stackup(int n, double board_thickness_nm) {
    if (board_thickness_nm <= 0) board_thickness_nm = 1.6 * kNmPerMm;
    const double cu = 0.035 * kNmPerMm;
    const double diel = n > 1 ? (board_thickness_nm - n * cu) / (n - 1) : board_thickness_nm - cu;
    std::vector<StackupLayer> out;
    for (int i = 1; i <= n; ++i) {
        StackupLayer c;
        c.name = "L" + std::to_string(i);
        c.type = LayerType::Copper;
        c.thickness_nm = cu;
        c.material = "copper";
        c.copper_index = i;
        out.push_back(c);
        if (i < n || n == 1) {
            StackupLayer d;
            d.name = n == 1 ? "core" : "D" + std::to_string(i);
            d.type = LayerType::Dielectric;
            d.thickness_nm = diel;
            d.material = "FR-4";
            out.push_back(d);
        }
    }
    return out;
}

void set_layer_copper(Board& board, CopperLayer& layer, const std::map<std::string, PolygonSet>& by_net) {
    std::map<int, PolygonSet> merged;
    for (const auto& [name, polys] : by_net) {
        const std::string net = name == "N/C" ? std::string() : name;
        auto& dst = merged[board.net_id(net)];
        if (dst.empty()) dst = polys;
        else dst = unite(dst, polys);
    }
    layer.nets.clear();
    for (auto& [id, polys] : merged) layer.nets.push_back({id, std::move(polys)});
}

void assign_drill_nets(Board& board) {
    const int n = static_cast<int>(board.copper.size());
    std::vector<LayerLookup> lookup(n);
    for (int l = 0; l < n; ++l) {
        for (const auto& nc : board.copper[l].nets) {
            for (const auto& p : nc.polygons) {
                lookup[l].polygons.push_back(&p);
                lookup[l].ids.push_back(static_cast<std::size_t>(nc.net));
            }
        }
        lookup[l].build();
    }
    int conflicts = 0;
    int unconnected = 0;
    for (auto& d : board.drills) {
        if (!d.plated || d.net != 0) continue;
        std::set<std::size_t> found;
        for (int l = clamp_layer(d.from_layer, n); l <= clamp_layer(d.to_layer, n); ++l)
            for (std::size_t net : lookup[l - 1].at(d.at))
                if (net != 0) found.insert(net);
        if (found.empty()) {
            ++unconnected;
            continue;
        }
        d.net = static_cast<int>(*found.begin());
        if (found.size() > 1) ++conflicts;
    }
    if (conflicts)
        board.warnings.push_back(std::to_string(conflicts) + " plated holes touch copper of more than one net");
    if (unconnected)
        board.warnings.push_back(std::to_string(unconnected) + " plated holes touch no copper with a net");
}

bool register_netlist(Netlist& nl, const Board& board, std::string& note) {
    std::vector<Point64> drills;
    for (const auto& d : board.drills)
        if (d.plated) drills.push_back(d.at);
    std::vector<Point64> holes;
    for (const auto& p : nl.points)
        if (p.hole_nm > 0) holes.push_back(p.at);
    if (drills.empty() || holes.empty()) {
        note = "the netlist has no holes to place it by; taken as it stands";
        return false;
    }
    // Votes for offsets drill - point (10 um bins), for y as given and y mirrored.
    const std::int64_t bin = 10000;
    std::set<std::pair<std::int64_t, std::int64_t>> drill_set;
    for (const auto& d : drills) drill_set.insert({d.x / bin, d.y / bin});
    auto matches = [&](int sy, Point64 off) {
        std::size_t n = 0;
        for (const auto& h : holes) {
            const std::int64_t x = h.x + off.x, y = sy * h.y + off.y;
            bool hit = false;
            for (std::int64_t dx = -1; dx <= 1 && !hit; ++dx)
                for (std::int64_t dy = -1; dy <= 1 && !hit; ++dy) hit = drill_set.count({x / bin + dx, y / bin + dy}) > 0;
            n += hit;
        }
        return n;
    };
    int best_sy = 1;
    Point64 best_off{0, 0};
    std::size_t best = matches(1, best_off);
    const std::size_t step = std::max<std::size_t>(1, holes.size() / 64);
    for (const int sy : {1, -1}) {
        std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> votes;
        for (std::size_t i = 0; i < holes.size(); i += step)
            for (const auto& d : drills) ++votes[{(d.x - holes[i].x) / bin, (d.y - sy * holes[i].y) / bin}];
        std::vector<std::pair<std::size_t, std::pair<std::int64_t, std::int64_t>>> top;
        for (const auto& [k, v] : votes) top.push_back({v, k});
        std::partial_sort(top.begin(), top.begin() + std::min<std::size_t>(5, top.size()), top.end(), std::greater<>());
        for (std::size_t t = 0; t < std::min<std::size_t>(5, top.size()); ++t) {
            // Refine the offset to the mean of the matched pairs near the voted bin.
            const Point64 guess{top[t].second.first * bin, top[t].second.second * bin};
            double sx = 0, syy = 0;
            std::size_t n = 0;
            for (std::size_t i = 0; i < holes.size(); i += step)
                for (const auto& d : drills) {
                    const double ex = double(d.x - holes[i].x - guess.x), ey = double(d.y - sy * holes[i].y - guess.y);
                    if (std::abs(ex) < 2.0 * bin && std::abs(ey) < 2.0 * bin) {
                        sx += ex;
                        syy += ey;
                        ++n;
                    }
                }
            const Point64 off = n ? Point64(guess.x + std::llround(sx / double(n)), guess.y + std::llround(syy / double(n))) : guess;
            const std::size_t m = matches(sy, off);
            if (m > best) {
                best = m;
                best_sy = sy;
                best_off = off;
            }
        }
    }
    if (2 * best < holes.size()) {
        note = "the netlist's holes do not land on the drills under any translation; taken as it stands";
        return false;
    }
    if (best_sy == 1 && best_off.x == 0 && best_off.y == 0) return true;
    for (auto& p : nl.points) p.at = Point64(p.at.x + best_off.x, best_sy * p.at.y + best_off.y);
    char text[200];
    std::snprintf(text, sizeof text, "the netlist is placed by an offset of (%.4f, %.4f) mm%s: %zu of its %zu holes land on drills",
                  best_off.x * 1e-6, best_off.y * 1e-6, best_sy < 0 ? " with y mirrored" : "", best, holes.size());
    note = text;
    return true;
}

void geometric_nets(Board& board, const std::vector<PolygonSet>& copper, const Netlist* netlist) {
    const int n = static_cast<int>(copper.size());
    std::vector<LayerLookup> lookup(n);
    std::vector<std::pair<int, const Polygon*>> nodes;
    for (int l = 0; l < n; ++l) {
        for (const auto& p : copper[l]) {
            lookup[l].polygons.push_back(&p);
            lookup[l].ids.push_back(nodes.size());
            nodes.emplace_back(l, &p);
        }
        lookup[l].build();
    }
    UnionFind uf(nodes.size());
    std::vector<long> drill_node(board.drills.size(), -1);
    for (std::size_t i = 0; i < board.drills.size(); ++i) {
        const auto& d = board.drills[i];
        if (!d.plated) continue;
        long first = -1;
        for (int l = clamp_layer(d.from_layer, n); l <= clamp_layer(d.to_layer, n); ++l) {
            for (std::size_t node : lookup[l - 1].at(d.at)) {
                if (first < 0) first = static_cast<long>(node);
                else uf.unite(node, static_cast<std::size_t>(first));
            }
        }
        drill_node[i] = first;
    }

    std::map<std::size_t, std::string> label;
    int conflicts = 0;
    if (netlist) {
        for (const auto& pt : netlist->points) {
            if (pt.net.empty() || pt.net == "N/C") continue;
            std::vector<int> layers;
            if (pt.access == 0) layers = {1, n};
            else layers = {clamp_layer(pt.access, n)};
            for (int l : layers) {
                for (std::size_t node : lookup[l - 1].at(pt.at)) {
                    const std::size_t root = uf.find(node);
                    auto [it, inserted] = label.emplace(root, pt.net);
                    if (!inserted && it->second != pt.net) ++conflicts;
                }
            }
        }
    }
    if (conflicts)
        board.warnings.push_back(std::to_string(conflicts) +
                                 " netlist points disagree with the geometric connectivity (shorts or missing copper)");

    std::map<std::size_t, int> root_net;
    int unnamed = 0;
    auto net_of = [&](std::size_t node) {
        const std::size_t root = uf.find(node);
        auto it = root_net.find(root);
        if (it != root_net.end()) return it->second;
        const auto l = label.find(root);
        const int id = board.net_id(l != label.end() ? l->second : "N$" + std::to_string(++unnamed));
        root_net[root] = id;
        return id;
    };
    for (int l = 0; l < n; ++l) {
        std::map<int, PolygonSet> by_net;
        for (std::size_t k = 0; k < lookup[l].polygons.size(); ++k)
            by_net[net_of(lookup[l].ids[k])].push_back(*lookup[l].polygons[k]);
        board.copper[l].nets.clear();
        for (auto& [id, polys] : by_net) board.copper[l].nets.push_back({id, std::move(polys)});
    }
    for (std::size_t i = 0; i < board.drills.size(); ++i)
        if (drill_node[i] >= 0) board.drills[i].net = net_of(static_cast<std::size_t>(drill_node[i]));
}

void cross_check_netlist(Board& board, const Netlist& netlist) {
    std::map<std::pair<std::string, std::string>, std::set<std::string>> pad_nets;
    for (const auto& p : board.pads) pad_nets[{p.component, p.pin}].insert(board.nets[p.net]);
    // IPC-D-356 fields are fixed width: designators keep 6 characters, net names 14.
    // KiCad marks a shortened net name with a "#<n>" suffix.
    auto candidates = [&](const NetlistPoint& pt) {
        std::set<std::string> nets;
        auto it = pad_nets.find({pt.refdes, pt.pin});
        if (it != pad_nets.end()) return it->second;
        if (pt.refdes.size() == 6)
            for (const auto& [key, names] : pad_nets)
                if (key.second == pt.pin && key.first.rfind(pt.refdes, 0) == 0) nets.insert(names.begin(), names.end());
        return nets;
    };
    auto same = [](const std::string& have_raw, const std::string& want_raw) {
        const std::string have = normalise_net_name(have_raw);
        std::string want = normalise_net_name(want_raw);
        if (have == want) return true;
        if (have.empty() && (want.empty() || want == "N/C")) return true;
        bool shortened = want.size() == 14;
        const std::size_t hash = want.rfind('#');
        if (hash != std::string::npos && hash > 0 && hash + 1 < want.size() &&
            want.find_first_not_of("0123456789", hash + 1) == std::string::npos) {
            want.resize(hash);
            shortened = true;
        }
        return shortened && have.size() > want.size() && have.ends_with(want);
    };
    int compared = 0, mismatched = 0, missing = 0;
    std::string example;
    for (const auto& pt : netlist.points) {
        if (pt.refdes.empty() || pt.pin.empty()) continue;  // vias and tooling holes
        const auto nets = candidates(pt);
        if (nets.empty()) {
            ++missing;
            continue;
        }
        ++compared;
        bool ok = false;
        for (const auto& name : nets) ok = ok || same(name, pt.net);
        if (!ok) {
            ++mismatched;
            if (example.empty()) example = pt.refdes + "-" + pt.pin + " (" + pt.net + ")";
        }
    }
    if (mismatched)
        board.warnings.push_back(std::to_string(mismatched) + " of " + std::to_string(compared) +
                                 " netlist points disagree with the X2 pad nets, first: " + example);
    if (missing)
        board.warnings.push_back(std::to_string(missing) + " netlist points have no matching pad");
}

void components_from_pads(Board& board) {
    struct Acc {
        double x = 0, y = 0;
        int count = 0;
        bool bottom_only = true;
    };
    std::map<std::string, Acc> acc;
    const int n = static_cast<int>(board.copper.size());
    for (const auto& p : board.pads) {
        if (p.component.empty()) continue;
        auto& a = acc[p.component];
        a.x += static_cast<double>(p.at.x);
        a.y += static_cast<double>(p.at.y);
        ++a.count;
        if (p.layer != n) a.bottom_only = false;
    }
    for (const auto& [ref, a] : acc) {
        Component c;
        c.refdes = ref;
        c.at = Point64(std::llround(a.x / a.count), std::llround(a.y / a.count));
        c.side = a.bottom_only && n > 1 ? Side::Bottom : Side::Top;
        c.position_source = "pads";
        board.components.push_back(std::move(c));
    }
}

void collect_tracks(Board& board, const std::vector<const std::vector<Feature>*>& layers) {
    board.tracks.clear();
    for (std::size_t l = 0; l < layers.size() && l < board.copper.size(); ++l) {
        if (!layers[l]) continue;
        LayerLookup lookup;
        for (const auto& nc : board.copper[l].nets)
            for (const auto& p : nc.polygons) {
                lookup.polygons.push_back(&p);
                lookup.ids.push_back(static_cast<std::size_t>(nc.net));
            }
        lookup.build();
        for (const auto& f : *layers[l]) {
            if (f.kind != FeatureKind::Stroke || !f.dark || f.is_pad() || !(f.stroke_width_nm > 0)) continue;
            Track t;
            t.layer = static_cast<int>(l) + 1;
            t.width_nm = f.stroke_width_nm;
            t.start = f.at;
            t.end = f.end;
            t.centre = f.arc_centre;
            t.ccw = f.arc_ccw;
            // The middle of the stroke lies in its own copper.
            Point64 mid((f.at.x + f.end.x) / 2, (f.at.y + f.end.y) / 2);
            if (f.arc_centre) {
                const double cx = double(f.arc_centre->x), cy = double(f.arc_centre->y);
                const double a0 = std::atan2(f.at.y - cy, f.at.x - cx), r = std::hypot(f.at.x - cx, f.at.y - cy);
                double sweep = std::atan2(f.end.y - cy, f.end.x - cx) - a0;
                if (f.arc_ccw && sweep <= 0) sweep += 2 * kPi;
                if (!f.arc_ccw && sweep >= 0) sweep -= 2 * kPi;
                mid = Point64(std::llround(cx + r * std::cos(a0 + sweep / 2)), std::llround(cy + r * std::sin(a0 + sweep / 2)));
            }
            const auto ids = lookup.at(mid);
            if (ids.size() == 1) t.net = static_cast<int>(ids[0]);
            else if (ids.size() > 1 && f.has_net) t.net = board.find_net(f.net == "N/C" ? std::string() : f.net);
            if (t.net < 0) t.net = 0;
            board.tracks.push_back(t);
        }
    }
}

void check_net_overlaps(Board& board) {
    for (const auto& layer : board.copper) {
        if (layer.nets.size() < 2) continue;
        double sum = 0.0;
        Paths64 all;
        for (const auto& nc : layer.nets) {
            sum += area_nm2(nc.polygons);
            const auto p = to_paths(nc.polygons);
            all.insert(all.end(), p.begin(), p.end());
        }
        const double overlap = sum - area_nm2(union_polygons(all));
        // 0.001 mm^2 absorbs the rounding of separately merged nets.
        if (overlap > 1.0e9)
            board.warnings.push_back("copper of different nets overlaps on layer " + layer.name + ": " +
                                     std::to_string(overlap / 1.0e12) + " mm^2");
    }
}

}  // namespace kairoq::detail
