#include "kairoq/feature.hpp"

#include <cmath>
#include <unordered_map>

namespace kairoq {

using namespace Clipper2Lib;

std::map<std::string, PolygonSet> render_by_net(const std::vector<Feature>& features) {
    std::map<std::string, PolygonSet> done;
    std::map<std::string, Paths64> pending;
    Paths64 clear;

    auto flush_dark = [&]() {
        for (auto& [net, paths] : pending) {
            if (paths.empty()) continue;
            Paths64 all = to_paths(done[net]);
            all.insert(all.end(), std::make_move_iterator(paths.begin()),
                       std::make_move_iterator(paths.end()));
            done[net] = union_polygons(all);
            paths.clear();
        }
    };
    auto flush_clear = [&]() {
        if (clear.empty()) return;
        const PolygonSet cut = union_polygons(clear);
        for (auto& [net, polys] : done) polys = difference(polys, cut);
        clear.clear();
    };

    for (const auto& f : features) {
        if (f.shape.empty()) continue;
        if (f.dark) {
            if (!clear.empty()) flush_clear();
            auto& p = pending[f.net];
            p.insert(p.end(), f.shape.begin(), f.shape.end());
        } else {
            flush_dark();
            clear.insert(clear.end(), f.shape.begin(), f.shape.end());
        }
    }
    flush_dark();
    flush_clear();
    for (auto it = done.begin(); it != done.end();) {
        if (it->second.empty()) it = done.erase(it);
        else ++it;
    }
    return done;
}

PolygonSet render_all(const std::vector<Feature>& features) {
    std::vector<Feature> stripped;
    stripped.reserve(features.size());
    for (const auto& f : features) {
        Feature g;
        g.dark = f.dark;
        g.shape = f.shape;
        stripped.push_back(std::move(g));
    }
    auto by_net = render_by_net(stripped);
    return by_net.empty() ? PolygonSet{} : by_net.begin()->second;
}

namespace {

struct Piece {
    Path64 points;
    bool used = false;
};

struct PointHash {
    std::size_t operator()(const std::pair<std::int64_t, std::int64_t>& k) const {
        return std::hash<std::int64_t>()(k.first * 1000003LL) ^ std::hash<std::int64_t>()(k.second);
    }
};

}  // namespace

OutlineResult outline_from_features(const std::vector<Feature>& features, double tolerance_nm,
                                    double join_distance_nm) {
    OutlineResult result;
    Paths64 loops;
    std::vector<Piece> pieces;

    for (const auto& f : features) {
        if (f.kind == FeatureKind::Region) {
            for (const auto& p : f.shape) loops.push_back(p);
            continue;
        }
        if (f.kind != FeatureKind::Stroke) continue;
        Path64 pts = f.arc_centre ? arc_points(f.at, f.end, *f.arc_centre, f.arc_ccw, tolerance_nm, true)
                                  : Path64{f.at, f.end};
        if (f.arc_centre && f.at == f.end) {
            pts.pop_back();
            loops.push_back(std::move(pts));
            continue;
        }
        if (f.at == f.end) continue;
        pieces.push_back({std::move(pts)});
    }

    // Endpoints are matched on a grid of join_distance; neighbouring cells are searched
    // so that points closer than join_distance always meet.
    const double cell = std::max(join_distance_nm, 1.0);
    using Key = std::pair<std::int64_t, std::int64_t>;
    std::unordered_multimap<Key, std::pair<std::size_t, bool>, PointHash> ends;
    auto key = [&](Point64 p) {
        return Key{static_cast<std::int64_t>(std::floor(p.x / cell)),
                   static_cast<std::int64_t>(std::floor(p.y / cell))};
    };
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        ends.emplace(key(pieces[i].points.front()), std::make_pair(i, false));
        ends.emplace(key(pieces[i].points.back()), std::make_pair(i, true));
    }
    auto find_next = [&](Point64 p) -> std::optional<std::pair<std::size_t, bool>> {
        const Key k = key(p);
        double best = join_distance_nm;
        std::optional<std::pair<std::size_t, bool>> found;
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                auto range = ends.equal_range(Key{k.first + dx, k.second + dy});
                for (auto it = range.first; it != range.second; ++it) {
                    const auto [idx, at_end] = it->second;
                    if (pieces[idx].used) continue;
                    const Point64 q = at_end ? pieces[idx].points.back() : pieces[idx].points.front();
                    const double d = std::hypot(static_cast<double>(q.x - p.x),
                                                static_cast<double>(q.y - p.y));
                    if (d <= best) {
                        best = d;
                        found = it->second;
                    }
                }
            }
        }
        return found;
    };

    for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (pieces[i].used) continue;
        pieces[i].used = true;
        Path64 chain = pieces[i].points;
        const Point64 first = chain.front();
        bool closed = false;
        while (true) {
            const Point64 tail = chain.back();
            if (chain.size() > 2 && std::hypot(static_cast<double>(tail.x - first.x),
                                               static_cast<double>(tail.y - first.y)) <=
                                        join_distance_nm) {
                closed = true;
                break;
            }
            auto next = find_next(tail);
            if (!next) break;
            Piece& p = pieces[next->first];
            p.used = true;
            if (next->second) {
                chain.insert(chain.end(), p.points.rbegin() + 1, p.points.rend());
            } else {
                chain.insert(chain.end(), p.points.begin() + 1, p.points.end());
            }
        }
        if (closed) {
            chain.pop_back();
            loops.push_back(std::move(chain));
        } else {
            ++result.open_chains;
        }
    }
    result.polygons = even_odd_polygons(loops);
    return result;
}

}  // namespace kairoq
