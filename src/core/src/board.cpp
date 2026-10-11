#include "kairoq/board.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace kairoq {

const char* to_string(LayerType t) {
    switch (t) {
    case LayerType::Copper: return "copper";
    case LayerType::Dielectric: return "dielectric";
    case LayerType::SolderMask: return "soldermask";
    case LayerType::Silkscreen: return "silkscreen";
    case LayerType::SolderPaste: return "solderpaste";
    case LayerType::Other: return "other";
    }
    return "other";
}

const char* to_string(Side s) {
    switch (s) {
    case Side::Top: return "top";
    case Side::Inner: return "inner";
    case Side::Bottom: return "bottom";
    }
    return "inner";
}

const char* to_string(DrillFunction f) {
    switch (f) {
    case DrillFunction::Via: return "via";
    case DrillFunction::Component: return "component";
    case DrillFunction::Mechanical: return "mechanical";
    case DrillFunction::Unknown: return "unknown";
    }
    return "unknown";
}

LayerType layer_type_from_string(const std::string& s) {
    if (s == "copper") return LayerType::Copper;
    if (s == "dielectric") return LayerType::Dielectric;
    if (s == "soldermask") return LayerType::SolderMask;
    if (s == "silkscreen") return LayerType::Silkscreen;
    if (s == "solderpaste") return LayerType::SolderPaste;
    return LayerType::Other;
}

Side side_from_string(const std::string& s) {
    if (s == "top") return Side::Top;
    if (s == "bottom") return Side::Bottom;
    return Side::Inner;
}

DrillFunction drill_function_from_string(const std::string& s) {
    if (s == "via") return DrillFunction::Via;
    if (s == "component") return DrillFunction::Component;
    if (s == "mechanical") return DrillFunction::Mechanical;
    return DrillFunction::Unknown;
}

int Board::net_id(const std::string& name) {
    if (nets.empty()) nets.push_back("");
    if (name.empty()) return 0;
    const int found = find_net(name);
    if (found >= 0) return found;
    nets.push_back(name);
    return static_cast<int>(nets.size()) - 1;
}

int Board::find_net(const std::string& name) const {
    const auto it = std::find(nets.begin(), nets.end(), name);
    return it == nets.end() ? -1 : static_cast<int>(it - nets.begin());
}

const StackupLayer* Board::stackup_for_copper(int index) const {
    for (const auto& l : stackup)
        if (l.type == LayerType::Copper && l.copper_index == index) return &l;
    return nullptr;
}

double Board::thickness_nm() const {
    double t = 0.0;
    for (const auto& l : stackup) t += l.thickness_nm;
    return t;
}

Board crop_board(const Board& b, const Rect64& r) {
    if (!(r.right > r.left && r.bottom > r.top)) throw std::runtime_error("crop: empty rectangle");
    Path64 ring{{r.left, r.top}, {r.right, r.top}, {r.right, r.bottom}, {r.left, r.bottom}};
    make_positive(ring);
    const PolygonSet box{Polygon{ring, {}}};
    auto within = [&](const Rect64& q) { return q.left >= r.left && q.right <= r.right && q.top >= r.top && q.bottom <= r.bottom; };
    auto inside = [&](Point64 p) { return p.x >= r.left && p.x <= r.right && p.y >= r.top && p.y <= r.bottom; };
    Board out = b;
    out.outline = intersection(b.outline, box);
    if (out.outline.empty()) throw std::runtime_error("crop: the rectangle does not overlap the board");
    for (auto& cl : out.copper) {
        std::vector<NetCopper> kept;
        for (auto& nc : cl.nets) {
            nc.polygons = intersection(nc.polygons, box);
            if (!nc.polygons.empty()) kept.push_back(std::move(nc));
        }
        cl.nets = std::move(kept);
    }
    for (auto& m : out.masks) m.openings = intersection(m.openings, box);
    out.pads.clear();
    for (const auto& p : b.pads) {
        Pad q = p;
        q.shape = intersection(p.shape, box);
        if (!q.shape.empty()) out.pads.push_back(std::move(q));
    }
    out.drills.clear();
    std::size_t cut = 0;
    for (const auto& d : b.drills) {
        const auto rad = static_cast<std::int64_t>(std::ceil(d.diameter_nm / 2));
        Rect64 q(d.at.x - rad, d.at.y - rad, d.at.x + rad, d.at.y + rad);
        if (d.slot_end) q = Rect64(std::min(d.at.x, d.slot_end->x) - rad, std::min(d.at.y, d.slot_end->y) - rad,
                                   std::max(d.at.x, d.slot_end->x) + rad, std::max(d.at.y, d.slot_end->y) + rad);
        if (within(q)) out.drills.push_back(d);
        else if (inside(d.at)) ++cut;
    }
    out.tracks.clear();
    for (const auto& t : b.tracks)
        if (inside(t.start) && inside(t.end)) out.tracks.push_back(t);
    out.components.clear();
    for (const auto& c : b.components)
        if (inside(c.at)) out.components.push_back(c);
    char line[200];
    std::snprintf(line, sizeof line, "cropped to %.3f, %.3f to %.3f, %.3f mm", to_mm(double(r.left)), to_mm(double(r.top)), to_mm(double(r.right)),
                  to_mm(double(r.bottom)));
    out.warnings.push_back(line);
    if (cut) out.warnings.push_back(std::to_string(cut) + " holes cut by the crop's edge are left out");
    return out;
}

std::vector<LayerZ> stackup_z(const Board& board) {
    std::vector<LayerZ> z(board.stackup.size());
    double top = board.thickness_nm();
    for (std::size_t i = 0; i < board.stackup.size(); ++i) {
        z[i].top_nm = top;
        top -= board.stackup[i].thickness_nm;
        z[i].bottom_nm = top;
    }
    return z;
}

}  // namespace kairoq
