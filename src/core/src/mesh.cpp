#include "kairoq/mesh.hpp"

#include <mapbox/earcut.hpp>

#include <array>
#include <bit>
#include <cstring>
#include <stdexcept>

namespace kairoq {

std::vector<std::uint32_t> triangulate(const Polygon& polygon) {
    using Pt = std::array<double, 2>;
    std::vector<std::vector<Pt>> rings;
    auto add = [&](const Path64& path) {
        std::vector<Pt> r;
        r.reserve(path.size());
        for (const auto& p : path) r.push_back({static_cast<double>(p.x), static_cast<double>(p.y)});
        rings.push_back(std::move(r));
    };
    add(polygon.outer);
    for (const auto& h : polygon.holes) add(h);
    return mapbox::earcut<std::uint32_t>(rings);
}

namespace {

void push_vertex(Mesh& m, double x, double y, double z, std::uint32_t id) {
    m.positions.push_back(static_cast<float>(x));
    m.positions.push_back(static_cast<float>(y));
    m.positions.push_back(static_cast<float>(z));
    m.ids.push_back(id);
}

void walls(Mesh& m, const Path64& ring, double z0, double z1, Point64 o, std::uint32_t id) {
    const std::size_t n = ring.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Point64 a = ring[i], b = ring[(i + 1) % n];
        const auto base = static_cast<std::uint32_t>(m.ids.size());
        push_vertex(m, to_mm(double(a.x - o.x)), to_mm(double(a.y - o.y)), z0, id);
        push_vertex(m, to_mm(double(b.x - o.x)), to_mm(double(b.y - o.y)), z0, id);
        push_vertex(m, to_mm(double(b.x - o.x)), to_mm(double(b.y - o.y)), z1, id);
        push_vertex(m, to_mm(double(a.x - o.x)), to_mm(double(a.y - o.y)), z1, id);
        m.indices.insert(m.indices.end(), {base, base + 1, base + 2, base, base + 2, base + 3});
    }
}

}  // namespace

void extrude(Mesh& m, const PolygonSet& polygons, double z0_nm, double z1_nm, Point64 o, std::uint32_t id) {
    const double z0 = to_mm(z0_nm), z1 = to_mm(z1_nm);
    for (const auto& poly : polygons) {
        const auto tri = triangulate(poly);
        for (double z : {z1, z0}) {
            const auto base = static_cast<std::uint32_t>(m.ids.size());
            auto add_ring = [&](const Path64& r) {
                for (const auto& p : r) push_vertex(m, to_mm(double(p.x - o.x)), to_mm(double(p.y - o.y)), z, id);
            };
            add_ring(poly.outer);
            for (const auto& h : poly.holes) add_ring(h);
            if (z == z1) {
                for (auto t : tri) m.indices.push_back(base + t);
            } else {
                // Bottom face with reversed winding.
                for (std::size_t k = 0; k + 2 < tri.size(); k += 3)
                    m.indices.insert(m.indices.end(), {base + tri[k], base + tri[k + 2], base + tri[k + 1]});
            }
        }
        walls(m, poly.outer, z0, z1, o, id);
        for (const auto& h : poly.holes) walls(m, h, z0, z1, o, id);
    }
}

Point64 display_origin(const Board& board) {
    if (board.outline.empty()) return Point64(0, 0);
    const Rect64 b = bounds(board.outline);
    return Point64((b.left + b.right) / 2, (b.top + b.bottom) / 2);
}

std::vector<std::string> mesh_keys(const Board& board) {
    std::vector<std::string> keys;
    for (const auto& l : board.copper) keys.push_back("copper/" + std::to_string(l.index));
    for (std::size_t i = 0; i < board.stackup.size(); ++i)
        if (board.stackup[i].type == LayerType::Dielectric) keys.push_back("dielectric/" + std::to_string(i));
    for (const auto& m : board.masks) keys.push_back(std::string("mask/") + to_string(m.side));
    if (!board.drills.empty()) keys.push_back("drills");
    return keys;
}

namespace {

// Copper layer index directly above (smaller index) a stackup position, or 0.
int copper_above(const Board& b, std::size_t pos) {
    for (std::size_t i = pos; i-- > 0;)
        if (b.stackup[i].type == LayerType::Copper) return b.stackup[i].copper_index;
    return 0;
}

int copper_below(const Board& b, std::size_t pos) {
    for (std::size_t i = pos + 1; i < b.stackup.size(); ++i)
        if (b.stackup[i].type == LayerType::Copper) return b.stackup[i].copper_index;
    return 0;
}

Path64 drill_outline(const Drill& d, double tol) {
    if (d.slot_end) return capsule(d.at, *d.slot_end, d.diameter_nm / 2, tol);
    return circle_path(d.at, d.diameter_nm / 2, tol);
}

std::size_t copper_stackup_position(const Board& b, int index) {
    for (std::size_t i = 0; i < b.stackup.size(); ++i)
        if (b.stackup[i].type == LayerType::Copper && b.stackup[i].copper_index == index) return i;
    throw std::runtime_error("copper layer " + std::to_string(index) + " is not in the stackup");
}

}  // namespace

Mesh board_mesh(const Board& b, const std::string& key) {
    Mesh m;
    const Point64 o = display_origin(b);
    const auto z = stackup_z(b);
    const double tol = b.source.arc_tolerance_nm > 0 ? b.source.arc_tolerance_nm * 5 : 5000.0;

    if (key.rfind("copper/", 0) == 0) {
        const int index = std::stoi(key.substr(7));
        const auto pos = copper_stackup_position(b, index);
        for (const auto& l : b.copper)
            if (l.index == index)
                for (const auto& nc : l.nets)
                    extrude(m, nc.polygons, z[pos].bottom_nm, z[pos].top_nm, o, static_cast<std::uint32_t>(nc.net));
        return m;
    }
    if (key.rfind("dielectric/", 0) == 0) {
        const auto pos = static_cast<std::size_t>(std::stoul(key.substr(11)));
        if (pos >= b.stackup.size()) throw std::runtime_error("no stackup layer " + key);
        const int above = copper_above(b, pos);
        const int below = copper_below(b, pos);
        Paths64 holes;
        for (const auto& d : b.drills) {
            // A hole passes through a dielectric between two copper layers it spans.
            const bool through = above == 0 || below == 0 ? (d.from_layer <= 1 && d.to_layer >= static_cast<int>(b.copper.size()))
                                                          : (d.from_layer <= above && d.to_layer >= below);
            if (through) holes.push_back(drill_outline(d, tol));
        }
        extrude(m, difference(b.outline, union_polygons(holes)), z[pos].bottom_nm, z[pos].top_nm, o, 0);
        return m;
    }
    if (key.rfind("mask/", 0) == 0) {
        const Side side = side_from_string(key.substr(5));
        double z0 = 0, z1 = 0;
        // The mask layer of that side in the stackup, else a 10 um film on the copper.
        std::optional<std::size_t> pos;
        for (std::size_t i = 0; i < b.stackup.size(); ++i) {
            if (b.stackup[i].type != LayerType::SolderMask) continue;
            const bool top = copper_below(b, i) != 0 && copper_above(b, i) == 0;
            if ((side == Side::Top) == top) pos = i;
        }
        if (pos && b.stackup[*pos].thickness_nm > 0) {
            z0 = z[*pos].bottom_nm;
            z1 = z[*pos].top_nm;
        } else if (side == Side::Top) {
            z0 = b.thickness_nm();
            z1 = z0 + 10000;
        } else {
            z1 = 0;
            z0 = -10000;
        }
        for (const auto& mk : b.masks)
            if (mk.side == side) extrude(m, difference(b.outline, mk.openings), z0, z1, o, 0);
        return m;
    }
    if (key == "drills") {
        for (const auto& d : b.drills) {
            if (!d.plated) continue;
            const double top = z[copper_stackup_position(b, std::max(d.from_layer, 1))].top_nm;
            const double bottom = z[copper_stackup_position(b, std::min<int>(std::max(d.to_layer, 1), b.copper.size()))].bottom_nm;
            walls(m, drill_outline(d, tol), to_mm(bottom), to_mm(top), o, static_cast<std::uint32_t>(d.net));
        }
        return m;
    }
    throw std::runtime_error("unknown mesh " + key);
}

std::string encode_mesh(const Mesh& m) {
    static_assert(std::endian::native == std::endian::little, "mesh encoding assumes a little-endian host");
    std::string out;
    const auto nv = static_cast<std::uint32_t>(m.ids.size());
    const auto ni = static_cast<std::uint32_t>(m.indices.size());
    out.append(reinterpret_cast<const char*>(&nv), 4);
    out.append(reinterpret_cast<const char*>(&ni), 4);
    out.append(reinterpret_cast<const char*>(m.positions.data()), m.positions.size() * sizeof(float));
    out.append(reinterpret_cast<const char*>(m.indices.data()), m.indices.size() * sizeof(std::uint32_t));
    out.append(reinterpret_cast<const char*>(m.ids.data()), m.ids.size() * sizeof(std::uint32_t));
    return out;
}

}  // namespace kairoq
