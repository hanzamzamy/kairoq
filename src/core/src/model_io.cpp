#include "kairoq/model_io.hpp"

#include "kairoq/archive.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <bit>
#include <cstring>
#include <filesystem>
#include <map>
#include <stdexcept>

namespace kairoq {

namespace fs = std::filesystem;
using nlohmann::json;

namespace {

constexpr const char* kFormat = "kairoq-board";

template <typename T>
void put(std::string& out, T v) {
    if constexpr (std::endian::native == std::endian::big) {
        unsigned char b[sizeof(T)];
        std::memcpy(b, &v, sizeof(T));
        std::reverse(b, b + sizeof(T));
        out.append(reinterpret_cast<const char*>(b), sizeof(T));
    } else {
        out.append(reinterpret_cast<const char*>(&v), sizeof(T));
    }
}

template <typename T>
T get(std::string_view data, std::size_t& pos, bool little) {
    if (pos + sizeof(T) > data.size()) throw std::runtime_error("truncated WKB data");
    unsigned char b[sizeof(T)];
    std::memcpy(b, data.data() + pos, sizeof(T));
    pos += sizeof(T);
    if (little != (std::endian::native == std::endian::little)) std::reverse(b, b + sizeof(T));
    T v;
    std::memcpy(&v, b, sizeof(T));
    return v;
}

void put_ring(std::string& out, const Path64& ring) {
    put<std::uint32_t>(out, static_cast<std::uint32_t>(ring.size() + 1));
    for (const auto& p : ring) {
        put<double>(out, to_mm(static_cast<double>(p.x)));
        put<double>(out, to_mm(static_cast<double>(p.y)));
    }
    put<double>(out, to_mm(static_cast<double>(ring.front().x)));
    put<double>(out, to_mm(static_cast<double>(ring.front().y)));
}

Path64 get_ring(std::string_view data, std::size_t& pos, bool little) {
    const auto n = get<std::uint32_t>(data, pos, little);
    Path64 ring;
    ring.reserve(n);
    for (std::uint32_t i = 0; i < n; ++i) {
        const double x = get<double>(data, pos, little);
        const double y = get<double>(data, pos, little);
        ring.emplace_back(nm_from_mm(x), nm_from_mm(y));
    }
    if (ring.size() > 1 && ring.front() == ring.back()) ring.pop_back();
    return ring;
}

Polygon get_polygon(std::string_view data, std::size_t& pos) {
    const bool little = data.at(pos++) == 1;
    const auto type = get<std::uint32_t>(data, pos, little);
    if (type != 3) throw std::runtime_error("WKB: expected a Polygon");
    const auto rings = get<std::uint32_t>(data, pos, little);
    Polygon poly;
    for (std::uint32_t r = 0; r < rings; ++r) {
        Path64 ring = get_ring(data, pos, little);
        if (r == 0) poly.outer = std::move(ring);
        else poly.holes.push_back(std::move(ring));
    }
    return poly;
}

// Byte size of encode_wkb(polygons), without encoding.
std::size_t wkb_size(const PolygonSet& polygons) {
    std::size_t n = 9;
    for (const auto& p : polygons) {
        n += 9 + 4 + 16 * (p.outer.size() + 1);
        for (const auto& h : p.holes) n += 4 + 16 * (h.size() + 1);
    }
    return n;
}

json point_json(Point64 p) { return json::array({to_mm(double(p.x)), to_mm(double(p.y))}); }

Point64 point_from(const json& j) { return Point64(nm_from_mm(j.at(0).get<double>()), nm_from_mm(j.at(1).get<double>())); }

template <typename T>
void put_optional(json& j, const char* key, const std::optional<T>& v) {
    if (v) j[key] = *v;
}

// Files of a model: name -> content.
using FileMap = std::map<std::string, std::string>;

FileMap read_model_files(const std::string& path) {
    FileMap files;
    if (fs::is_directory(path)) {
        for (const auto& e : fs::recursive_directory_iterator(path))
            if (e.is_regular_file())
                files[fs::relative(e.path(), path).generic_string()] = read_file(e.path().string());
    } else {
        for (auto& e : read_zip(path)) files[e.name] = std::move(e.data);
    }
    return files;
}

}  // namespace

std::string encode_wkb(const PolygonSet& polygons) {
    std::string out;
    out.push_back(1);
    put<std::uint32_t>(out, 6);
    put<std::uint32_t>(out, static_cast<std::uint32_t>(polygons.size()));
    for (const auto& p : polygons) {
        out.push_back(1);
        put<std::uint32_t>(out, 3);
        put<std::uint32_t>(out, static_cast<std::uint32_t>(1 + p.holes.size()));
        put_ring(out, p.outer);
        for (const auto& h : p.holes) put_ring(out, h);
    }
    return out;
}

PolygonSet decode_wkb(std::string_view data) {
    std::size_t pos = 0;
    if (data.empty()) return {};
    const bool little = data.at(pos++) == 1;
    const auto type = get<std::uint32_t>(data, pos, little);
    PolygonSet out;
    if (type == 3) {
        pos = 0;
        out.push_back(get_polygon(data, pos));
    } else if (type == 6) {
        const auto n = get<std::uint32_t>(data, pos, little);
        for (std::uint32_t i = 0; i < n; ++i) out.push_back(get_polygon(data, pos));
    } else {
        throw std::runtime_error("WKB: expected a Polygon or MultiPolygon");
    }
    return out;
}

json board_json(const Board& b, bool geometry) {
    json j;
    j["format"] = kFormat;
    j["version"] = kModelVersion;
    j["name"] = b.name;
    j["units"] = "mm";
    j["source"] = {{"importer", b.source.importer},
                   {"generator", b.source.generator},
                   {"files", b.source.files},
                   {"arc_tolerance_mm", to_mm(b.source.arc_tolerance_nm)}};

    json layers = json::array();
    for (const auto& l : b.stackup) {
        json e{{"name", l.name}, {"type", to_string(l.type)}, {"thickness_mm", to_mm(l.thickness_nm)}};
        if (!l.material.empty()) e["material"] = l.material;
        put_optional(e, "epsilon_r", l.epsilon_r);
        put_optional(e, "loss_tangent", l.loss_tangent);
        put_optional(e, "conductivity_s_per_m", l.conductivity_s_per_m);
        if (l.copper_index) e["copper_index"] = l.copper_index;
        layers.push_back(e);
    }
    j["stackup"] = {{"source", b.stackup_source}, {"layers", layers}};

    std::size_t offset = 0;
    json copper = json::array();
    for (const auto& l : b.copper) {
        json e{{"index", l.index}, {"name", l.name}, {"side", to_string(l.side)}};
        json nets = json::array();
        offset = 0;
        for (const auto& nc : l.nets) {
            json ne{{"net", nc.net}};
            if (geometry) {
                const std::size_t len = wkb_size(nc.polygons);
                ne["offset"] = offset;
                ne["length"] = len;
                offset += len;
            }
            nets.push_back(ne);
        }
        if (geometry) e["file"] = "copper/L" + std::to_string(l.index) + ".wkb";
        e["nets"] = nets;
        copper.push_back(e);
    }
    j["copper_layers"] = copper;

    json masks = json::array();
    for (const auto& m : b.masks) {
        json e{{"side", to_string(m.side)}};
        if (geometry) e["file"] = std::string("mask/") + to_string(m.side) + ".wkb";
        masks.push_back(e);
    }
    j["mask_layers"] = masks;
    if (geometry) j["outline"] = {{"file", "outline.wkb"}};

    j["nets"] = b.nets;
    j["net_source"] = b.net_source;

    json comps = json::array();
    for (const auto& c : b.components)
        comps.push_back({{"refdes", c.refdes},
                         {"value", c.value},
                         {"package", c.package},
                         {"position_mm", point_json(c.at)},
                         {"rotation_deg", c.rotation_deg},
                         {"side", to_string(c.side)},
                         {"position_source", c.position_source}});
    j["components"] = comps;

    json pads = json::array();
    std::map<int, std::size_t> pad_offset;
    for (const auto& p : b.pads) {
        json e{{"component", p.component},
               {"pin", p.pin},
               {"function", p.function},
               {"net", p.net},
               {"layer", p.layer},
               {"position_mm", point_json(p.at)}};
        if (geometry) {
            const std::size_t len = wkb_size(p.shape);
            e["file"] = "pads/L" + std::to_string(p.layer) + ".wkb";
            e["offset"] = pad_offset[p.layer];
            e["length"] = len;
            pad_offset[p.layer] += len;
        }
        pads.push_back(e);
    }
    j["pads"] = pads;

    json drills = json::array();
    for (const auto& d : b.drills) {
        json e{{"position_mm", point_json(d.at)},
               {"diameter_mm", to_mm(d.diameter_nm)},
               {"plated", d.plated},
               {"function", to_string(d.function)},
               {"from_layer", d.from_layer},
               {"to_layer", d.to_layer},
               {"net", d.net}};
        if (d.slot_end) e["slot_end_mm"] = point_json(*d.slot_end);
        drills.push_back(e);
    }
    j["drills"] = drills;
    json tracks = json::array();
    for (const auto& t : b.tracks) {
        json e{{"layer", t.layer}, {"net", t.net}, {"width_mm", to_mm(t.width_nm)}, {"start_mm", point_json(t.start)},
               {"end_mm", point_json(t.end)}};
        if (t.centre) {
            e["centre_mm"] = point_json(*t.centre);
            e["ccw"] = t.ccw;
        }
        tracks.push_back(e);
    }
    j["tracks"] = tracks;
    j["warnings"] = b.warnings;
    return j;
}

void save_board(const Board& b, const std::string& path) {
    std::vector<ArchiveEntry> entries;
    entries.push_back({"board.json", board_json(b, true).dump(1)});
    entries.push_back({"outline.wkb", encode_wkb(b.outline)});
    for (const auto& l : b.copper) {
        std::string data;
        for (const auto& nc : l.nets) data += encode_wkb(nc.polygons);
        entries.push_back({"copper/L" + std::to_string(l.index) + ".wkb", std::move(data)});
    }
    std::map<int, std::string> pad_files;
    for (const auto& p : b.pads) pad_files[p.layer] += encode_wkb(p.shape);
    for (auto& [layer, data] : pad_files) entries.push_back({"pads/L" + std::to_string(layer) + ".wkb", std::move(data)});
    for (const auto& m : b.masks) entries.push_back({std::string("mask/") + to_string(m.side) + ".wkb", encode_wkb(m.openings)});

    if (path.size() > 4 && path.substr(path.size() - 4) == ".zip") {
        write_zip(path, entries);
        return;
    }
    if (fs::exists(path) && !fs::is_directory(path)) throw std::runtime_error(path + " exists and is not a directory");
    if (fs::exists(fs::path(path) / "board.json")) {
        // Replace a previous model; leave anything else in the directory alone.
        for (const char* sub : {"copper", "pads", "mask"}) fs::remove_all(fs::path(path) / sub);
    }
    for (const auto& e : entries) {
        const fs::path p = fs::path(path) / e.name;
        fs::create_directories(p.parent_path());
        write_file(p.string(), e.data);
    }
}

bool is_board_model(const std::string& path) {
    if (fs::is_directory(path)) return fs::exists(fs::path(path) / "board.json");
    if (!is_zip_file(path)) return false;
    try {
        for (const auto& e : read_zip(path))
            if (e.name == "board.json") return true;
    } catch (const std::exception&) {
    }
    return false;
}

Board load_board(const std::string& path) {
    const FileMap files = read_model_files(path);
    auto file = [&](const std::string& name) -> const std::string& {
        auto it = files.find(name);
        if (it == files.end()) throw std::runtime_error("model file " + name + " is missing from " + path);
        return it->second;
    };
    const json j = json::parse(file("board.json"));
    if (j.value("format", "") != kFormat) throw std::runtime_error(path + " is not a " + std::string(kFormat) + " model");
    if (j.value("version", 0) > kModelVersion)
        throw std::runtime_error(path + " has model version " + std::to_string(j.value("version", 0)) +
                                 "; this build reads up to " + std::to_string(kModelVersion));
    auto slice = [&](const std::string& name, std::size_t offset, std::size_t length) {
        const std::string& data = file(name);
        if (offset + length > data.size()) throw std::runtime_error("byte range outside " + name);
        return decode_wkb(std::string_view(data).substr(offset, length));
    };

    Board b;
    b.name = j.value("name", "");
    const auto& src = j.at("source");
    b.source.importer = src.value("importer", "");
    b.source.generator = src.value("generator", "");
    b.source.files = src.value("files", std::vector<std::string>{});
    b.source.arc_tolerance_nm = src.value("arc_tolerance_mm", 0.0) * kNmPerMm;

    b.stackup_source = j.at("stackup").value("source", "");
    for (const auto& e : j.at("stackup").at("layers")) {
        StackupLayer l;
        l.name = e.value("name", "");
        l.type = layer_type_from_string(e.value("type", ""));
        l.thickness_nm = e.value("thickness_mm", 0.0) * kNmPerMm;
        l.material = e.value("material", "");
        if (e.contains("epsilon_r")) l.epsilon_r = e["epsilon_r"].get<double>();
        if (e.contains("loss_tangent")) l.loss_tangent = e["loss_tangent"].get<double>();
        if (e.contains("conductivity_s_per_m")) l.conductivity_s_per_m = e["conductivity_s_per_m"].get<double>();
        l.copper_index = e.value("copper_index", 0);
        b.stackup.push_back(std::move(l));
    }
    for (const auto& e : j.at("copper_layers")) {
        CopperLayer l;
        l.index = e.at("index").get<int>();
        l.name = e.value("name", "");
        l.side = side_from_string(e.value("side", "inner"));
        const std::string f = e.at("file").get<std::string>();
        for (const auto& ne : e.at("nets"))
            l.nets.push_back({ne.at("net").get<int>(), slice(f, ne.at("offset").get<std::size_t>(), ne.at("length").get<std::size_t>())});
        b.copper.push_back(std::move(l));
    }
    for (const auto& e : j.value("mask_layers", json::array())) {
        MaskLayer m;
        m.side = side_from_string(e.value("side", "top"));
        m.openings = decode_wkb(file(e.at("file").get<std::string>()));
        b.masks.push_back(std::move(m));
    }
    b.outline = decode_wkb(file(j.at("outline").at("file").get<std::string>()));
    b.nets = j.at("nets").get<std::vector<std::string>>();
    b.net_source = j.value("net_source", "");
    for (const auto& e : j.value("components", json::array())) {
        Component c;
        c.refdes = e.value("refdes", "");
        c.value = e.value("value", "");
        c.package = e.value("package", "");
        c.at = point_from(e.at("position_mm"));
        c.rotation_deg = e.value("rotation_deg", 0.0);
        c.side = side_from_string(e.value("side", "top"));
        c.position_source = e.value("position_source", "");
        b.components.push_back(std::move(c));
    }
    for (const auto& e : j.value("pads", json::array())) {
        Pad p;
        p.component = e.value("component", "");
        p.pin = e.value("pin", "");
        p.function = e.value("function", "");
        p.net = e.value("net", 0);
        p.layer = e.value("layer", 0);
        p.at = point_from(e.at("position_mm"));
        if (e.contains("file")) p.shape = slice(e["file"], e.at("offset").get<std::size_t>(), e.at("length").get<std::size_t>());
        b.pads.push_back(std::move(p));
    }
    for (const auto& e : j.value("drills", json::array())) {
        Drill d;
        d.at = point_from(e.at("position_mm"));
        if (e.contains("slot_end_mm")) d.slot_end = point_from(e["slot_end_mm"]);
        d.diameter_nm = e.value("diameter_mm", 0.0) * kNmPerMm;
        d.plated = e.value("plated", true);
        d.function = drill_function_from_string(e.value("function", "unknown"));
        d.from_layer = e.value("from_layer", 1);
        d.to_layer = e.value("to_layer", 1);
        d.net = e.value("net", 0);
        b.drills.push_back(d);
    }
    for (const auto& e : j.value("tracks", json::array())) {  // version 2
        Track t;
        t.layer = e.value("layer", 0);
        t.net = e.value("net", 0);
        t.width_nm = e.value("width_mm", 0.0) * kNmPerMm;
        t.start = point_from(e.at("start_mm"));
        t.end = point_from(e.at("end_mm"));
        if (e.contains("centre_mm")) {
            t.centre = point_from(e["centre_mm"]);
            t.ccw = e.value("ccw", false);
        }
        b.tracks.push_back(t);
    }
    b.warnings = j.value("warnings", std::vector<std::string>{});
    return b;
}

}  // namespace kairoq
