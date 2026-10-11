// IPC-2581 importer. Element names are matched without namespace prefixes.
#include "kairoq/ipc2581.hpp"

#include "build.hpp"
#include "kairoq/archive.hpp"
#include "kairoq/feature.hpp"

#include <pugixml.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <map>
#include <set>
#include <stdexcept>

namespace kairoq {

namespace fs = std::filesystem;

namespace {

const char* local_name(const pugi::xml_node& n) {
    const char* name = n.name();
    const char* colon = std::strrchr(name, ':');
    return colon ? colon + 1 : name;
}

bool is(const pugi::xml_node& n, const char* name) { return std::strcmp(local_name(n), name) == 0; }

pugi::xml_node child(const pugi::xml_node& n, const char* name) {
    for (auto c : n.children())
        if (is(c, name)) return c;
    return {};
}

std::vector<pugi::xml_node> children(const pugi::xml_node& n, const char* name) {
    std::vector<pugi::xml_node> out;
    for (auto c : n.children())
        if (is(c, name)) out.push_back(c);
    return out;
}

std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::string load_xml_text(const std::string& path) {
    if (!is_zip_file(path)) return read_file(path);
    for (auto& e : read_zip(path)) {
        if (e.data.find("<IPC-2581") != std::string::npos || e.data.find(":IPC-2581") != std::string::npos)
            return std::move(e.data);
    }
    throw std::runtime_error(path + ": no IPC-2581 document in the archive");
}

enum class FillMode { Fill, Hollow };

struct LineDesc {
    double width_nm = 0.0;
    std::string end = "ROUND";
};

class Importer {
    struct LayerInfo {
        std::string function;
        std::string side;
        bool negative = false;
        std::string span_from, span_to;
    };

public:
    explicit Importer(const ImportOptions& options) : tol_(options.geometry.arc_tolerance_nm) {}

    Board run(const std::string& path) {
        const std::string text = load_xml_text(path);
        pugi::xml_document doc;
        const auto res = doc.load_buffer(text.data(), text.size());
        if (!res) throw std::runtime_error(path + ": XML error: " + res.description());
        const pugi::xml_node root = doc.first_child();
        if (!is(root, "IPC-2581")) throw std::runtime_error(path + " is not an IPC-2581 document");

        board_.nets.push_back("");
        board_.source.importer = "ipc2581";
        board_.source.files.push_back(fs::path(path).filename().string());
        board_.source.arc_tolerance_nm = tol_;
        read_generator(root);

        const auto ecad = child(root, "Ecad");
        const auto header = child(ecad, "CadHeader");
        const std::string units = upper(header.attribute("units").as_string("MILLIMETER"));
        if (units == "INCH") unit_nm_ = 25.4 * kNmPerMm;
        else if (units == "MICRON") unit_nm_ = 1000.0;
        else unit_nm_ = kNmPerMm;

        read_dictionaries(child(root, "Content"));
        for (auto spec : children(header, "Spec")) specs_[spec.attribute("name").as_string()] = spec;

        const auto cad = child(ecad, "CadData");
        read_layers(cad);
        const auto step = child(cad, "Step");
        if (!step) throw std::runtime_error(path + ": no Step in CadData");
        board_.name = step.attribute("name").as_string();
        read_profile(step);
        read_components(step);
        read_layer_features(step);
        detail::assign_drill_nets(board_);
        detail::check_net_overlaps(board_);
        board_.net_source = "ipc2581";
        return std::move(board_);
    }

private:
    // ------------------------------------------------------------ header data

    void read_generator(const pugi::xml_node& root) {
        for (auto h : children(root, "HistoryRecord")) {
            const auto sp = child(child(h, "FileRevision"), "SoftwarePackage");
            if (sp) {
                board_.source.generator = std::string(sp.attribute("vendor").as_string()) + " " + sp.attribute("name").as_string() +
                                          " " + sp.attribute("revision").as_string();
                while (!board_.source.generator.empty() && board_.source.generator.back() == ' ') board_.source.generator.pop_back();
            }
        }
    }

    void read_dictionaries(const pugi::xml_node& content) {
        for (auto e : children(child(content, "DictionaryStandard"), "EntryStandard"))
            standard_[e.attribute("id").as_string()] = e.first_child();
        for (auto e : children(child(content, "DictionaryUser"), "EntryUser"))
            user_[e.attribute("id").as_string()] = e.first_child();
        for (auto e : children(child(content, "DictionaryLineDesc"), "EntryLineDesc"))
            line_desc_[e.attribute("id").as_string()] = child(e, "LineDesc");
        for (auto e : children(child(content, "DictionaryFillDesc"), "EntryFillDesc"))
            fill_desc_[e.attribute("id").as_string()] = child(e, "FillDesc");
    }

    static LayerType type_of(const std::string& function) {
        static const std::set<std::string> copper{"CONDUCTOR", "PLANE", "SIGNAL", "MIXED", "CONDFILM", "CONDFOIL"};
        if (copper.count(function)) return LayerType::Copper;
        if (function.rfind("DIEL", 0) == 0) return LayerType::Dielectric;
        if (function == "SOLDERMASK") return LayerType::SolderMask;
        if (function == "SILKSCREEN" || function == "LEGEND") return LayerType::Silkscreen;
        if (function == "SOLDERPASTE") return LayerType::SolderPaste;
        return LayerType::Other;
    }

    void read_layers(const pugi::xml_node& cad) {
        std::vector<std::string> copper_order;
        for (auto l : children(cad, "Layer")) {
            LayerInfo info;
            info.function = upper(l.attribute("layerFunction").as_string());
            info.side = upper(l.attribute("side").as_string());
            info.negative = upper(l.attribute("polarity").as_string("POSITIVE")) == "NEGATIVE";
            if (auto span = child(l, "Span")) {
                info.span_from = span.attribute("fromLayer").as_string();
                info.span_to = span.attribute("toLayer").as_string();
            }
            const std::string name = l.attribute("name").as_string();
            layers_[name] = info;
            if (type_of(info.function) == LayerType::Copper) copper_order.push_back(name);
        }

        // Stackup order when present, otherwise the Layer order.
        std::vector<std::pair<int, pugi::xml_node>> stack;
        for (auto group : children(child(cad, "Stackup"), "StackupGroup"))
            for (auto sl : children(group, "StackupLayer")) stack.emplace_back(sl.attribute("sequence").as_int(), sl);
        std::stable_sort(stack.begin(), stack.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        if (!stack.empty()) {
            board_.stackup_source = "ipc2581";
            int copper = 0;
            for (const auto& [seq, sl] : stack) {
                const std::string name = sl.attribute("layerOrGroupRef").as_string();
                auto it = layers_.find(name);
                if (it == layers_.end()) continue;
                const LayerType type = type_of(it->second.function);
                const double t = sl.attribute("thickness").as_double() * unit_nm_;
                if (type != LayerType::Copper && type != LayerType::Dielectric && !(type == LayerType::SolderMask && t > 0))
                    continue;
                StackupLayer l;
                l.name = name;
                l.type = type;
                l.thickness_nm = t;
                if (type == LayerType::Copper) {
                    l.copper_index = ++copper;
                    copper_index_[name] = copper;
                }
                apply_spec(l, child(sl, "SpecRef").attribute("id").as_string());
                board_.stackup.push_back(std::move(l));
            }
        }
        if (copper_index_.empty()) {
            for (std::size_t i = 0; i < copper_order.size(); ++i) copper_index_[copper_order[i]] = static_cast<int>(i) + 1;
            board_.stackup = detail::default_stackup(static_cast<int>(copper_order.size()), 0.0);
            for (std::size_t i = 0, k = 0; i < board_.stackup.size(); ++i)
                if (board_.stackup[i].type == LayerType::Copper) board_.stackup[i].name = copper_order[k++];
            board_.stackup_source = "default";
            board_.warnings.push_back("no stackup in the IPC-2581 file; assumed 1.6 mm FR-4 with 35 um copper");
        }
        if (copper_index_.empty()) throw std::runtime_error("no copper layers in the IPC-2581 file");
        const int n = static_cast<int>(copper_index_.size());
        board_.copper.resize(n);
        for (const auto& [name, idx] : copper_index_) {
            auto& cl = board_.copper[idx - 1];
            cl.index = idx;
            cl.name = name;
            cl.side = idx == 1 ? Side::Top : (idx == n ? Side::Bottom : Side::Inner);
        }
    }

    void apply_spec(StackupLayer& l, const std::string& id) {
        auto it = specs_.find(id);
        if (it == specs_.end()) return;
        for (auto c : it->second.children()) {
            const std::string type = upper(c.attribute("type").as_string());
            const auto prop = child(c, "Property");
            if (is(c, "General") && type == "MATERIAL") {
                const std::string text = prop.attribute("text").as_string();
                if (!text.empty()) l.material = text;
            } else if (is(c, "Dielectric") && type == "DIELECTRIC_CONSTANT") {
                l.epsilon_r = prop.attribute("value").as_double();
            } else if (is(c, "Dielectric") && type == "LOSS_TANGENT") {
                l.loss_tangent = prop.attribute("value").as_double();
            } else if (is(c, "Conductor") && type == "CONDUCTIVITY") {
                l.conductivity_s_per_m = prop.attribute("value").as_double();
            }
        }
    }

    // ------------------------------------------------------------ geometry

    double L(const pugi::xml_node& n, const char* attr) const { return n.attribute(attr).as_double() * unit_nm_; }

    Point64 point(const pugi::xml_node& n, const char* ax = "x", const char* ay = "y") const {
        return Point64(std::llround(L(n, ax)), std::llround(L(n, ay)));
    }

    Placement placement_of(const pugi::xml_node& parent) const {
        Placement t;
        if (auto x = child(parent, "Xform")) {
            t.rotation_deg = x.attribute("rotation").as_double();
            t.mirror_x = x.attribute("mirror").as_bool();
            t.scale = x.attribute("scale").as_double(1.0);
            t.offset = Point64(std::llround(L(x, "xOffset")), std::llround(L(x, "yOffset")));
        }
        if (auto loc = child(parent, "Location")) {
            const Point64 p = point(loc);
            t.offset = Point64(t.offset.x + p.x, t.offset.y + p.y);
        }
        return t;
    }

    LineDesc line_desc_of(const pugi::xml_node& shape) const {
        LineDesc d;
        pugi::xml_node ld = child(shape, "LineDesc");
        if (!ld) {
            if (auto ref = child(shape, "LineDescRef")) {
                auto it = line_desc_.find(ref.attribute("id").as_string());
                if (it != line_desc_.end()) ld = it->second;
            }
        }
        if (ld) {
            d.width_nm = L(ld, "lineWidth");
            d.end = upper(ld.attribute("lineEnd").as_string("ROUND"));
        }
        return d;
    }

    FillMode fill_of(const pugi::xml_node& shape, FillMode def) const {
        pugi::xml_node fd = child(shape, "FillDesc");
        if (!fd) {
            if (auto ref = child(shape, "FillDescRef")) {
                auto it = fill_desc_.find(ref.attribute("id").as_string());
                if (it != fill_desc_.end()) fd = it->second;
            }
        }
        if (!fd) return def;
        return upper(fd.attribute("fillProperty").as_string("FILL")) == "HOLLOW" ? FillMode::Hollow : FillMode::Fill;
    }

    Paths64 stroke_open(const Path64& path, const LineDesc& d) const {
        if (d.width_nm <= 0 || path.size() < 2) return {};
        if (d.end == "ROUND") return stroke_polyline(path, d.width_nm / 2, tol_);
        Clipper2Lib::ClipperOffset off(2.0, tol_);
        off.AddPath(path, Clipper2Lib::JoinType::Round,
                    d.end == "SQUARE" ? Clipper2Lib::EndType::Square : Clipper2Lib::EndType::Butt);
        Paths64 out;
        off.Execute(d.width_nm / 2, out);
        return out;
    }

    // Polygon/Polyline steps as a point path. Arcs of filled outlines keep their area.
    Path64 poly_path(const pugi::xml_node& n, bool filled = true) const {
        Path64 path;
        for (auto c : n.children()) {
            if (is(c, "PolyBegin")) {
                path.push_back(point(c));
            } else if (is(c, "PolyStepSegment")) {
                path.push_back(point(c));
            } else if (is(c, "PolyStepCurve")) {
                const Point64 end = point(c);
                const Point64 centre = point(c, "centerX", "centerY");
                const bool cw = c.attribute("clockwise").as_bool();
                if (path.empty()) continue;
                const Path64 arc = arc_points(path.back(), end, centre, !cw, tol_, filled);
                path.insert(path.end(), arc.begin() + 1, arc.end());
            }
        }
        return path;
    }

    static Path64 closed_ring(Path64 p) {
        if (p.size() > 1 && p.front() == p.back()) p.pop_back();
        make_positive(p);
        return p;
    }

    Paths64 rect_round(double w, double h, double r, bool ur, bool ul, bool lr, bool ll) const {
        r = std::min({r, w / 2, h / 2});
        Paths64 parts{{Point64(std::llround(-w / 2), std::llround(-h / 2 + r)), Point64(std::llround(w / 2), std::llround(-h / 2 + r)),
                       Point64(std::llround(w / 2), std::llround(h / 2 - r)), Point64(std::llround(-w / 2), std::llround(h / 2 - r))},
                      {Point64(std::llround(-w / 2 + r), std::llround(-h / 2)), Point64(std::llround(w / 2 - r), std::llround(-h / 2)),
                       Point64(std::llround(w / 2 - r), std::llround(h / 2)), Point64(std::llround(-w / 2 + r), std::llround(h / 2))}};
        const std::array<std::pair<bool, Point64>, 4> corners{
            std::make_pair(ur, Point64(std::llround(w / 2 - r), std::llround(h / 2 - r))),
            std::make_pair(ul, Point64(std::llround(-w / 2 + r), std::llround(h / 2 - r))),
            std::make_pair(lr, Point64(std::llround(w / 2 - r), std::llround(-h / 2 + r))),
            std::make_pair(ll, Point64(std::llround(-w / 2 + r), std::llround(-h / 2 + r)))};
        const std::array<Point64, 4> square{Point64(std::llround(w / 2), std::llround(h / 2)), Point64(std::llround(-w / 2), std::llround(h / 2)),
                                            Point64(std::llround(w / 2), std::llround(-h / 2)), Point64(std::llround(-w / 2), std::llround(-h / 2))};
        for (std::size_t i = 0; i < 4; ++i) {
            if (corners[i].first && r > 0) {
                parts.push_back(circle_path(corners[i].second, r, tol_));
            } else {
                // Square corner: the quadrant between the corner centre and the corner.
                const Point64 c = corners[i].second, s = square[i];
                Path64 q{c, Point64(s.x, c.y), s, Point64(c.x, s.y)};
                make_positive(q);
                parts.push_back(q);
            }
        }
        return to_paths(union_polygons(parts));
    }

    // Filled area of a shape element, before placement.
    Paths64 shape(const pugi::xml_node& n) {
        auto D = [&](const char* a) { return L(n, a); };
        auto hollow = [&](const Path64& ring) {
            const LineDesc d = line_desc_of(n);
            return d.width_nm > 0 ? stroke_ring(ring, d.width_nm / 2, tol_) : Paths64{};
        };
        auto filled_or_hollow = [&](Path64 ring) -> Paths64 {
            make_positive(ring);
            return fill_of(n, FillMode::Fill) == FillMode::Fill ? Paths64{ring} : hollow(ring);
        };
        const char* name = local_name(n);
        if (!std::strcmp(name, "Circle")) return filled_or_hollow(circle_path(Point64(0, 0), D("diameter") / 2, tol_));
        if (!std::strcmp(name, "RectCenter")) {
            const double w = D("width"), h = D("height");
            return filled_or_hollow({Point64(std::llround(-w / 2), std::llround(-h / 2)), Point64(std::llround(w / 2), std::llround(-h / 2)),
                                     Point64(std::llround(w / 2), std::llround(h / 2)), Point64(std::llround(-w / 2), std::llround(h / 2))});
        }
        if (!std::strcmp(name, "RectRound")) {
            return rect_round(D("width"), D("height"), D("radius"), n.attribute("upperRight").as_bool(), n.attribute("upperLeft").as_bool(),
                              n.attribute("lowerRight").as_bool(), n.attribute("lowerLeft").as_bool());
        }
        if (!std::strcmp(name, "RectCham")) {
            const double w = D("width"), h = D("height"), c = std::min({D("chamfer"), w / 2, h / 2});
            const bool ur = n.attribute("upperRight").as_bool(), ul = n.attribute("upperLeft").as_bool();
            const bool lr = n.attribute("lowerRight").as_bool(), ll = n.attribute("lowerLeft").as_bool();
            Path64 p;
            auto P = [&](double x, double y) { p.emplace_back(std::llround(x), std::llround(y)); };
            if (ll) { P(-w / 2, -h / 2 + c); P(-w / 2 + c, -h / 2); } else P(-w / 2, -h / 2);
            if (lr) { P(w / 2 - c, -h / 2); P(w / 2, -h / 2 + c); } else P(w / 2, -h / 2);
            if (ur) { P(w / 2, h / 2 - c); P(w / 2 - c, h / 2); } else P(w / 2, h / 2);
            if (ul) { P(-w / 2 + c, h / 2); P(-w / 2, h / 2 - c); } else P(-w / 2, h / 2);
            return filled_or_hollow(p);
        }
        if (!std::strcmp(name, "Oval")) {
            const double w = D("width"), h = D("height");
            if (w == h) return filled_or_hollow(circle_path(Point64(0, 0), w / 2, tol_));
            const bool horizontal = w > h;
            const double r = std::min(w, h) / 2, d = (std::max(w, h) - std::min(w, h)) / 2;
            const Point64 a = horizontal ? Point64(std::llround(-d), 0LL) : Point64(0LL, std::llround(-d));
            const Point64 b = horizontal ? Point64(std::llround(d), 0LL) : Point64(0LL, std::llround(d));
            return filled_or_hollow(capsule(a, b, r, tol_));
        }
        if (!std::strcmp(name, "Ellipse")) {
            const double a = D("width") / 2, b = D("height") / 2;
            const int k = circle_segments(std::max(a, b), tol_);
            Path64 p;
            for (int i = 0; i < k; ++i) {
                const double t = 2 * kPi * i / k;
                p.emplace_back(std::llround(a * std::cos(t)), std::llround(b * std::sin(t)));
            }
            return filled_or_hollow(p);
        }
        if (!std::strcmp(name, "Diamond")) {
            const double w = D("width") / 2, h = D("height") / 2;
            return filled_or_hollow({Point64(std::llround(w), 0LL), Point64(0LL, std::llround(h)), Point64(std::llround(-w), 0LL), Point64(0LL, std::llround(-h))});
        }
        if (!std::strcmp(name, "Octagon") || !std::strcmp(name, "Hexagon")) {
            const int k = name[0] == 'O' ? 8 : 6;
            // "length": distance across flats.
            const double r = D("length") / 2 / std::cos(kPi / k);
            Path64 p;
            for (int i = 0; i < k; ++i) {
                const double t = kPi / k + 2 * kPi * i / k;
                p.emplace_back(std::llround(r * std::cos(t)), std::llround(r * std::sin(t)));
            }
            return filled_or_hollow(p);
        }
        if (!std::strcmp(name, "Triangle")) {
            const double b = D("base"), h = D("height");
            return filled_or_hollow({Point64(std::llround(-b / 2), std::llround(-h / 2)), Point64(std::llround(b / 2), std::llround(-h / 2)),
                                     Point64(0LL, std::llround(h / 2))});
        }
        if (!std::strcmp(name, "Donut")) {
            const double od = D("outerDiameter"), id = D("innerDiameter");
            Paths64 ring{circle_path(Point64(0, 0), od / 2, tol_)};
            if (id > 0) {
                Path64 hole = circle_path(Point64(0, 0), id / 2, tol_);
                std::reverse(hole.begin(), hole.end());
                ring.push_back(hole);
            }
            if (upper(n.attribute("shape").as_string("ROUND")) != "ROUND") warn("Donut of shape " + std::string(n.attribute("shape").as_string()) + " drawn round");
            return ring;
        }
        if (!std::strcmp(name, "Contour")) {
            Paths64 outer;
            Paths64 cut;
            for (auto c : n.children()) {
                if (is(c, "Polygon")) outer.push_back(closed_ring(poly_path(c)));
                else if (is(c, "Cutout")) cut.push_back(closed_ring(poly_path(c)));
            }
            if (cut.empty()) return outer;
            return to_paths(difference(union_polygons(outer), union_polygons(cut)));
        }
        if (!std::strcmp(name, "Polygon")) {
            const bool filled = fill_of(n, FillMode::Fill) == FillMode::Fill;
            return filled_or_hollow(closed_ring(poly_path(n, filled)));
        }
        if (!std::strcmp(name, "Polyline")) return stroke_open(poly_path(n, false), line_desc_of(n));
        if (!std::strcmp(name, "Line")) {
            return stroke_open({point(n, "startX", "startY"), point(n, "endX", "endY")}, line_desc_of(n));
        }
        if (!std::strcmp(name, "Arc")) {
            const Point64 s = point(n, "startX", "startY"), e = point(n, "endX", "endY"), c = point(n, "centerX", "centerY");
            const Path64 pts = arc_points(s, e, c, !n.attribute("clockwise").as_bool(), tol_);
            const LineDesc d = line_desc_of(n);
            if (s == e) {
                Path64 ring(pts.begin(), pts.end() - 1);
                return d.width_nm > 0 ? stroke_ring(ring, d.width_nm / 2, tol_) : Paths64{};
            }
            return stroke_open(pts, d);
        }
        if (!std::strcmp(name, "UserSpecial")) {
            Paths64 all;
            for (auto c : n.children()) {
                const Paths64 s = shape(c);
                all.insert(all.end(), s.begin(), s.end());
            }
            return all;
        }
        if (!std::strcmp(name, "StandardPrimitiveRef")) {
            auto it = standard_.find(n.attribute("id").as_string());
            if (it == standard_.end()) {
                warn(std::string("unknown standard primitive ") + n.attribute("id").as_string());
                return {};
            }
            return shape(it->second);
        }
        if (!std::strcmp(name, "UserPrimitiveRef")) {
            auto it = user_.find(n.attribute("id").as_string());
            if (it == user_.end()) {
                warn(std::string("unknown user primitive ") + n.attribute("id").as_string());
                return {};
            }
            return shape(it->second);
        }
        static const std::set<std::string> ignore{"Location", "Xform", "LineDesc", "LineDescRef", "FillDesc", "FillDescRef",
                                                  "PinRef", "ColorRef", "Color", "NonstandardAttribute", "Text"};
        if (!ignore.count(name)) warn(std::string("unsupported shape ") + name);
        return {};
    }

    // ------------------------------------------------------------ board content

    void read_profile(const pugi::xml_node& step) {
        const auto profile = child(step, "Profile");
        Paths64 outer, cut;
        for (auto c : profile.children()) {
            if (is(c, "Polygon")) outer.push_back(closed_ring(poly_path(c)));
            else if (is(c, "Cutout")) cut.push_back(closed_ring(poly_path(c)));
        }
        if (outer.empty()) {
            warn("no Profile polygon; the outline is the bounding box of the copper");
            return;
        }
        board_.outline = difference(union_polygons(outer), union_polygons(cut));
    }

    void read_components(const pugi::xml_node& step) {
        const int n = static_cast<int>(board_.copper.size());
        for (auto c : children(step, "Component")) {
            Component comp;
            comp.refdes = c.attribute("refDes").as_string();
            comp.package = c.attribute("packageRef").as_string();
            comp.value = c.attribute("part").as_string();
            const std::string layer = c.attribute("layerRef").as_string();
            auto it = copper_index_.find(layer);
            comp.side = it != copper_index_.end() && it->second == n && n > 1 ? Side::Bottom : Side::Top;
            const Placement t = placement_of(c);
            comp.at = t.offset;
            // A mirrored placement's rotation is in the mirrored frame (mirror about y, then
            // rotate); the model's is counter-clockwise as seen from the top.
            comp.rotation_deg = t.mirror_x ? std::fmod(540.0 - t.rotation_deg, 360.0) : t.rotation_deg;
            comp.position_source = "ipc2581";
            board_.components.push_back(std::move(comp));
        }
    }

    // Features of one Set, appended to out.
    // The centre lines of round-ended lines, arcs and polylines, as strokes without a
    // shape of their own (the shape is drawn by the feature itself): they become tracks.
    void stroke_segments(const pugi::xml_node& n, const Placement& t, const Feature& base, std::vector<Feature>& out) const {
        if (is(n, "UserSpecial")) {
            for (auto c : n.children()) stroke_segments(c, t, base, out);
            return;
        }
        const bool line = is(n, "Line"), arc = is(n, "Arc"), poly = is(n, "Polyline");
        if (!line && !arc && !poly) return;
        const LineDesc d = line_desc_of(n);
        if (!(d.width_nm > 0) || upper(d.end) != "ROUND") return;
        const bool mirrored = t.mirror_x;
        auto add = [&](Point64 a, Point64 b, std::optional<Point64> centre, bool ccw) {
            if (a == b) return;
            Feature f;
            f.kind = FeatureKind::Stroke;
            f.dark = base.dark;
            f.has_net = base.has_net;
            f.net = base.net;
            f.function = base.function;
            f.at = transform(a, t);
            f.end = transform(b, t);
            if (centre) {
                f.arc_centre = transform(*centre, t);
                f.arc_ccw = mirrored ? !ccw : ccw;
            }
            f.stroke_width_nm = d.width_nm * (t.scale > 0 ? t.scale : 1.0);
            out.push_back(std::move(f));
        };
        if (line) add(point(n, "startX", "startY"), point(n, "endX", "endY"), std::nullopt, false);
        if (arc) add(point(n, "startX", "startY"), point(n, "endX", "endY"), point(n, "centerX", "centerY"), !n.attribute("clockwise").as_bool());
        if (poly) {
            std::optional<Point64> prev;
            for (auto c : n.children()) {
                if (is(c, "PolyBegin")) prev = point(c);
                else if (is(c, "PolyStepSegment") && prev) {
                    add(*prev, point(c), std::nullopt, false);
                    prev = point(c);
                } else if (is(c, "PolyStepCurve") && prev) {
                    add(*prev, point(c), point(c, "centerX", "centerY"), !c.attribute("clockwise").as_bool());
                    prev = point(c);
                }
            }
        }
    }

    void set_features(const pugi::xml_node& set, std::vector<Feature>& out) {
        const std::string net = set.attribute("net").as_string();
        const bool dark = upper(set.attribute("polarity").as_string("POSITIVE")) != "NEGATIVE";
        const std::string usage = set.attribute("padUsage").as_string();
        for (auto c : set.children()) {
            if (is(c, "Features")) {
                const Placement t = placement_of(c);
                for (auto s : c.children()) {
                    if (is(s, "Location") || is(s, "Xform")) continue;
                    Feature f;
                    f.kind = FeatureKind::Region;
                    f.dark = dark;
                    f.has_net = true;
                    f.net = net;
                    f.function = usage;
                    f.shape = transform(shape(s), t);
                    f.at = t.offset;
                    out.push_back(std::move(f));
                    const Feature base = out.back();
                    stroke_segments(s, t, base, out);
                }
            } else if (is(c, "Pad")) {
                const Placement t = placement_of(c);
                Feature f;
                f.kind = FeatureKind::Flash;
                f.dark = dark;
                f.has_net = true;
                f.net = net;
                f.function = usage.empty() ? "pad" : usage;
                f.at = t.offset;
                for (auto s : c.children()) {
                    if (is(s, "PinRef")) {
                        f.pad = true;
                        f.component = s.attribute("componentRef").as_string();
                        f.pin = s.attribute("pin").as_string();
                        continue;
                    }
                    const Paths64 p = transform(shape(s), t);
                    f.shape.insert(f.shape.end(), p.begin(), p.end());
                }
                out.push_back(std::move(f));
            }
        }
    }

    void read_drills(const pugi::xml_node& set, const LayerInfo& info) {
        const int n = static_cast<int>(board_.copper.size());
        int from = 1, to = n;
        if (!info.span_from.empty()) {
            auto a = copper_index_.find(info.span_from), b = copper_index_.find(info.span_to);
            if (a != copper_index_.end() && b != copper_index_.end()) {
                from = std::min(a->second, b->second);
                to = std::max(a->second, b->second);
            }
        }
        const std::string net = set.attribute("net").as_string();
        for (auto c : set.children()) {
            Drill d;
            d.from_layer = from;
            d.to_layer = to;
            d.net = board_.net_id(net);
            std::string plating;
            if (is(c, "Hole")) {
                d.at = point(c);
                d.diameter_nm = L(c, "diameter");
                plating = upper(c.attribute("platingStatus").as_string("PLATED"));
            } else if (is(c, "SlotCavity")) {
                plating = upper(c.attribute("platingStatus").as_string("PLATED"));
                const Placement t = placement_of(c);
                pugi::xml_node s;
                for (auto k : c.children())
                    if (!is(k, "Location") && !is(k, "Xform")) s = k;
                double w = L(s, "width"), h = L(s, "height");
                if (!is(s, "Oval")) warn(std::string("slot of shape ") + local_name(s) + " drawn as an oval");
                const double dia = std::min(w, h), half = (std::max(w, h) - dia) / 2;
                const Point64 a = w >= h ? Point64(std::llround(-half), 0LL) : Point64(0LL, std::llround(-half));
                const Point64 b = w >= h ? Point64(std::llround(half), 0LL) : Point64(0LL, std::llround(half));
                d.at = transform(a, t);
                d.slot_end = transform(b, t);
                d.diameter_nm = dia;
            } else {
                continue;
            }
            d.plated = plating != "NONPLATED";
            // IPC-2581 does not state the function of non-plated holes.
            d.function = plating == "VIA" ? DrillFunction::Via : (d.plated ? DrillFunction::Component : DrillFunction::Unknown);
            if (!d.plated) d.net = 0;
            board_.drills.push_back(d);
        }
    }

    void read_layer_features(const pugi::xml_node& step) {
        const int n = static_cast<int>(board_.copper.size());
        std::vector<std::vector<Feature>> copper(n);
        std::map<Side, std::vector<Feature>> masks;
        for (auto lf : children(step, "LayerFeature")) {
            const std::string name = lf.attribute("layerRef").as_string();
            auto it = layers_.find(name);
            if (it == layers_.end()) continue;
            const LayerInfo& info = it->second;
            const LayerType type = type_of(info.function);
            if (info.function == "DRILL" || info.function == "ROUT") {
                for (auto set : children(lf, "Set")) read_drills(set, info);
            } else if (type == LayerType::Copper) {
                if (info.negative) warn("negative copper layer " + name + " is not supported");
                const int idx = copper_index_.at(name);
                for (auto set : children(lf, "Set")) set_features(set, copper[idx - 1]);
            } else if (type == LayerType::SolderMask) {
                const Side side = info.side == "BOTTOM" ? Side::Bottom : Side::Top;
                for (auto set : children(lf, "Set")) set_features(set, masks[side]);
            }
        }

        std::vector<std::future<std::map<std::string, PolygonSet>>> rendered;
        for (int i = 0; i < n; ++i)
            rendered.push_back(std::async(std::launch::async, [&copper, i]() { return render_by_net(copper[i]); }));
        for (int i = 0; i < n; ++i) detail::set_layer_copper(board_, board_.copper[i], rendered[i].get());
        {
            std::vector<const std::vector<Feature>*> layers;
            for (int i = 0; i < n; ++i) layers.push_back(&copper[i]);
            detail::collect_tracks(board_, layers);
        }

        for (int i = 0; i < n; ++i) {
            for (const auto& f : copper[i]) {
                if (!f.is_pad() || !f.dark || f.shape.empty()) continue;
                Pad p;
                p.component = f.component;
                p.pin = f.pin;
                p.function = f.function;
                p.net = board_.net_id(f.net);
                p.layer = i + 1;
                p.at = f.at;
                p.shape = union_polygons(f.shape);
                board_.pads.push_back(std::move(p));
            }
        }
        for (auto& [side, features] : masks) board_.masks.push_back({side, render_all(features)});

        if (board_.outline.empty()) {
            Rect64 r(INT64_MAX, INT64_MAX, INT64_MIN, INT64_MIN);
            for (const auto& l : board_.copper)
                for (const auto& nc : l.nets) {
                    const Rect64 b = bounds(nc.polygons);
                    r.left = std::min(r.left, b.left);
                    r.top = std::min(r.top, b.top);
                    r.right = std::max(r.right, b.right);
                    r.bottom = std::max(r.bottom, b.bottom);
                }
            board_.outline = {Polygon{{Point64(r.left, r.top), Point64(r.right, r.top), Point64(r.right, r.bottom), Point64(r.left, r.bottom)}, {}}};
        }
    }

    void warn(const std::string& msg) {
        if (warned_.insert(msg).second && board_.warnings.size() < 100) board_.warnings.push_back(msg);
    }

    double tol_;
    double unit_nm_ = kNmPerMm;
    Board board_;
    std::map<std::string, pugi::xml_node> standard_, user_, line_desc_, fill_desc_, specs_;
    std::map<std::string, LayerInfo> layers_;
    std::map<std::string, int> copper_index_;
    std::set<std::string> warned_;
};

}  // namespace

Board import_ipc2581(const std::string& path, const ImportOptions& options) { return Importer(options).run(path); }

bool looks_like_ipc2581(const std::string& path) {
    try {
        if (is_zip_file(path)) {
            for (const auto& e : read_zip(path))
                if (e.data.find("<IPC-2581") != std::string::npos) return true;
            return false;
        }
        std::ifstream in(path, std::ios::binary);
        std::string head(4096, '\0');
        in.read(head.data(), static_cast<std::streamsize>(head.size()));
        head.resize(static_cast<std::size_t>(in.gcount()));
        return head.find("<IPC-2581") != std::string::npos || head.find(":IPC-2581") != std::string::npos;
    } catch (const std::exception&) {
        return false;
    }
}

}  // namespace kairoq
