// Importer for a set of fabrication files.
#include "kairoq/fabrication.hpp"

#include "build.hpp"
#include "kairoq/archive.hpp"
#include "kairoq/drill.hpp"
#include "kairoq/gerber.hpp"

#include <algorithm>
#include <cstdio>
#include <cctype>
#include <filesystem>
#include <future>
#include <map>
#include <regex>
#include <tuple>

namespace kairoq {

namespace fs = std::filesystem;

namespace {

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

struct InputFile {
    std::string name;  // file name, used for roles and reporting
    std::string text;
};

std::vector<InputFile> collect_inputs(const std::vector<std::string>& inputs) {
    std::vector<InputFile> out;
    auto add_path = [&](const fs::path& p) {
        if (is_zip_file(p.string())) {
            for (auto& e : read_zip(p.string())) out.push_back({fs::path(e.name).filename().string(), std::move(e.data)});
        } else {
            out.push_back({p.filename().string(), read_file(p.string())});
        }
    };
    for (const auto& in : inputs) {
        const fs::path p(in);
        if (fs::is_directory(p)) {
            std::vector<fs::path> files;
            for (const auto& e : fs::directory_iterator(p))
                if (e.is_regular_file()) files.push_back(e.path());
            std::sort(files.begin(), files.end());
            for (const auto& f : files) add_path(f);
        } else {
            add_path(p);
        }
    }
    return out;
}

enum class Role { Copper, Mask, Profile, Drill, Other };

struct GerberRole {
    Role role = Role::Other;
    int copper_index = 0;  // 0: resolve later (bottom layer from an extension)
    bool bottom = false;
};

GerberRole role_from_function(const std::vector<std::string>& ff) {
    GerberRole r;
    if (ff.empty()) return r;
    const std::string kind = lower(ff[0]);
    if (kind == "copper" && ff.size() >= 2 && ff[1].size() > 1 && (ff[1][0] == 'L' || ff[1][0] == 'l')) {
        r.role = Role::Copper;
        r.copper_index = std::stoi(ff[1].substr(1));
    } else if (kind == "soldermask") {
        r.role = Role::Mask;
        r.bottom = ff.size() > 1 && lower(ff[1]).rfind("bot", 0) == 0;
    } else if (kind == "profile") {
        r.role = Role::Profile;
    } else if (kind == "plated" || kind == "nonplated") {
        r.role = Role::Drill;
    }
    return r;
}

// Protel-style extensions, for Gerber files without X2 attributes.
GerberRole role_from_extension(const std::string& name) {
    GerberRole r;
    const std::string ext = lower(fs::path(name).extension().string());
    static const std::regex inner(R"(\.g(\d+))");
    std::smatch m;
    if (ext == ".gtl") {
        r.role = Role::Copper;
        r.copper_index = 1;
    } else if (ext == ".gbl") {
        r.role = Role::Copper;
        r.bottom = true;
    } else if (std::regex_match(ext, m, inner)) {
        r.role = Role::Copper;
        r.copper_index = std::stoi(m[1]) + 1;
    } else if (ext == ".gts" || ext == ".gbs") {
        r.role = Role::Mask;
        r.bottom = ext == ".gbs";
    } else if (ext == ".gko" || ext == ".gm1" || ext == ".gml" || ext == ".gm") {
        r.role = Role::Profile;
    }
    return r;
}

DrillFunction drill_function(const std::string& x2) {
    const std::string f = lower(x2);
    if (f.find("via") != std::string::npos) return DrillFunction::Via;
    if (f.find("component") != std::string::npos || f.find("castellated") != std::string::npos ||
        f.find("pressfit") != std::string::npos)
        return DrillFunction::Component;
    if (f.find("mechanical") != std::string::npos || f.find("slot") != std::string::npos)
        return DrillFunction::Mechanical;
    return DrillFunction::Unknown;
}

struct ParsedGerber {
    std::string name;
    GerberRole role;
    GerberImage image;
};

}  // namespace

Board import_fabrication(const std::vector<std::string>& inputs, const ImportOptions& options) {
    const auto files = collect_inputs(inputs);
    Board board;
    board.nets.push_back("");
    board.source.importer = "fabrication";
    board.source.arc_tolerance_nm = options.geometry.arc_tolerance_nm;

    std::optional<JobFile> job;
    std::vector<std::pair<std::string, DrillFile>> excellon;
    std::optional<Netlist> netlist;
    std::vector<PlacementRow> placement;
    std::optional<std::vector<StackupLayer>> yaml_stackup;
    std::vector<const InputFile*> gerbers;

    for (const auto& f : files) {
        const std::string ext = lower(fs::path(f.name).extension().string());
        if (ext == ".gbrjob") {
            if (job) board.warnings.push_back("more than one job file; using the first");
            else job = parse_job_file(f.text);
        } else if (ext == ".yaml" || ext == ".yml") {
            if (lower(f.name).find("stackup") != std::string::npos) yaml_stackup = parse_stackup_yaml(f.text);
            else continue;
        } else if (looks_like_gerber(f.text)) {
            gerbers.push_back(&f);
        } else if (looks_like_excellon(f.text)) {
            excellon.emplace_back(f.name, parse_excellon(f.text));
        } else if (looks_like_ipc356(f.text)) {
            netlist = parse_ipc356(f.text);
        } else if (ext == ".csv" && looks_like_placement_csv(f.text)) {
            placement = parse_placement_csv(f.text);
        } else {
            board.warnings.push_back("ignored file " + f.name);
            continue;
        }
        board.source.files.push_back(f.name);
    }
    if (!options.stackup_yaml_path.empty()) {
        yaml_stackup = parse_stackup_yaml(read_file(options.stackup_yaml_path));
        board.source.files.push_back(fs::path(options.stackup_yaml_path).filename().string());
    }

    // Parse the Gerber files in parallel; each is independent.
    std::vector<std::future<ParsedGerber>> futures;
    for (const InputFile* f : gerbers) {
        futures.push_back(std::async(std::launch::async, [f, &options]() {
            ParsedGerber pg;
            pg.name = f->name;
            pg.image = parse_gerber(f->text, options.geometry);
            return pg;
        }));
    }
    std::vector<ParsedGerber> parsed;
    for (auto& fu : futures) parsed.push_back(fu.get());

    std::map<std::string, std::vector<std::string>> job_functions;
    if (job) {
        for (const auto& e : job->files) job_functions[e.path] = e.function;
        board.source.generator = job->generator;
        board.name = job->project;
    }
    int max_index = 0;
    for (auto& pg : parsed) {
        for (const auto& w : pg.image.warnings) board.warnings.push_back(pg.name + ": " + w);
        auto ff = pg.image.file_function();
        if (ff.empty())
            if (auto it = job_functions.find(pg.name); it != job_functions.end()) ff = it->second;
        pg.role = role_from_function(ff);
        if (pg.role.role == Role::Other && ff.empty()) pg.role = role_from_extension(pg.name);
        if (pg.role.role == Role::Copper) max_index = std::max(max_index, pg.role.copper_index);
        if (board.source.generator.empty()) {
            auto it = pg.image.file_attributes.find("GenerationSoftware");
            if (it != pg.image.file_attributes.end()) {
                std::string g;
                for (const auto& v : it->second) g += (g.empty() ? "" : " ") + v;
                board.source.generator = g;
            }
        }
    }
    int layer_count = job && job->layer_count > 0 ? job->layer_count : max_index;
    for (auto& pg : parsed)
        if (pg.role.role == Role::Copper && pg.role.copper_index == 0) pg.role.copper_index = std::max(layer_count, max_index + 1);
    for (const auto& pg : parsed)
        if (pg.role.role == Role::Copper) layer_count = std::max(layer_count, pg.role.copper_index);
    if (layer_count == 0) throw std::runtime_error("no copper layers found among the inputs");
    if (board.name.empty() && !files.empty()) board.name = fs::path(files.front().name).stem().string();

    // Stackup.
    if (yaml_stackup) {
        board.stackup = *yaml_stackup;
        board.stackup_source = "yaml";
    } else if (job && !job->stackup.empty()) {
        board.stackup = job->stackup;
        board.stackup_source = "gbrjob";
    } else {
        board.stackup = detail::default_stackup(layer_count, job ? job->board_thickness_nm : 0.0);
        board.stackup_source = "default";
        board.warnings.push_back("no stackup in the inputs; assumed 1.6 mm FR-4 with 35 um copper");
    }
    int stackup_copper = 0;
    for (const auto& l : board.stackup)
        if (l.type == LayerType::Copper) ++stackup_copper;
    if (stackup_copper != layer_count)
        throw std::runtime_error("the stackup has " + std::to_string(stackup_copper) + " copper layers, the Gerber files " +
                                 std::to_string(layer_count));

    board.copper.resize(layer_count);
    for (int i = 0; i < layer_count; ++i) {
        auto& cl = board.copper[i];
        cl.index = i + 1;
        cl.side = i == 0 ? Side::Top : (i == layer_count - 1 ? Side::Bottom : Side::Inner);
        const auto* s = board.stackup_for_copper(i + 1);
        cl.name = s && !s->name.empty() ? s->name : "L" + std::to_string(i + 1);
    }

    // Copper, pads, masks, outline.
    bool all_x2_nets = true;
    std::vector<bool> have_layer(layer_count, false);
    std::vector<PolygonSet> plain_copper(layer_count);
    std::vector<std::future<std::map<std::string, PolygonSet>>> rendered(layer_count);
    std::vector<const ParsedGerber*> layer_source(layer_count, nullptr);
    std::vector<Feature> profile;
    std::vector<std::pair<std::string, DrillFile>> gerber_drills;
    for (const auto& pg : parsed) {
        switch (pg.role.role) {
        case Role::Copper: {
            const int idx = pg.role.copper_index - 1;
            if (have_layer[idx]) {
                board.warnings.push_back("more than one file for copper layer " + std::to_string(idx + 1) + "; using " +
                                         layer_source[idx]->name);
                continue;
            }
            have_layer[idx] = true;
            layer_source[idx] = &pg;
            // A layer without objects says nothing about net attributes.
            if (!pg.image.has_net_attributes && !pg.image.features.empty()) all_x2_nets = false;
            rendered[idx] = std::async(std::launch::async, [&pg]() { return render_by_net(pg.image.features); });
            break;
        }
        case Role::Mask: {
            MaskLayer m;
            m.side = pg.role.bottom ? Side::Bottom : Side::Top;
            m.openings = render_all(pg.image.features);
            board.masks.push_back(std::move(m));
            break;
        }
        case Role::Profile:
            profile.insert(profile.end(), pg.image.features.begin(), pg.image.features.end());
            break;
        case Role::Drill:
            gerber_drills.emplace_back(pg.name, drill_from_gerber(pg.image));
            break;
        case Role::Other:
            break;
        }
    }
    for (int i = 0; i < layer_count; ++i)
        if (!have_layer[i]) throw std::runtime_error("no Gerber file for copper layer " + std::to_string(i + 1));

    std::vector<std::map<std::string, PolygonSet>> by_net(layer_count);
    for (int i = 0; i < layer_count; ++i) by_net[i] = rendered[i].get();

    if (profile.empty()) {
        board.warnings.push_back("no board profile; the outline is the bounding box of the copper");
    } else {
        auto outline = outline_from_features(profile, options.geometry.arc_tolerance_nm, options.outline_join_nm);
        board.outline = std::move(outline.polygons);
        if (outline.open_chains)
            board.warnings.push_back(std::to_string(outline.open_chains) + " open chains in the board profile were ignored");
    }

    // Drills: Excellon preferred; Gerber drill files otherwise.
    auto& drill_files = excellon.empty() ? gerber_drills : excellon;
    if (!excellon.empty() && !gerber_drills.empty())
        board.warnings.push_back("both Excellon and Gerber drill files given; using the Excellon files");
    for (auto& [name, df] : drill_files) {
        for (const auto& w : df.warnings) board.warnings.push_back(name + ": " + w);
        bool plated = true;
        if (df.plated) plated = *df.plated;
        else if (lower(name).find("npth") != std::string::npos) plated = false;
        else if (lower(name).find("pth") == std::string::npos)
            board.warnings.push_back(name + ": plating not stated; assumed plated");
        const int from = df.from_layer > 0 ? std::min(df.from_layer, df.to_layer) : 1;
        const int to = df.to_layer > 0 ? std::max(df.from_layer, df.to_layer) : layer_count;
        for (const auto& h : df.hits) {
            const auto& tool = df.tools[h.tool];
            Drill d;
            d.at = h.at;
            d.slot_end = h.slot_end;
            d.diameter_nm = tool.diameter_nm;
            d.plated = plated;
            d.function = drill_function(tool.function);
            if (d.function == DrillFunction::Unknown && !plated) d.function = DrillFunction::Mechanical;
            d.from_layer = from;
            d.to_layer = to;
            board.drills.push_back(d);
        }
    }
    // Drill and placement files may be written relative to an auxiliary origin, usually the
    // board's lower-left corner: when none of their points lies on the board and all do
    // relative to that corner, they are moved there.
    auto register_points = [&](std::vector<Point64*> pts, const std::string& what) {
        if (board.outline.empty() || pts.empty()) return;
        auto inside = [&](Point64 d) {
            std::size_t n = 0;
            for (const Point64* p : pts) n += contains(board.outline, Point64(p->x + d.x, p->y + d.y));
            return n;
        };
        if (inside({0, 0}) > 0) return;
        const Rect64 r = bounds(board.outline);
        const Point64 corner(r.left, r.top);  // lowest x and y
        if (inside(corner) != pts.size()) return;
        for (Point64* p : pts) *p = Point64(p->x + corner.x, p->y + corner.y);
        char line[200];
        std::snprintf(line, sizeof line, "%s coordinates are relative to the board's lower-left corner; moved by (%.3f, %.3f) mm", what.c_str(),
                      corner.x * 1e-6, corner.y * 1e-6);
        board.warnings.push_back(line);
    };
    {
        std::vector<Point64*> pts;
        for (auto& d : board.drills) {
            pts.push_back(&d.at);
            if (d.slot_end) pts.push_back(&*d.slot_end);
        }
        register_points(pts, "drill");
    }

    if (all_x2_nets) {
        board.net_source = "x2";
        for (int i = 0; i < layer_count; ++i) detail::set_layer_copper(board, board.copper[i], by_net[i]);
        // Pins of outer-layer pads by component and position: KiCad marks inner-layer
        // pads of through-hole pins with .C only.
        std::map<std::tuple<std::string, std::int64_t, std::int64_t>, std::string> outer_pins;
        for (int i : {0, layer_count - 1})
            for (const auto& f : layer_source[i]->image.features)
                if (f.is_pad() && f.kind == FeatureKind::Flash) outer_pins[{f.component, f.at.x, f.at.y}] = f.pin;
        for (int i = 0; i < layer_count; ++i) {
            for (const auto& f : layer_source[i]->image.features) {
                if (!f.dark || f.shape.empty()) continue;
                std::string pin = f.pin;
                if (!f.is_pad()) {
                    const bool inner_pad = f.kind == FeatureKind::Flash && !f.component.empty() &&
                                           f.function.find("Pad") != std::string::npos;
                    if (!inner_pad) continue;
                    auto it = outer_pins.find({f.component, f.at.x, f.at.y});
                    if (it != outer_pins.end()) pin = it->second;
                }
                Pad p;
                p.component = f.component;
                p.pin = pin;
                p.function = f.function;
                p.net = board.net_id(f.net == "N/C" ? std::string() : f.net);
                p.layer = i + 1;
                p.shape = union_polygons(f.shape);
                if (f.kind == FeatureKind::Flash) {
                    p.at = f.at;
                } else {
                    const Rect64 b = bounds(p.shape);
                    p.at = Point64((b.left + b.right) / 2, (b.top + b.bottom) / 2);
                }
                board.pads.push_back(std::move(p));
            }
        }
        detail::assign_drill_nets(board);
        if (netlist) detail::cross_check_netlist(board, *netlist);
    } else {
        board.net_source = netlist ? "ipc356" : "geometric";
        board.warnings.push_back("copper without X2 net attributes; nets derived from connectivity");
        for (int i = 0; i < layer_count; ++i) {
            Paths64 all;
            for (const auto& [net, polys] : by_net[i]) {
                const auto p = to_paths(polys);
                all.insert(all.end(), p.begin(), p.end());
            }
            plain_copper[i] = union_polygons(all);
        }
        if (netlist) {
            std::string note;
            if (!detail::register_netlist(*netlist, board, note) || !note.empty()) board.warnings.push_back(note);
        }
        detail::geometric_nets(board, plain_copper, netlist ? &*netlist : nullptr);
    }
    {
        std::vector<const std::vector<Feature>*> layers;
        for (int i = 0; i < layer_count; ++i) layers.push_back(&layer_source[i]->image.features);
        detail::collect_tracks(board, layers);
    }

    if (board.outline.empty()) {
        Rect64 r(INT64_MAX, INT64_MAX, INT64_MIN, INT64_MIN);
        for (const auto& l : board.copper)
            for (const auto& nc : l.nets) {
                const Rect64 b = bounds(nc.polygons);
                r.left = std::min(r.left, b.left);
                r.top = std::min(r.top, b.top);
                r.right = std::max(r.right, b.right);
                r.bottom = std::max(r.bottom, b.bottom);
            }
        board.outline = {Polygon{{Point64(r.left, r.top), Point64(r.right, r.top), Point64(r.right, r.bottom),
                                  Point64(r.left, r.bottom)},
                                 {}}};
    }

    if (!placement.empty()) {
        std::vector<Point64*> pts;
        for (auto& row : placement) pts.push_back(&row.at);
        register_points(pts, "placement");
        for (const auto& row : placement) {
            Component c;
            c.refdes = row.refdes;
            c.value = row.value;
            c.package = row.package;
            c.at = row.at;
            c.rotation_deg = row.rotation_deg;
            c.side = row.side;
            c.position_source = "placement";
            board.components.push_back(std::move(c));
        }
    } else {
        detail::components_from_pads(board);
    }
    detail::check_net_overlaps(board);
    return board;
}

}  // namespace kairoq
