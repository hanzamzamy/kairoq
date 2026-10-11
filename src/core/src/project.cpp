// Project file reader and DC result summary.
#include "kairoq/archive.hpp"
#include "kairoq/dc.hpp"
#include "kairoq/elmer.hpp"
#include "kairoq/emi.hpp"
#include "kairoq/eye.hpp"
#include "kairoq/network.hpp"
#include "kairoq/openems.hpp"
#include "kairoq/pdn.hpp"
#include "kairoq/si.hpp"
#include "kairoq/timedomain.hpp"
#include "kairoq/thermal.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace fs = std::filesystem;

namespace {

double number(const YAML::Node& n, const char* key, double def) { return n[key] ? n[key].as<double>() : def; }

void check_keys(const YAML::Node& n, std::initializer_list<const char*> allowed, const std::string& where) {
    if (!n.IsMap()) throw std::runtime_error(where + " must be a mapping");
    for (const auto& kv : n) {
        const std::string k = kv.first.as<std::string>();
        if (std::find_if(allowed.begin(), allowed.end(), [&](const char* a) { return k == a; }) == allowed.end())
            throw std::runtime_error(where + ": unknown key '" + k + "'");
    }
}

// Area selection shared by ports and heat sources.
PortSpec parse_selection(const YAML::Node& n, const std::string& where, const std::string& default_name, bool component_needs_net) {
    PortSpec p;
    p.name = n["name"] ? n["name"].as<std::string>() : default_name;
    if (n["pad"]) p.pads.push_back(n["pad"].as<std::string>());
    if (n["pads"])
        for (const auto& v : n["pads"]) p.pads.push_back(v.as<std::string>());
    if (n["component"]) {
        p.component = n["component"].as<std::string>();
        if (component_needs_net && !n["net"]) throw std::runtime_error(where + ": component needs a net");
        if (n["net"]) p.net = n["net"].as<std::string>();
    }
    if (n["at_mm"]) {
        const auto v = n["at_mm"].as<std::vector<double>>();
        if (v.size() != 2) throw std::runtime_error(where + ": at_mm needs [x, y]");
        p.at = Point64(nm_from_mm(v[0]), nm_from_mm(v[1]));
        p.layer = n["layer"] ? n["layer"].as<int>() : 0;
        p.diameter_nm = number(n, "diameter_mm", 0.0) * kNmPerMm;
    }
    return p;
}

PortSpec parse_port(const YAML::Node& n, std::size_t i) {
    const std::string where = "port " + std::to_string(i + 1);
    check_keys(n, {"name", "pad", "pads", "component", "net", "at_mm", "layer", "diameter_mm", "voltage_v", "current_a"}, where);
    PortSpec p = parse_selection(n, where, "P" + std::to_string(i + 1), true);
    if (n["voltage_v"]) p.voltage_v = n["voltage_v"].as<double>();
    if (n["current_a"]) p.current_a = n["current_a"].as<double>();
    return p;
}

void parse_mesh(const YAML::Node& m, MeshSettings& mesh, const std::string& where) {
    check_keys(m, {"boundary_tolerance_um", "max_edge_mm", "min_size_um", "grading", "min_angle_deg", "min_edge_um"}, where);
    mesh.boundary_tolerance_nm = number(m, "boundary_tolerance_um", mesh.boundary_tolerance_nm / 1000.0) * 1000.0;
    mesh.max_edge_nm = number(m, "max_edge_mm", mesh.max_edge_nm / kNmPerMm) * kNmPerMm;
    mesh.min_angle_deg = number(m, "min_angle_deg", mesh.min_angle_deg);
    mesh.min_edge_nm = number(m, "min_edge_um", mesh.min_edge_nm / 1000.0) * 1000.0;
    mesh.min_size_nm = number(m, "min_size_um", mesh.min_size_nm / 1000.0) * 1000.0;
    mesh.grading = number(m, "grading", mesh.grading);
}

// Rectangles and polygons in millimetres: [{rect_mm: [x0, y0, x1, y1]}, {polygon_mm: [[x, y], ...]}].
PolygonSet parse_area(const YAML::Node& a, const std::string& where) {
    Paths64 paths;
    for (const auto& r : a) {
        check_keys(r, {"rect_mm", "polygon_mm"}, where + " area");
        Path64 path;
        if (r["rect_mm"]) {
            const auto v = r["rect_mm"].as<std::vector<double>>();
            if (v.size() != 4) throw std::runtime_error(where + ": rect_mm needs [x0, y0, x1, y1]");
            const double x0 = std::min(v[0], v[2]), x1 = std::max(v[0], v[2]);
            const double y0 = std::min(v[1], v[3]), y1 = std::max(v[1], v[3]);
            path = {Point64(nm_from_mm(x0), nm_from_mm(y0)), Point64(nm_from_mm(x1), nm_from_mm(y0)), Point64(nm_from_mm(x1), nm_from_mm(y1)),
                    Point64(nm_from_mm(x0), nm_from_mm(y1))};
        } else if (r["polygon_mm"]) {
            for (const auto& pt : r["polygon_mm"]) {
                const auto v = pt.as<std::vector<double>>();
                if (v.size() != 2) throw std::runtime_error(where + ": polygon_mm needs [x, y] points");
                path.push_back(Point64(nm_from_mm(v[0]), nm_from_mm(v[1])));
            }
            if (path.size() < 3) throw std::runtime_error(where + ": a polygon needs three points");
        } else {
            throw std::runtime_error(where + ": an area entry needs rect_mm or polygon_mm");
        }
        paths.push_back(path);
    }
    return union_polygons(paths);
}

SurfaceSpec parse_surface(const YAML::Node& n, const std::string& where) {
    check_keys(n, {"convection_w_m2k", "emissivity"}, where);
    SurfaceSpec sp;
    sp.convection_w_m2k = number(n, "convection_w_m2k", sp.convection_w_m2k);
    sp.emissivity = number(n, "emissivity", sp.emissivity);
    return sp;
}

}  // namespace

ThermalProject parse_thermal_project(std::string_view text, const std::string& base_dir) {
    const Project base = parse_project(text, base_dir);
    ThermalProject p;
    p.board = base.board;
    p.dc = base.dc;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node th = root["thermal"];
    if (!th) return p;
    check_keys(th, {"ambient_c", "top", "bottom", "copper_k", "dielectric", "layers", "mask_k", "outer_gaps", "heat_capacity",
                    "dielectric_slice_mm", "dielectric_first_slice_mm", "plating_thickness_um", "sources", "component_power",
                    "heatsinks", "transient", "coupling", "coupling_tolerance_k", "max_coupling_iterations", "coupling_acceleration", "target_error",
                    "max_passes", "max_triangles", "mesh", "solver", "elmer"},
               "thermal");
    ThermalSettings s;
    s.ambient_c = number(th, "ambient_c", s.ambient_c);
    if (th["top"]) s.top = parse_surface(th["top"], "thermal.top");
    if (th["bottom"]) s.bottom = parse_surface(th["bottom"], "thermal.bottom");
    s.copper_k = number(th, "copper_k", s.copper_k);
    if (const auto d = th["dielectric"]) {
        check_keys(d, {"k_xy", "k_z"}, "thermal.dielectric");
        s.dielectric_k_xy = number(d, "k_xy", s.dielectric_k_xy);
        s.dielectric_k_z = number(d, "k_z", s.dielectric_k_z);
    }
    if (const auto c = th["heat_capacity"]) {
        check_keys(c, {"copper", "dielectric", "mask"}, "thermal.heat_capacity");
        s.copper_heat_capacity = number(c, "copper", s.copper_heat_capacity);
        s.dielectric_heat_capacity = number(c, "dielectric", s.dielectric_heat_capacity);
        s.mask_heat_capacity = number(c, "mask", s.mask_heat_capacity);
    }
    if (const auto ls = th["layers"])
        for (const auto& kv : ls) {
            const std::string name = kv.first.as<std::string>();
            check_keys(kv.second, {"k_xy", "k_z", "heat_capacity"}, "thermal.layers." + name);
            if (kv.second["k_xy"] || kv.second["k_z"])
                s.layer_k[name] = {number(kv.second, "k_xy", s.dielectric_k_xy), number(kv.second, "k_z", s.dielectric_k_z)};
            if (kv.second["heat_capacity"]) s.layer_heat_capacity[name] = kv.second["heat_capacity"].as<double>();
        }
    s.mask_k = number(th, "mask_k", s.mask_k);
    if (const auto og = th["outer_gaps"]) {
        const std::string v = og.as<std::string>();
        if (v != "air" && v != "mask") throw std::runtime_error("thermal.outer_gaps must be air or mask");
        s.outer_gaps_mask = v == "mask";
    }
    s.dielectric_slice_nm = number(th, "dielectric_slice_mm", s.dielectric_slice_nm / 1e6) * 1e6;
    s.dielectric_first_slice_nm = number(th, "dielectric_first_slice_mm", s.dielectric_first_slice_nm / 1e6) * 1e6;
    if (!(s.dielectric_slice_nm > 0) || !(s.dielectric_first_slice_nm > 0))
        throw std::runtime_error("thermal.dielectric_slice_mm and dielectric_first_slice_mm must be positive");
    s.plating_thickness_nm = number(th, "plating_thickness_um", s.plating_thickness_nm / 1000.0) * 1000.0;
    if (const auto src = th["sources"])
        for (std::size_t i = 0; i < src.size(); ++i) {
            const auto n = src[i];
            const std::string where = "thermal source " + std::to_string(i + 1);
            check_keys(n, {"name", "pad", "pads", "component", "net", "at_mm", "layer", "diameter_mm", "power_w", "power_schedule",
                           "temperature_c", "junction_k_per_w", "heat_capacity_j_per_k"},
                       where);
            HeatSpec h;
            h.area = parse_selection(n, where, "S" + std::to_string(i + 1), false);
            if (n["power_w"]) h.power_w = n["power_w"].as<double>();
            if (const auto ps = n["power_schedule"])
                for (const auto& pt : ps) {
                    const auto v = pt.as<std::vector<double>>();
                    if (v.size() != 2) throw std::runtime_error(where + ": power_schedule needs [time_s, power_w] points");
                    h.power_schedule.push_back({v[0], v[1]});
                }
            h.heat_capacity_j_per_k = number(n, "heat_capacity_j_per_k", 0.0);
            if (n["temperature_c"]) h.temperature_c = n["temperature_c"].as<double>();
            if (n["junction_k_per_w"]) h.junction_k_per_w = n["junction_k_per_w"].as<double>();
            s.sources.push_back(std::move(h));
        }
    if (const auto cp = th["component_power"]) {
        auto one = [&](const YAML::Node& n) {
            check_keys(n, {"file", "reference_column", "power_column", "junction_column", "heat_capacity_column", "ignore_unknown"},
                       "thermal.component_power");
            if (!n["file"]) throw std::runtime_error("thermal.component_power needs a file");
            ComponentPowerTable t;
            const fs::path f(n["file"].as<std::string>());
            t.path = f.is_absolute() ? f.string() : (fs::path(base_dir) / f).lexically_normal().string();
            if (n["reference_column"]) t.reference_column = n["reference_column"].as<std::string>();
            if (n["power_column"]) t.power_column = n["power_column"].as<std::string>();
            if (n["junction_column"]) t.junction_column = n["junction_column"].as<std::string>();
            if (n["heat_capacity_column"]) t.heat_capacity_column = n["heat_capacity_column"].as<std::string>();
            if (n["ignore_unknown"]) t.ignore_unknown = n["ignore_unknown"].as<bool>();
            p.component_power.push_back(t);
        };
        if (cp.IsSequence())
            for (const auto& n : cp) one(n);
        else
            one(cp);
    }
    if (const auto hk = th["heatsinks"])
        for (std::size_t i = 0; i < hk.size(); ++i) {
            const auto n = hk[i];
            const std::string where = "heat sink " + std::to_string(i + 1);
            check_keys(n, {"name", "face", "area", "tim", "temperature_c", "resistance_k_per_w", "heat_capacity_j_per_k"}, where);
            HeatsinkSpec h;
            h.name = n["name"] ? n["name"].as<std::string>() : "HS" + std::to_string(i + 1);
            const std::string face = n["face"] ? n["face"].as<std::string>() : "bottom";
            if (face != "top" && face != "bottom") throw std::runtime_error(where + ": face is top or bottom");
            const bool bottom = face == "bottom";
            if (const auto a = n["area"]) {
                h.regions.push_back({bottom, parse_area(a, where)});
            } else {
                h.regions.push_back({bottom, {}});
            }
            if (const auto t = n["tim"]) {
                check_keys(t, {"thickness_mm", "k", "resistance_k_mm2_per_w"}, where + " tim");
                if (t["resistance_k_mm2_per_w"]) {
                    h.tim_m2k_per_w = t["resistance_k_mm2_per_w"].as<double>() * 1e-6;
                } else {
                    if (!t["thickness_mm"] || !t["k"]) throw std::runtime_error(where + ": tim needs thickness_mm and k, or resistance_k_mm2_per_w");
                    const double k = t["k"].as<double>();
                    if (!(k > 0)) throw std::runtime_error(where + ": tim k must be positive");
                    h.tim_m2k_per_w = t["thickness_mm"].as<double>() * 1e-3 / k;
                }
            }
            if (n["temperature_c"]) h.temperature_c = n["temperature_c"].as<double>();
            if (n["resistance_k_per_w"]) h.resistance_k_per_w = n["resistance_k_per_w"].as<double>();
            h.heat_capacity_j_per_k = number(n, "heat_capacity_j_per_k", 0.0);
            s.heatsinks.push_back(std::move(h));
        }
    if (const auto tr = th["transient"]) {
        check_keys(tr, {"duration_s", "step_s", "max_step_s", "steps_per_doubling", "start", "snapshots_s", "dc_update_k", "step_tolerance"},
                   "thermal.transient");
        TransientSettings t;
        if (!tr["duration_s"]) throw std::runtime_error("thermal.transient needs duration_s");
        t.duration_s = tr["duration_s"].as<double>();
        t.step_s = number(tr, "step_s", t.step_s);
        t.max_step_s = number(tr, "max_step_s", t.max_step_s);
        t.steps_per_doubling = static_cast<int>(number(tr, "steps_per_doubling", t.steps_per_doubling));
        const std::string start = tr["start"] ? tr["start"].as<std::string>() : "ambient";
        if (start != "ambient" && start != "steady") throw std::runtime_error("thermal.transient.start is ambient or steady");
        t.start_steady = start == "steady";
        if (tr["snapshots_s"]) t.snapshot_times_s = tr["snapshots_s"].as<std::vector<double>>();
        t.dc_update_k = number(tr, "dc_update_k", t.dc_update_k);
        t.step_tolerance = number(tr, "step_tolerance", t.step_tolerance);
        s.transient = t;
    }
    s.coupling_tolerance_k = number(th, "coupling_tolerance_k", s.coupling_tolerance_k);
    s.max_coupling_iterations = static_cast<int>(number(th, "max_coupling_iterations", s.max_coupling_iterations));
    if (th["coupling_acceleration"]) s.coupling_acceleration = th["coupling_acceleration"].as<bool>();
    s.target_error = number(th, "target_error", s.target_error);
    s.max_passes = static_cast<int>(number(th, "max_passes", s.max_passes));
    s.max_triangles = static_cast<std::size_t>(number(th, "max_triangles", static_cast<double>(s.max_triangles)));
    if (th["mesh"]) parse_mesh(th["mesh"], s.mesh, "thermal.mesh");
    if (const auto sv = th["solver"]) {
        check_keys(sv, {"tolerance", "max_iterations"}, "thermal.solver");
        s.tolerance = number(sv, "tolerance", s.tolerance);
        s.max_iterations = static_cast<int>(number(sv, "max_iterations", s.max_iterations));
    }
    // Electrothermal coupling when the project has DC ports, unless switched off.
    const bool coupling = th["coupling"] ? th["coupling"].as<bool>() : true;
    if (coupling && !p.dc.ports.empty()) s.dc = p.dc;
    p.thermal = s;
    return p;
}

ThermalProject load_thermal_project(const std::string& path) {
    return parse_thermal_project(read_file(path), fs::absolute(path).parent_path().string());
}

EmiProject parse_emi_project(std::string_view text, const std::string& base_dir) {
    const Project base = parse_project(text, base_dir);
    EmiProject p;
    p.board = base.board;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node em = root["emi"];
    if (!em) return p;
    check_keys(em, {"loops", "loop_settings", "return_path", "plates", "switch_nodes", "coupling_settings", "dielectric_epsilon_r", "mask_epsilon_r", "openems"},
               "emi");
    EmiSettings s;
    s.dielectric_epsilon_r = number(em, "dielectric_epsilon_r", s.dielectric_epsilon_r);
    s.mask_epsilon_r = number(em, "mask_epsilon_r", s.mask_epsilon_r);
    if (const auto ls = em["loop_settings"]) {
        check_keys(ls, {"window_margin_mm", "passive_min_area_mm2", "min_size_mm", "max_triangles", "depth_layers", "via_radius_mm",
                        "estimate_error", "distance_m", "max_frequency_hz", "moment_tolerance", "max_window_doublings", "loop_basis",
                        "frequencies_hz", "plating_um", "field", "board_height_m", "antenna_heights_m", "dense_matrix_mb"},
                   "emi.loop_settings");
        auto& l = s.loop;
        l.window_margin_nm = number(ls, "window_margin_mm", l.window_margin_nm / 1e6) * 1e6;
        l.passive_min_area_mm2 = number(ls, "passive_min_area_mm2", l.passive_min_area_mm2);
        l.min_size_nm = number(ls, "min_size_mm", l.min_size_nm / 1e6) * 1e6;
        l.max_triangles = static_cast<std::size_t>(number(ls, "max_triangles", static_cast<double>(l.max_triangles)));
        l.depth_layers = static_cast<int>(number(ls, "depth_layers", l.depth_layers));
        l.via_radius_nm = number(ls, "via_radius_mm", l.via_radius_nm / 1e6) * 1e6;
        if (ls["estimate_error"]) l.estimate_error = ls["estimate_error"].as<bool>();
        l.distance_m = number(ls, "distance_m", l.distance_m);
        l.moment_tolerance = number(ls, "moment_tolerance", l.moment_tolerance);
        l.max_window_doublings = static_cast<int>(number(ls, "max_window_doublings", l.max_window_doublings));
        if (ls["frequencies_hz"]) l.frequencies_hz = ls["frequencies_hz"].as<std::vector<double>>();
        l.plating_nm = number(ls, "plating_um", l.plating_nm * 1e-3) * 1e3;
        if (ls["field"]) {
            l.field_site = ls["field"].as<std::string>();
            if (l.field_site != "free" && l.field_site != "ground" && l.field_site != "small_loop")
                throw std::runtime_error("emi.loop_settings.field must be free, ground or small_loop");
        }
        l.board_height_m = number(ls, "board_height_m", l.board_height_m);
        l.dense_matrix_mb = number(ls, "dense_matrix_mb", l.dense_matrix_mb);
        if (const auto h = ls["antenna_heights_m"]) {
            const auto v = h.as<std::vector<double>>();
            if (v.size() != 2 || !(v[0] > 0) || v[1] < v[0]) throw std::runtime_error("emi.loop_settings.antenna_heights_m: [lowest, highest]");
            l.antenna_min_height_m = v[0];
            l.antenna_max_height_m = v[1];
        }
        if (ls["loop_basis"]) l.loop_basis = ls["loop_basis"].as<bool>();
        l.max_frequency_hz = number(ls, "max_frequency_hz", l.max_frequency_hz);
    }
    if (const auto oe = em["openems"]) {
        check_keys(oe, {"f_max_hz", "min_cell_mm", "board_cell_mm", "air_margin_mm", "port_resistance_ohm", "far_field_hz",
                        "end_criterion_db", "max_timesteps", "estimate"},
                   "emi.openems");
        auto& o = s.openems;
        o.f_max_hz = number(oe, "f_max_hz", o.f_max_hz);
        o.min_cell_nm = number(oe, "min_cell_mm", o.min_cell_nm / 1e6) * 1e6;
        o.board_cell_nm = number(oe, "board_cell_mm", o.board_cell_nm / 1e6) * 1e6;
        o.air_margin_nm = number(oe, "air_margin_mm", o.air_margin_nm / 1e6) * 1e6;
        o.port_resistance_ohm = number(oe, "port_resistance_ohm", o.port_resistance_ohm);
        if (oe["far_field_hz"]) o.far_field_hz = oe["far_field_hz"].as<std::vector<double>>();
        o.end_criterion_db = number(oe, "end_criterion_db", o.end_criterion_db);
        o.max_timesteps = static_cast<std::size_t>(number(oe, "max_timesteps", static_cast<double>(o.max_timesteps)));
        if (oe["estimate"]) o.estimate = oe["estimate"].as<bool>();
        if (!(o.port_resistance_ohm > 0) || !(o.end_criterion_db < 0))
            throw std::runtime_error("emi.openems: port_resistance_ohm must be positive, end_criterion_db negative");
    }
    if (const auto loops = em["loops"])
        for (std::size_t i = 0; i < loops.size(); ++i) {
            const auto n = loops[i];
            const std::string where = "loop " + std::to_string(i + 1);
            check_keys(n, {"name", "elements", "stitches", "waveform"}, where);
            LoopSpec lp;
            lp.name = n["name"] ? n["name"].as<std::string>() : "LOOP" + std::to_string(i + 1);
            if (!n["elements"] || n["elements"].size() == 0) throw std::runtime_error(where + " needs elements");
            for (const auto& e : n["elements"]) {
                check_keys(e, {"from", "to", "radius_mm", "esl_nh"}, where + " element");
                if (!e["from"] || !e["to"]) throw std::runtime_error(where + ": an element needs from and to (REF.PIN)");
                LoopElement el;
                el.from_pad = e["from"].as<std::string>();
                el.to_pad = e["to"].as<std::string>();
                el.radius_nm = number(e, "radius_mm", el.radius_nm / 1e6) * 1e6;
                el.extra_inductance_nh = number(e, "esl_nh", 0.0);
                if (!(el.radius_nm > 0)) throw std::runtime_error(where + ": radius_mm must be positive");
                lp.elements.push_back(el);
            }
            if (const auto st = n["stitches"])
                for (const auto& e : st) {
                    check_keys(e, {"from", "to", "radius_mm", "c_nf", "esr_mohm", "esl_nh"}, where + " stitch");
                    if (!e["from"] || !e["to"]) throw std::runtime_error(where + ": a stitch needs from and to (REF.PIN)");
                    LoopStitch sc;
                    sc.from_pad = e["from"].as<std::string>();
                    sc.to_pad = e["to"].as<std::string>();
                    sc.radius_nm = number(e, "radius_mm", sc.radius_nm / 1e6) * 1e6;
                    sc.capacitance_f = number(e, "c_nf", 0.0) * 1e-9;
                    sc.resistance_ohm = number(e, "esr_mohm", 0.0) * 1e-3;
                    sc.extra_inductance_nh = number(e, "esl_nh", 0.0);
                    if (!(sc.radius_nm > 0) || sc.capacitance_f < 0 || sc.resistance_ohm < 0 || sc.extra_inductance_nh < 0)
                        throw std::runtime_error(where + ": a stitch's radius must be positive and its c_nf, esr_mohm, esl_nh not negative");
                    lp.stitches.push_back(sc);
                }
            if (const auto w = n["waveform"]) {
                check_keys(w, {"amplitude_a", "frequency_hz", "duty", "rise_ns", "fall_ns"}, where + " waveform");
                lp.waveform.amplitude_a = number(w, "amplitude_a", 0.0);
                lp.waveform.frequency_hz = number(w, "frequency_hz", 0.0);
                lp.waveform.duty = number(w, "duty", lp.waveform.duty);
                lp.waveform.rise_s = number(w, "rise_ns", 0.0) * 1e-9;
                lp.waveform.fall_s = w["fall_ns"] ? w["fall_ns"].as<double>() * 1e-9 : -1.0;
                if (!(lp.waveform.amplitude_a > 0) || !(lp.waveform.frequency_hz > 0))
                    throw std::runtime_error(where + ": the waveform needs amplitude_a and frequency_hz");
            }
            s.loops.push_back(std::move(lp));
        }
    if (const auto rp = em["return_path"]) {
        ReturnPathSettings r;
        if (rp.IsMap()) {
            check_keys(rp, {"nets", "reference_nets", "plane_fraction", "min_area_mm2", "pad_clearance_mm"}, "emi.return_path");
            if (rp["nets"]) r.nets = rp["nets"].as<std::vector<std::string>>();
            if (rp["reference_nets"]) r.reference_nets = rp["reference_nets"].as<std::vector<std::string>>();
            r.plane_fraction = number(rp, "plane_fraction", r.plane_fraction);
            r.min_area_mm2 = number(rp, "min_area_mm2", r.min_area_mm2);
            r.pad_clearance_nm = number(rp, "pad_clearance_mm", r.pad_clearance_nm / 1e6) * 1e6;
            s.return_path = r;
        } else if (rp.as<bool>()) {
            s.return_path = r;
        }
    }
    if (const auto pl = em["plates"])
        for (std::size_t i = 0; i < pl.size(); ++i) {
            const auto n = pl[i];
            const std::string where = "plate " + std::to_string(i + 1);
            check_keys(n, {"name", "face", "area", "gap_mm", "epsilon_r"}, where);
            PlateSpec pt;
            pt.name = n["name"] ? n["name"].as<std::string>() : "PLATE" + std::to_string(i + 1);
            const std::string face = n["face"] ? n["face"].as<std::string>() : "bottom";
            if (face != "top" && face != "bottom") throw std::runtime_error(where + ": face is top or bottom");
            pt.bottom = face == "bottom";
            if (n["area"]) pt.area = parse_area(n["area"], where);
            pt.gap_nm = number(n, "gap_mm", 0.0) * 1e6;
            pt.gap_epsilon_r = number(n, "epsilon_r", 1.0);
            if (pt.gap_nm < 0 || !(pt.gap_epsilon_r > 0)) throw std::runtime_error(where + ": gap_mm must not be negative, epsilon_r positive");
            s.plates.push_back(std::move(pt));
        }
    if (const auto sw = em["switch_nodes"])
        for (std::size_t i = 0; i < sw.size(); ++i) {
            const auto n = sw[i];
            const std::string where = "switch node " + std::to_string(i + 1);
            check_keys(n, {"net", "plate", "dv_dt_v_per_ns"}, where);
            if (!n["net"] || !n["plate"]) throw std::runtime_error(where + " needs net and plate");
            s.switch_nodes.push_back({n["net"].as<std::string>(), n["plate"].as<std::string>(), number(n, "dv_dt_v_per_ns", 0.0)});
        }
    if (const auto cs = em["coupling_settings"]) {
        check_keys(cs, {"window_margin_mm", "max_panels", "min_size_mm"}, "emi.coupling_settings");
        s.coupling.min_size_nm = number(cs, "min_size_mm", 0.0) * 1e6;
        s.coupling.window_margin_nm = number(cs, "window_margin_mm", 0.0) * 1e6;
        s.coupling.max_panels = static_cast<std::size_t>(number(cs, "max_panels", double(s.coupling.max_panels)));
        if (s.coupling.max_panels < 10) throw std::runtime_error("emi.coupling_settings.max_panels: at least 10");
    }
    p.emi = s;
    return p;
}

namespace {

DielectricModel dielectric_from(const YAML::Node& n, const std::string& where, DielectricModel m) {
    check_keys(n, {"epsilon_r", "loss_tangent", "frequency_hz", "wideband"}, where);
    m.epsilon_r = number(n, "epsilon_r", m.epsilon_r);
    m.loss_tangent = number(n, "loss_tangent", m.loss_tangent);
    m.frequency_hz = number(n, "frequency_hz", m.frequency_hz);
    if (n["wideband"]) m.wideband = n["wideband"].as<bool>();
    if (!(m.epsilon_r >= 1) || m.loss_tangent < 0 || !(m.frequency_hz > 0))
        throw std::runtime_error(where + ": epsilon_r at least 1, loss_tangent not negative, frequency_hz positive");
    return m;
}

}  // namespace

SiProject parse_si_project(std::string_view text, const std::string& base_dir) {
    const Project base = parse_project(text, base_dir);
    SiProject p;
    p.board = base.board;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node si = root["si"];
    if (!si) return p;
    check_keys(si, {"nets", "reference_nets", "plane_fraction", "pairs", "pair_by_name", "targets", "reference_hz", "dispersion", "materials", "copper",
                    "etch_um", "etch_face", "cut", "min_finding_mm", "findings", "channels", "network", "crosstalk", "simulations", "eyes", "openems"},
               "si");
    SiSettings s;
    if (si["nets"]) s.nets = si["nets"].as<std::vector<std::string>>();
    if (si["reference_nets"]) s.reference_nets = si["reference_nets"].as<std::vector<std::string>>();
    s.plane_fraction = number(si, "plane_fraction", s.plane_fraction);
    if (si["pair_by_name"]) s.pair_by_name = si["pair_by_name"].as<bool>();
    if (const auto pairs = si["pairs"])
        for (std::size_t i = 0; i < pairs.size(); ++i) {
            const std::string where = "si pair " + std::to_string(i + 1);
            check_keys(pairs[i], {"p", "n"}, where);
            if (!pairs[i]["p"] || !pairs[i]["n"]) throw std::runtime_error(where + " needs p and n");
            s.pairs.push_back({pairs[i]["p"].as<std::string>(), pairs[i]["n"].as<std::string>()});
        }
    if (const auto targets = si["targets"])
        for (std::size_t i = 0; i < targets.size(); ++i) {
            const std::string where = "si target " + std::to_string(i + 1);
            const auto t = targets[i];
            check_keys(t, {"nets", "ohm", "differential", "tolerance"}, where);
            SiTarget k;
            if (t["nets"]) k.nets = t["nets"].as<std::vector<std::string>>();
            k.ohm = number(t, "ohm", 0.0);
            if (t["differential"]) k.differential = t["differential"].as<bool>();
            k.tolerance = number(t, "tolerance", k.tolerance);
            if (!(k.ohm > 0) || !(k.tolerance > 0)) throw std::runtime_error(where + ": ohm and tolerance must be positive");
            s.targets.push_back(k);
        }
    s.reference_hz = number(si, "reference_hz", s.reference_hz);
    if (si["dispersion"]) s.dispersion = si["dispersion"].as<bool>();
    if (const auto m = si["materials"]) {
        for (const auto& kv : m) {
            const std::string name = kv.first.as<std::string>();
            if (name == "default") s.default_dielectric = dielectric_from(kv.second, "si.materials.default", s.default_dielectric);
            else if (name == "mask") {
                check_keys(kv.second, {"epsilon_r", "loss_tangent"}, "si.materials.mask");
                s.mask_epsilon_r = number(kv.second, "epsilon_r", s.mask_epsilon_r);
                s.mask_loss_tangent = number(kv.second, "loss_tangent", s.mask_loss_tangent);
            } else {
                s.layer_materials[name] = dielectric_from(kv.second, "si.materials." + name, s.default_dielectric);
            }
        }
    }
    if (const auto c = si["copper"]) {
        check_keys(c, {"conductivity_s_per_m", "roughness"}, "si.copper");
        s.copper.conductivity_s_per_m = number(c, "conductivity_s_per_m", s.copper.conductivity_s_per_m);
        if (const auto r = c["roughness"]) {
            check_keys(r, {"model", "rms_um", "radius_um", "ratio"}, "si.copper.roughness");
            const std::string model = r["model"] ? r["model"].as<std::string>() : "none";
            if (model == "hammerstad") {
                s.copper.roughness = ConductorModel::Roughness::Hammerstad;
                s.copper.rms_m = number(r, "rms_um", 0.0) * 1e-6;
            } else if (model == "huray") {
                s.copper.roughness = ConductorModel::Roughness::Huray;
                s.copper.huray_radius_m = number(r, "radius_um", 0.0) * 1e-6;
                s.copper.huray_ratio = number(r, "ratio", 0.0);
            } else if (model != "none") {
                throw std::runtime_error("si.copper.roughness: model is none, hammerstad or huray");
            }
        }
    }
    s.etch_nm = number(si, "etch_um", 0.0) * 1000.0;
    if (const auto f = si["etch_face"]) {
        if (!f.IsMap()) throw std::runtime_error("si etch_face: a mapping of copper layer names to top or bottom");
        for (const auto& kv : f) {
            const std::string face = kv.second.as<std::string>();
            if (face != "top" && face != "bottom") throw std::runtime_error("si etch_face: " + kv.first.as<std::string>() + " is top or bottom");
            s.etch_bottom[kv.first.as<std::string>()] = face == "bottom";
        }
    }
    s.min_finding_nm = number(si, "min_finding_mm", s.min_finding_nm / 1e6) * 1e6;
    if (const auto f = si["findings"]) {
        check_keys(f, {"rise_ps", "fraction"}, "si.findings");
        s.finding_rise_s = number(f, "rise_ps", s.finding_rise_s * 1e12) * 1e-12;
        s.finding_rise_fraction = number(f, "fraction", s.finding_rise_fraction);
    }
    if (const auto x = si["crosstalk"]) {
        check_keys(x, {"rise_ps", "threshold"}, "si.crosstalk");
        s.rise_time_s = number(x, "rise_ps", s.rise_time_s * 1e12) * 1e-12;
        s.crosstalk_threshold = number(x, "threshold", s.crosstalk_threshold);
        if (!(s.rise_time_s > 0) || !(s.crosstalk_threshold > 0)) throw std::runtime_error("si.crosstalk: rise_ps and threshold must be positive");
    }
    if (const auto c = si["cut"]) {
        check_keys(c, {"step_mm", "lateral_mm", "neighbour_mm", "large_copper_mm", "target_error", "neighbour_heights", "reference_merge"}, "si.cut");
        s.neighbour_heights = number(c, "neighbour_heights", s.neighbour_heights);
        s.reference_merge = number(c, "reference_merge", s.reference_merge);
        s.step_nm = number(c, "step_mm", s.step_nm / 1e6) * 1e6;
        s.lateral_nm = number(c, "lateral_mm", s.lateral_nm / 1e6) * 1e6;
        s.neighbour_nm = number(c, "neighbour_mm", s.neighbour_nm / 1e6) * 1e6;
        s.large_copper_nm = number(c, "large_copper_mm", s.large_copper_nm / 1e6) * 1e6;
        s.target_error = number(c, "target_error", s.target_error);
        if (!(s.step_nm > 0) || !(s.lateral_nm > 0) || !(s.target_error > 0)) throw std::runtime_error("si.cut: values must be positive");
    }
    p.si = s;
    return p;
}

ChannelProject parse_channel_project(std::string_view text) {
    ChannelProject p;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node si = root["si"];
    if (!si) return p;
    if (const auto n = si["network"]) {
        check_keys(n, {"f_max_hz", "points", "vias", "via_model", "rise_ps", "tdr_ns"}, "si.network");
        p.network.f_max_hz = number(n, "f_max_hz", p.network.f_max_hz);
        p.network.points = static_cast<std::size_t>(number(n, "points", static_cast<double>(p.network.points)));
        if (n["vias"]) p.network.vias = n["vias"].as<bool>();
        if (n["via_model"]) {
            p.network.via_model = n["via_model"].as<std::string>();
            if (p.network.via_model != "physical" && p.network.via_model != "coaxial")
                throw std::runtime_error("si.network.via_model: physical or coaxial");
        }
        p.network.rise_time_s = number(n, "rise_ps", p.network.rise_time_s * 1e12) * 1e-12;
        p.network.tdr_duration_s = number(n, "tdr_ns", 0.0) * 1e-9;
        if (!(p.network.f_max_hz > 0) || p.network.points < 2 || !(p.network.rise_time_s > 0))
            throw std::runtime_error("si.network: f_max_hz and rise_ps must be positive, points at least 2");
    }
    if (const auto chans = si["channels"])
        for (std::size_t i = 0; i < chans.size(); ++i) {
            const auto c = chans[i];
            const std::string where = "si channel " + std::to_string(i + 1);
            check_keys(c, {"name", "ports", "differential", "impedance_ohm", "parts", "loads", "drive"}, where);
            ChannelSpec ch;
            ch.name = c["name"] ? c["name"].as<std::string>() : "CHANNEL" + std::to_string(i + 1);
            if (!c["ports"] || c["ports"].size() == 0) throw std::runtime_error(where + " needs ports (REF.PIN)");
            ch.ports = c["ports"].as<std::vector<std::string>>();
            if (const auto d = c["differential"])
                for (const auto& pr : d) {
                    const auto v = pr.as<std::vector<std::size_t>>();
                    if (v.size() != 2 || v[0] < 1 || v[1] < 1) throw std::runtime_error(where + ": differential pairs are [p, n] port numbers from 1");
                    ch.differential.emplace_back(v[0] - 1, v[1] - 1);
                }
            ch.impedance_ohm = number(c, "impedance_ohm", ch.impedance_ohm);
            if (const auto dr = c["drive"])
                for (const auto& v : dr) {
                    const auto k = v.as<std::size_t>();
                    if (k < 1) throw std::runtime_error(where + ": drive lists port numbers from 1");
                    ch.drive.push_back(k - 1);
                }
            if (const auto parts = c["parts"])
                for (const auto& pt : parts) {
                    check_keys(pt, {"from", "to", "r_ohm", "l_nh", "c_nf", "touchstone"}, where + " part");
                    if (!pt["from"] || !pt["to"]) throw std::runtime_error(where + ": a part needs from and to (REF.PIN)");
                    ChannelPart q;
                    q.from_pad = pt["from"].as<std::string>();
                    q.to_pad = pt["to"].as<std::string>();
                    q.r_ohm = number(pt, "r_ohm", 0.0);
                    q.l_h = number(pt, "l_nh", 0.0) * 1e-9;
                    q.c_f = number(pt, "c_nf", 0.0) * 1e-9;
                    if (pt["touchstone"]) q.touchstone = pt["touchstone"].as<std::string>();
                    ch.parts.push_back(q);
                }
            if (const auto loads = c["loads"])
                for (const auto& ld : loads) {
                    check_keys(ld, {"pad", "r_ohm", "c_pf"}, where + " load");
                    if (!ld["pad"]) throw std::runtime_error(where + ": a load needs a pad");
                    ch.loads.push_back({ld["pad"].as<std::string>(), number(ld, "r_ohm", 0.0), number(ld, "c_pf", 0.0) * 1e-12});
                }
            p.channels.push_back(std::move(ch));
        }
    return p;
}

std::vector<SimulationSpec> parse_simulations(std::string_view text) {
    std::vector<SimulationSpec> out;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node sims = root["si"] ? root["si"]["simulations"] : YAML::Node();
    if (!sims) return out;
    auto ibis_of = [](const YAML::Node& n, const std::string& where, std::string& file, std::string& model, int& corner, bool& package) {
        check_keys(n, {"file", "model", "corner", "package"}, where + " ibis");
        if (!n["file"] || !n["model"]) throw std::runtime_error(where + ": ibis needs file and model");
        file = n["file"].as<std::string>();
        model = n["model"].as<std::string>();
        const std::string c = n["corner"] ? n["corner"].as<std::string>() : "typ";
        if (c != "typ" && c != "min" && c != "max") throw std::runtime_error(where + ": corner is typ, min or max");
        corner = c == "typ" ? 0 : c == "min" ? 1 : 2;
        if (n["package"]) package = n["package"].as<bool>();
    };
    for (std::size_t i = 0; i < sims.size(); ++i) {
        const auto n = sims[i];
        const std::string where = "si simulation " + std::to_string(i + 1);
        check_keys(n, {"name", "channel", "drivers", "receivers", "duration_ns", "oversample"}, where);
        SimulationSpec sp;
        sp.name = n["name"] ? n["name"].as<std::string>() : "SIM" + std::to_string(i + 1);
        if (!n["channel"]) throw std::runtime_error(where + " needs a channel");
        sp.channel = n["channel"].as<std::string>();
        sp.duration_s = number(n, "duration_ns", 0.0) * 1e-9;
        sp.oversample = static_cast<std::size_t>(number(n, "oversample", double(sp.oversample)));
        for (const auto& d : n["drivers"]) {
            check_keys(d, {"port", "pattern", "bit_rate", "bits", "delay_ps", "v_low", "v_high", "r_ohm", "rise_ps", "fall_ps", "ibis"},
                       where + " driver");
            DriverSpec ds;
            if (!d["port"]) throw std::runtime_error(where + ": a driver needs a port (from 1)");
            ds.port = d["port"].as<std::size_t>() - 1;
            if (d["pattern"]) ds.stimulus.pattern = d["pattern"].as<std::string>();
            ds.stimulus.bit_rate = number(d, "bit_rate", ds.stimulus.bit_rate);
            ds.stimulus.bits = static_cast<std::size_t>(number(d, "bits", double(ds.stimulus.bits)));
            ds.stimulus.delay_s = number(d, "delay_ps", ds.stimulus.delay_s * 1e12) * 1e-12;
            ds.v_low = number(d, "v_low", ds.v_low);
            ds.v_high = number(d, "v_high", ds.v_high);
            ds.r_ohm = number(d, "r_ohm", ds.r_ohm);
            ds.rise_s = number(d, "rise_ps", ds.rise_s * 1e12) * 1e-12;
            ds.fall_s = d["fall_ps"] ? d["fall_ps"].as<double>() * 1e-12 : -1;
            if (d["ibis"]) ibis_of(d["ibis"], where + " driver", ds.ibis_file, ds.ibis_model, ds.corner, ds.package);
            sp.drivers.push_back(ds);
        }
        for (const auto& r : n["receivers"]) {
            check_keys(r, {"port", "r_ohm", "v_term", "c_pf", "ibis"}, where + " receiver");
            ReceiverSpec rs;
            if (!r["port"]) throw std::runtime_error(where + ": a receiver needs a port (from 1)");
            rs.port = r["port"].as<std::size_t>() - 1;
            rs.r_ohm = number(r, "r_ohm", 0.0);
            rs.v_term = number(r, "v_term", 0.0);
            rs.c_f = number(r, "c_pf", 0.0) * 1e-12;
            if (r["ibis"]) ibis_of(r["ibis"], where + " receiver", rs.ibis_file, rs.ibis_model, rs.corner, rs.package);
            sp.receivers.push_back(rs);
        }
        out.push_back(std::move(sp));
    }
    return out;
}

PdnSettings parse_pdn(std::string_view text) {
    PdnSettings s;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node n = root["pdn"];
    if (!n) throw std::runtime_error("the project has no pdn section");
    check_keys(n, {"power_net", "ground_net", "f_min_hz", "f_max_hz", "points", "target_ohm", "loads", "capacitors", "other_capacitors",
                   "capacitor_prefix", "regulator", "plane_fraction", "via_reach_mm", "dielectric"},
               "pdn");
    if (!n["power_net"] || !n["ground_net"]) throw std::runtime_error("pdn needs power_net and ground_net");
    s.power_net = n["power_net"].as<std::string>();
    s.ground_net = n["ground_net"].as<std::string>();
    s.f_min_hz = number(n, "f_min_hz", s.f_min_hz);
    s.f_max_hz = number(n, "f_max_hz", s.f_max_hz);
    s.points = static_cast<std::size_t>(number(n, "points", double(s.points)));
    if (n["target_ohm"]) s.target_ohm = n["target_ohm"].as<double>();
    if (n["loads"]) s.loads = n["loads"].as<std::vector<std::string>>();
    auto capacitor = [](const YAML::Node& c, const std::string& where) {
        check_keys(c, {"component", "c_uf", "esr_mohm", "esl_nh", "mount_nh"}, where);
        PdnCapacitor k;
        if (c["component"]) k.component = c["component"].as<std::string>();
        if (!c["c_uf"]) throw std::runtime_error(where + " needs c_uf");
        k.c_f = c["c_uf"].as<double>() * 1e-6;
        k.esr_ohm = number(c, "esr_mohm", 0.0) * 1e-3;
        k.esl_h = number(c, "esl_nh", 0.0) * 1e-9;
        if (c["mount_nh"]) k.mount_h = c["mount_nh"].as<double>() * 1e-9;
        return k;
    };
    if (n["capacitors"])
        for (std::size_t i = 0; i < n["capacitors"].size(); ++i) {
            auto c = capacitor(n["capacitors"][i], "pdn capacitor " + std::to_string(i + 1));
            if (c.component.empty()) throw std::runtime_error("pdn capacitor " + std::to_string(i + 1) + " needs component");
            s.capacitors.push_back(c);
        }
    if (n["other_capacitors"]) s.other_capacitors = capacitor(n["other_capacitors"], "pdn other_capacitors");
    if (n["capacitor_prefix"]) s.capacitor_prefix = n["capacitor_prefix"].as<std::string>();
    if (const auto r = n["regulator"]) {
        check_keys(r, {"component", "r_mohm", "l_nh"}, "pdn regulator");
        PdnSettings::Regulator g;
        if (!r["component"]) throw std::runtime_error("pdn regulator needs component");
        g.component = r["component"].as<std::string>();
        g.r_ohm = number(r, "r_mohm", g.r_ohm * 1e3) * 1e-3;
        g.l_h = number(r, "l_nh", g.l_h * 1e9) * 1e-9;
        s.regulator = g;
    }
    s.plane_fraction = number(n, "plane_fraction", s.plane_fraction);
    s.via_reach_nm = number(n, "via_reach_mm", s.via_reach_nm * 1e-6) * 1e6;
    if (const auto d = n["dielectric"]) {
        check_keys(d, {"epsilon_r", "loss_tangent", "frequency_hz", "wideband"}, "pdn dielectric");
        s.dielectric.epsilon_r = number(d, "epsilon_r", s.dielectric.epsilon_r);
        s.dielectric.loss_tangent = number(d, "loss_tangent", s.dielectric.loss_tangent);
        s.dielectric.frequency_hz = number(d, "frequency_hz", s.dielectric.frequency_hz);
        if (d["wideband"]) s.dielectric.wideband = d["wideband"].as<bool>();
    }
    return s;
}

std::vector<EyeSpec> parse_eyes(std::string_view text) {
    std::vector<EyeSpec> out;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node eyes = root["si"] ? root["si"]["eyes"] : YAML::Node();
    if (!eyes) return out;
    for (std::size_t i = 0; i < eyes.size(); ++i) {
        const auto n = eyes[i];
        const std::string where = "si eye " + std::to_string(i + 1);
        check_keys(n, {"name", "channel", "tx", "rx", "bit_rate", "swing_v", "r_source_ohm", "rise_ps", "r_term_ohm", "c_term_pf", "pattern",
                       "bits", "samples_per_ui", "noise_mv", "rj_ps", "dj_ps", "ber", "ffe", "ffe_main", "ctle", "dfe"},
                   where);
        EyeSpec e;
        e.name = n["name"] ? n["name"].as<std::string>() : "EYE" + std::to_string(i + 1);
        if (!n["channel"] || !n["tx"] || !n["rx"]) throw std::runtime_error(where + " needs channel, tx and rx (ports from 1)");
        e.channel = n["channel"].as<std::string>();
        const auto tx = n["tx"].as<std::vector<std::size_t>>(), rx = n["rx"].as<std::vector<std::size_t>>();
        if (tx.empty() || tx.size() > 2 || rx.empty() || rx.size() > 2 || std::count(tx.begin(), tx.end(), 0u) || std::count(rx.begin(), rx.end(), 0u))
            throw std::runtime_error(where + ": tx and rx are [p] or [p, n], ports from 1");
        e.tx_p = tx[0] - 1;
        if (tx.size() == 2) e.tx_n = tx[1] - 1;
        e.rx_p = rx[0] - 1;
        if (rx.size() == 2) e.rx_n = rx[1] - 1;
        e.bit_rate = number(n, "bit_rate", e.bit_rate);
        e.swing_v = number(n, "swing_v", e.swing_v);
        e.r_source_ohm = number(n, "r_source_ohm", e.r_source_ohm);
        e.rise_s = number(n, "rise_ps", e.rise_s * 1e12) * 1e-12;
        e.r_term_ohm = number(n, "r_term_ohm", e.r_term_ohm);
        e.c_term_f = number(n, "c_term_pf", 0.0) * 1e-12;
        if (n["pattern"]) e.pattern = n["pattern"].as<std::string>();
        e.bits = static_cast<std::size_t>(number(n, "bits", double(e.bits)));
        e.samples_per_ui = static_cast<std::size_t>(number(n, "samples_per_ui", double(e.samples_per_ui)));
        e.noise_rms_v = number(n, "noise_mv", 0.0) * 1e-3;
        e.rj_rms_s = number(n, "rj_ps", 0.0) * 1e-12;
        e.dj_pp_s = number(n, "dj_ps", 0.0) * 1e-12;
        e.ber = number(n, "ber", e.ber);
        if (n["ffe"]) e.ffe = n["ffe"].as<std::vector<double>>();
        const double main = number(n, "ffe_main", 1.0);
        if (main < 1) throw std::runtime_error(where + ": ffe_main counts from 1");
        e.ffe_main = static_cast<std::size_t>(main) - 1;
        if (const auto c = n["ctle"]) {
            check_keys(c, {"dc_gain_db", "f_zero_hz", "f_pole1_hz", "f_pole2_hz"}, where + " ctle");
            e.ctle.dc_gain_db = number(c, "dc_gain_db", 0.0);
            e.ctle.f_zero_hz = number(c, "f_zero_hz", 0.0);
            e.ctle.f_pole1_hz = number(c, "f_pole1_hz", 0.0);
            e.ctle.f_pole2_hz = number(c, "f_pole2_hz", 0.0);
            if (!(e.ctle.f_zero_hz > 0)) throw std::runtime_error(where + ": the ctle needs f_zero_hz");
        }
        if (const auto d = n["dfe"]) {
            check_keys(d, {"taps", "limit"}, where + " dfe");
            e.dfe_taps = static_cast<std::size_t>(number(d, "taps", 0.0));
            e.dfe_limit = number(d, "limit", 1.0);
        }
        if (!(e.bit_rate > 0) || !(e.r_source_ohm > 0) || !(e.r_term_ohm > 0) || e.samples_per_ui < 4 || !(e.ber > 0 && e.ber < 0.5) ||
            e.noise_rms_v < 0 || e.rj_rms_s < 0 || e.dj_pp_s < 0)
            throw std::runtime_error(where + ": bit_rate, r_source_ohm and r_term_ohm positive, samples_per_ui at least 4, ber in (0, 0.5), "
                                             "noise and jitter not negative");
        if (e.pattern == "edge" || e.pattern == "pulse") throw std::runtime_error(where + ": the pattern is a PRBS or a string of 0 and 1");
        if (!e.ffe.empty() && e.ffe_main >= e.ffe.size()) throw std::runtime_error(where + ": ffe_main beyond the taps");
        out.push_back(std::move(e));
    }
    return out;
}

OpenEmsChannelSettings parse_openems_channel(std::string_view text) {
    OpenEmsChannelSettings s;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node n = root["si"] ? root["si"]["openems"] : YAML::Node();
    if (!n) return s;
    check_keys(n, {"f_max_hz", "points", "margin_mm", "min_cell_mm", "board_cell_mm", "air_margin_mm", "far_field_hz", "end_criterion_db",
                   "max_timesteps"},
               "si openems");
    s.f_max_hz = number(n, "f_max_hz", s.f_max_hz);
    s.points = static_cast<std::size_t>(number(n, "points", double(s.points)));
    s.margin_nm = number(n, "margin_mm", s.margin_nm * 1e-6) * 1e6;
    s.min_cell_nm = number(n, "min_cell_mm", 0.0) * 1e6;
    s.board_cell_nm = number(n, "board_cell_mm", 0.0) * 1e6;
    s.air_margin_nm = number(n, "air_margin_mm", 0.0) * 1e6;
    if (n["far_field_hz"]) s.far_field_hz = n["far_field_hz"].as<std::vector<double>>();
    s.end_criterion_db = number(n, "end_criterion_db", s.end_criterion_db);
    s.max_timesteps = static_cast<std::size_t>(number(n, "max_timesteps", double(s.max_timesteps)));
    return s;
}

ElmerSettings parse_elmer_settings(std::string_view text, const std::string& section) {
    ElmerSettings e;
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node n = root[section] ? root[section]["elmer"] : YAML::Node();
    if (!n) return e;
    check_keys(n, {"max_edge_mm", "min_size_mm", "grading", "copper_divisions", "solver"}, section + ".elmer");
    e.max_edge_nm = number(n, "max_edge_mm", e.max_edge_nm * 1e-6) * 1e6;
    e.min_size_nm = number(n, "min_size_mm", e.min_size_nm * 1e-6) * 1e6;
    e.grading = number(n, "grading", e.grading);
    e.copper_divisions = static_cast<int>(number(n, "copper_divisions", e.copper_divisions));
    if (n["solver"]) e.solver = n["solver"].as<std::string>();
    if (!(e.max_edge_nm > 0) || !(e.min_size_nm > 0) || !(e.grading > 0)) throw std::runtime_error(section + ".elmer: sizes and grading positive");
    return e;
}

SiProject load_si_project(const std::string& path) {
    return parse_si_project(read_file(path), fs::absolute(path).parent_path().string());
}

EmiProject load_emi_project(const std::string& path) {
    return parse_emi_project(read_file(path), fs::absolute(path).parent_path().string());
}

namespace {

// One CSV record per line; fields may be quoted ("R1, R2"), with "" for a quote.
std::vector<std::vector<std::string>> read_csv(const std::string& text) {
    // The delimiter is the one of comma, semicolon and tab most frequent in the first line.
    const std::string first = text.substr(0, text.find('\n'));
    char delim = ',';
    std::size_t best = 0;
    for (char c : {',', ';', '\t'}) {
        const auto k = static_cast<std::size_t>(std::count(first.begin(), first.end(), c));
        if (k > best) {
            best = k;
            delim = c;
        }
    }
    std::vector<std::vector<std::string>> rows;
    std::vector<std::string> row;
    std::string field;
    bool quoted = false, any = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (quoted) {
            if (c == '"' && i + 1 < text.size() && text[i + 1] == '"') {
                field += '"';
                ++i;
            } else if (c == '"') {
                quoted = false;
            } else {
                field += c;
            }
            continue;
        }
        if (c == '"') {
            quoted = true;
            any = true;
        } else if (c == delim) {
            row.push_back(field);
            field.clear();
            any = true;
        } else if (c == '\n' || c == '\r') {
            if (c == '\r' && i + 1 < text.size() && text[i + 1] == '\n') ++i;
            if (any || !field.empty()) {
                row.push_back(field);
                rows.push_back(row);
            }
            row.clear();
            field.clear();
            any = false;
        } else {
            field += c;
            any = true;
        }
    }
    if (any || !field.empty()) {
        row.push_back(field);
        rows.push_back(row);
    }
    return rows;
}

std::string trim(const std::string& v) {
    const auto a = v.find_first_not_of(" \t");
    if (a == std::string::npos) return "";
    return v.substr(a, v.find_last_not_of(" \t") - a + 1);
}

std::string lower(std::string v) {
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return v;
}

// A number with an optional unit suffix of the given base unit and an SI prefix (m, k):
// "0.25", "250mW", "0.25 W".
std::optional<double> parse_quantity(const std::string& text, const std::string& unit) {
    std::string v = trim(text);
    if (v.empty()) return std::nullopt;
    std::size_t used = 0;
    double x = 0;
    try {
        x = std::stod(v, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("'" + text + "' is not a number");
    }
    std::string rest = trim(v.substr(used));
    double scale = 1.0;
    const std::string u = lower(unit);
    std::string r = lower(rest);
    if (!u.empty() && r.size() >= u.size() && r.compare(r.size() - u.size(), u.size(), u) == 0) r = trim(r.substr(0, r.size() - u.size()));
    if (r == "m") scale = 1e-3;
    else if (r == "u" || r == "\xc2\xb5") scale = 1e-6;
    else if (r == "k") scale = 1e3;
    else if (!r.empty()) throw std::runtime_error("'" + text + "': unknown unit");
    return x * scale;
}

}  // namespace

std::vector<HeatSpec> read_component_power(const Board& board, const ComponentPowerTable& t, std::vector<std::string>& warnings) {
    const auto rows = read_csv(read_file(t.path));
    const std::string where = fs::path(t.path).filename().string();
    if (rows.empty()) throw std::runtime_error(where + " is empty");
    const auto& head = rows[0];
    auto column = [&](const std::string& name, std::initializer_list<const char*> defaults, bool required) -> long {
        auto find = [&](const std::string& n) -> long {
            for (std::size_t i = 0; i < head.size(); ++i)
                if (lower(trim(head[i])) == lower(n)) return static_cast<long>(i);
            return -1;
        };
        if (!name.empty()) {
            const long i = find(name);
            if (i < 0) throw std::runtime_error(where + ": no column '" + name + "'");
            return i;
        }
        for (const char* d : defaults)
            if (const long i = find(d); i >= 0) return i;
        if (required) throw std::runtime_error(where + ": no column for " + std::string(*defaults.begin()) + "; name it in the project file");
        return -1;
    };
    const long cref = column(t.reference_column, {"Reference", "Designator", "Ref", "RefDes", "References"}, true);
    const long cpow = column(t.power_column, {"Power", "Power_W", "Power (W)", "P", "P_W"}, true);
    const long cjun = t.junction_column.empty() ? -1 : column(t.junction_column, {}, true);
    const long ccap = t.heat_capacity_column.empty() ? -1 : column(t.heat_capacity_column, {}, true);
    std::set<std::string> on_board;
    for (const auto& pad : board.pads)
        if (!pad.component.empty()) on_board.insert(pad.component);
    std::vector<HeatSpec> out;
    std::set<std::string> seen;
    std::vector<std::string> unknown;
    for (std::size_t r = 1; r < rows.size(); ++r) {
        const auto& row = rows[r];
        auto cell = [&](long c) { return c >= 0 && static_cast<std::size_t>(c) < row.size() ? row[static_cast<std::size_t>(c)] : std::string(); };
        const std::string line = where + " line " + std::to_string(r + 1);
        std::optional<double> power;
        try {
            power = parse_quantity(cell(cpow), "W");
        } catch (const std::exception& e) {
            throw std::runtime_error(line + ": " + e.what());
        }
        if (!power) continue;
        std::optional<double> rjc, cap;
        try {
            if (cjun >= 0) rjc = parse_quantity(cell(cjun), "K/W");
            if (ccap >= 0) cap = parse_quantity(cell(ccap), "J/K");
        } catch (const std::exception& e) {
            throw std::runtime_error(line + ": " + e.what());
        }
        // References separated by commas, semicolons or spaces.
        std::string refs = cell(cref);
        for (char& c : refs)
            if (c == ',' || c == ';' || c == '\t') c = ' ';
        std::istringstream in(refs);
        std::string ref;
        bool any = false;
        while (in >> ref) {
            any = true;
            if (!on_board.count(ref)) {
                unknown.push_back(ref);
                continue;
            }
            if (!seen.insert(ref).second) throw std::runtime_error(line + ": " + ref + " is listed twice");
            HeatSpec h;
            h.area.name = ref;
            h.area.component = ref;
            h.power_w = *power;
            h.junction_k_per_w = rjc;
            h.heat_capacity_j_per_k = cap.value_or(0.0);
            out.push_back(std::move(h));
        }
        if (!any) throw std::runtime_error(line + ": a power without a reference");
    }
    if (!unknown.empty()) {
        std::string list;
        for (std::size_t i = 0; i < unknown.size() && i < 10; ++i) list += (i ? ", " : "") + unknown[i];
        if (unknown.size() > 10) list += ", ...";
        const std::string msg = where + ": " + std::to_string(unknown.size()) + " components not on the board (" + list + ")";
        if (!t.ignore_unknown) throw std::runtime_error(msg);
        warnings.push_back(msg);
    }
    return out;
}

ThermalSettings resolve_thermal_project(const Board& board, const ThermalProject& p, std::vector<std::string>& warnings) {
    if (!p.thermal) throw std::runtime_error("the project has no thermal section");
    ThermalSettings s = *p.thermal;
    std::set<std::string> named;
    for (const auto& h : s.sources) named.insert(h.area.name);
    for (const auto& t : p.component_power)
        for (auto& h : read_component_power(board, t, warnings)) {
            // A source given in the project file takes precedence over the table.
            if (named.count(h.area.name)) {
                warnings.push_back(h.area.name + ": the source in the project file replaces the table's");
                continue;
            }
            named.insert(h.area.name);
            s.sources.push_back(std::move(h));
        }
    return s;
}

std::string thermal_summary(const Board& b, const ThermalResult& r) {
    std::ostringstream s;
    char line[256];
    std::snprintf(line, sizeof line, "mesh         %zu triangles, %zu unknowns, estimated error %.3f after %zu passes\n", r.triangles,
                  r.unknowns, r.estimated_error, r.passes.size());
    s << line;
    if (r.dc) {
        std::snprintf(line, sizeof line, "coupling     %d iterations, last change %.3g K\n", r.coupling_iterations, r.coupling_change_k);
        s << line;
    }
    std::snprintf(line, sizeof line, "time         mesh %.2f s, assembly %.2f s, solve %.2f s, DC %.2f s\n", r.seconds_mesh,
                  r.seconds_assemble, r.seconds_solve, r.seconds_dc);
    s << line;
    std::snprintf(line, sizeof line, "maximum      %.2f C on L%d at (%.3f, %.3f) mm\n", r.max_c, r.max_layer, to_mm(r.max_x_nm),
                  to_mm(r.max_y_nm));
    s << line;
    std::snprintf(line, sizeof line, "heat in      %.6g W from sources, %.6g W Joule heat\n", r.heat_sources_w, r.heat_joule_w);
    s << line;
    std::snprintf(line, sizeof line, "heat out     %.6g W top, %.6g W bottom, %.6g W to heat sinks, %.6g W to fixed temperatures\n",
                  r.heat_top_w, r.heat_bottom_w, r.heat_sinks_w, r.heat_fixed_w);
    s << line;
    if (!r.heatsinks.empty()) s << "heat sinks\n";
    for (const auto& h : r.heatsinks) {
        std::snprintf(line, sizeof line, "  %-16s %10.4f W  at %8.2f C%s\n", h.name.c_str(), h.heat_w, h.temperature_c,
                      h.temperature_set ? "  (fixed)" : "");
        s << line;
    }
    if (!r.sources.empty()) s << "sources\n";
    for (const auto& src : r.sources) {
        std::snprintf(line, sizeof line, "  %-16s %10.4f W  mean %8.2f C  max %8.2f C", src.name.c_str(), src.power_w, src.mean_c,
                      src.max_c);
        s << line;
        if (src.junction_c) {
            std::snprintf(line, sizeof line, "  junction %8.2f C", *src.junction_c);
            s << line;
        }
        s << (src.temperature_set ? "  (fixed)\n" : "\n");
    }
    if (r.transient) {
        const auto& t = *r.transient;
        std::snprintf(line, sizeof line, "transient    %.6g s in %d steps, %d matrix set-ups, %d DC solutions; values at the end above\n",
                      t.time_s.empty() ? 0.0 : t.time_s.back(), t.steps, t.matrix_setups, t.dc_solutions);
        s << line;
        std::snprintf(line, sizeof line, "  time       steady solutions %.2f s, multigrid set-ups %.2f s, %ld solver iterations in all\n",
                      t.seconds_steady, t.seconds_setup, t.solver_iterations);
        s << line;
        const auto peak = [](const std::vector<double>& v, double& at, const std::vector<double>& time) {
            double m = -INFINITY;
            for (std::size_t i = 0; i < v.size(); ++i)
                if (v[i] > m) {
                    m = v[i];
                    at = time[i];
                }
            return m;
        };
        double at = 0;
        const double mx = peak(t.max_c, at, t.time_s);
        std::snprintf(line, sizeof line, "  highest    %.2f C at %.4g s\n", mx, at);
        s << line;
        for (std::size_t a = 0; a < r.sources.size() && a < t.source_max_c.size(); ++a) {
            const bool junction = !t.source_junction_c[a].empty() && !std::isnan(t.source_junction_c[a].back());
            const double m = peak(junction ? t.source_junction_c[a] : t.source_max_c[a], at, t.time_s);
            std::snprintf(line, sizeof line, "  %-16s highest %s %8.2f C at %.4g s\n", r.sources[a].name.c_str(),
                          junction ? "junction" : "temperature", m, at);
            s << line;
        }
    }
    if (r.dc) s << dc_summary(b, *r.dc);
    for (const auto& w : r.warnings) s << "warning: " << w << "\n";
    return s.str();
}


Project parse_project(std::string_view text, const std::string& base_dir) {
    const YAML::Node root = YAML::Load(std::string(text));
    check_keys(root, {"board", "dc", "thermal", "emi", "si", "pdn"}, "project");
    Project p;
    if (root["board"]) {
        auto add = [&](const std::string& path) {
            const fs::path q(path);
            p.board.push_back(q.is_absolute() ? path : (fs::path(base_dir) / q).lexically_normal().string());
        };
        if (root["board"].IsSequence())
            for (const auto& v : root["board"]) add(v.as<std::string>());
        else
            add(root["board"].as<std::string>());
    }
    if (const auto dc = root["dc"]) {
        check_keys(dc, {"temperature_c", "copper_conductivity_s_per_m", "plating_thickness_um", "plating_conductivity_s_per_m",
                        "target_error", "max_passes", "max_triangles", "mesh", "solver", "ports", "elmer"},
                   "dc");
        DcSettings& s = p.dc;
        s.temperature_c = number(dc, "temperature_c", s.temperature_c);
        if (dc["copper_conductivity_s_per_m"]) s.copper_conductivity_s_per_m = dc["copper_conductivity_s_per_m"].as<double>();
        s.plating_thickness_nm = number(dc, "plating_thickness_um", s.plating_thickness_nm / 1000.0) * 1000.0;
        s.plating_conductivity_s_per_m = number(dc, "plating_conductivity_s_per_m", s.plating_conductivity_s_per_m);
        s.target_error = number(dc, "target_error", s.target_error);
        s.max_passes = static_cast<int>(number(dc, "max_passes", s.max_passes));
        s.max_triangles = static_cast<std::size_t>(number(dc, "max_triangles", static_cast<double>(s.max_triangles)));
        if (const auto m = dc["mesh"]) parse_mesh(m, s.mesh, "dc.mesh");
        if (const auto sv = dc["solver"]) {
            check_keys(sv, {"tolerance", "max_iterations"}, "dc.solver");
            s.tolerance = number(sv, "tolerance", s.tolerance);
            s.max_iterations = static_cast<int>(number(sv, "max_iterations", s.max_iterations));
        }
        if (const auto ports = dc["ports"])
            for (std::size_t i = 0; i < ports.size(); ++i) s.ports.push_back(parse_port(ports[i], i));
    }
    return p;
}

Project load_project(const std::string& path) {
    return parse_project(read_file(path), fs::absolute(path).parent_path().string());
}

std::string dc_summary(const Board& b, const DcResult& r) {
    std::ostringstream s;
    char line[256];
    std::snprintf(line, sizeof line, "mesh         %zu islands, %zu triangles, %zu unknowns\n", r.islands, r.triangles, r.unknowns);
    s << line;
    std::snprintf(line, sizeof line, "solver       %d iterations, residual %.2e\n", r.iterations, r.residual);
    s << line;
    std::snprintf(line, sizeof line, "accuracy     estimated error %.3f (energy norm) after %zu passes\n", r.estimated_error,
                  r.passes.size());
    s << line;
    std::snprintf(line, sizeof line, "time         mesh %.2f s, assembly %.2f s, solve %.2f s\n", r.seconds_mesh, r.seconds_assemble,
                  r.seconds_solve);
    s << line;
    s << "ports\n";
    for (const auto& p : r.ports) {
        std::snprintf(line, sizeof line, "  %-16s %12.6f V %s %12.6f A%s\n", p.name.c_str(), p.voltage_v, p.voltage_set ? "(set)" : "     ",
                      p.current_a, p.voltage_set ? "" : " (set)");
        s << line;
    }
    std::snprintf(line, sizeof line, "power        %.6g W in copper, %.6g W in hole plating; %.6g W delivered by the ports\n",
                  r.power_copper_w, r.power_barrels_w, r.power_ports_w);
    s << line;
    if (!r.barrels.empty()) {
        auto top = r.barrels;
        std::sort(top.begin(), top.end(), [](const auto& a, const auto& c) { return std::abs(a.current_a) > std::abs(c.current_a); });
        s << "hole currents (largest)\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(5, top.size()); ++i) {
            const auto& h = top[i];
            const auto& d = b.drills[h.drill];
            std::snprintf(line, sizeof line, "  (%9.3f, %9.3f) mm  L%d-L%d  %-14s %10.4f A  %.3g mOhm\n", to_mm(double(d.at.x)),
                          to_mm(double(d.at.y)), h.layer_a, h.layer_b, b.nets[h.net].c_str(), std::abs(h.current_a),
                          h.resistance_ohm * 1e3);
            s << line;
        }
    }
    if (!r.hot_spots.empty()) {
        s << "current density (largest)\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(5, r.hot_spots.size()); ++i) {
            const auto& h = r.hot_spots[i];
            std::snprintf(line, sizeof line, "  (%9.3f, %9.3f) mm  L%d  %-14s %9.3f A/mm2 (%.4f A/mm of width)\n", to_mm(h.x_nm),
                          to_mm(h.y_nm), h.layer, b.nets[h.net].c_str(), h.density_a_per_mm2, h.sheet_a_per_mm);
            s << line;
        }
    }
    for (const auto& w : r.warnings) s << "warning: " << w << "\n";
    return s.str();
}

}  // namespace kairoq
