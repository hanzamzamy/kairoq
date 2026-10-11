#include "kairoq/server.hpp"
#include "kairoq/thermal.hpp"

#include "kairoq/dc.hpp"
#include "kairoq/emi.hpp"
#include "kairoq/si.hpp"
#include "kairoq/embedded.hpp"
#include "kairoq/mesh.hpp"
#include "kairoq/model_io.hpp"

#include <httplib.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace kairoq {

namespace {

const char* mime_type(const std::string& path) {
    auto ends = [&](const char* s) { return path.size() >= std::strlen(s) && path.compare(path.size() - std::strlen(s), std::string::npos, s) == 0; };
    if (ends(".html")) return "text/html; charset=utf-8";
    if (ends(".js")) return "text/javascript; charset=utf-8";
    if (ends(".css")) return "text/css; charset=utf-8";
    if (ends(".json")) return "application/json";
    if (ends(".svg")) return "image/svg+xml";
    return "application/octet-stream";
}

void open_url(const std::string& url) {
#if defined(_WIN32)
    const std::string cmd = "start \"\" \"" + url + "\"";
#elif defined(__APPLE__)
    const std::string cmd = "open '" + url + "' >/dev/null 2>&1 &";
#else
    const std::string cmd = "xdg-open '" + url + "' >/dev/null 2>&1 &";
#endif
    if (std::system(cmd.c_str()) != 0) std::cerr << "could not open a browser; open " << url << "\n";
}

nlohmann::json display_json(const Board& board) {
    nlohmann::json j = board_json(board, false);
    const Point64 o = display_origin(board);
    j["display"]["origin_mm"] = {to_mm(double(o.x)), to_mm(double(o.y))};
    const auto z = stackup_z(board);
    nlohmann::json zs = nlohmann::json::array();
    for (const auto& lz : z) zs.push_back({to_mm(lz.bottom_nm), to_mm(lz.top_nm)});
    j["display"]["stackup_z_mm"] = zs;
    j["display"]["meshes"] = mesh_keys(board);
    double area = 0.0;
    for (const auto& p : board.outline) area += area_nm2(p);
    j["display"]["outline_area_mm2"] = area / 1.0e12;
    return j;
}

// DC solution of one copper layer, for display. Little endian:
// uint32 vertex count V, uint32 triangle count T, float32 positions [3V] (mm, relative
// to the display origin, z at the top of the copper), uint32 indices [3T], float32
// voltage [V] (V), float32 current density [T] (A/mm^2), float32 power density [T]
// (W/mm^2), uint32 net [T].
std::string encode_field(const Board& board, const DcLayerField& f) {
    const Point64 o = display_origin(board);
    const auto zs = stackup_z(board);
    double ztop = 0;
    for (std::size_t i = 0; i < board.stackup.size(); ++i)
        if (board.stackup[i].type == LayerType::Copper && board.stackup[i].copper_index == f.layer) ztop = zs[i].top_nm;
    std::string out;
    auto put32 = [&](std::uint32_t v) { out.append(reinterpret_cast<const char*>(&v), 4); };
    auto putf = [&](float v) { out.append(reinterpret_cast<const char*>(&v), 4); };
    const auto nv = static_cast<std::uint32_t>(f.x.size()), nt = static_cast<std::uint32_t>(f.triangles.size());
    put32(nv);
    put32(nt);
    for (std::size_t i = 0; i < nv; ++i) {
        putf(static_cast<float>(to_mm(f.x[i] - double(o.x))));
        putf(static_cast<float>(to_mm(f.y[i] - double(o.y))));
        putf(static_cast<float>(to_mm(ztop)));
    }
    for (const auto& t : f.triangles)
        for (auto v : t) put32(v);
    for (double v : f.voltage) putf(static_cast<float>(v));
    const double t_m = f.thickness_nm * 1e-9;
    for (std::size_t t = 0; t < nt; ++t) putf(static_cast<float>(std::hypot(f.jx[t], f.jy[t]) / t_m * 1e-6));
    for (std::size_t t = 0; t < nt; ++t) {
        const auto& tr = f.triangles[t];
        const double area_mm2 = 0.5 * std::abs((f.x[tr[1]] - f.x[tr[0]]) * (f.y[tr[2]] - f.y[tr[0]]) -
                                               (f.x[tr[2]] - f.x[tr[0]]) * (f.y[tr[1]] - f.y[tr[0]])) * 1e-12;
        putf(static_cast<float>(area_mm2 > 0 ? f.power_w[t] / area_mm2 : 0.0));
    }
    for (int n : f.net) put32(static_cast<std::uint32_t>(n));
    return out;
}

nlohmann::json dc_json(const Board& board, const DcResult& r) {
    nlohmann::json j;
    j["layers"] = nlohmann::json::array();
    double vmin = INFINITY, vmax = -INFINITY, jmax = 0, pmax = 0;
    for (const auto& f : r.layers) {
        j["layers"].push_back(f.layer);
        for (double v : f.voltage) {
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
        for (std::size_t t = 0; t < f.jx.size(); ++t) jmax = std::max(jmax, std::hypot(f.jx[t], f.jy[t]) / (f.thickness_nm * 1e-9) * 1e-6);
    }
    for (const auto& f : r.layers)
        for (std::size_t t = 0; t < f.triangles.size(); ++t) {
            const auto& tr = f.triangles[t];
            const double a = 0.5 * std::abs((f.x[tr[1]] - f.x[tr[0]]) * (f.y[tr[2]] - f.y[tr[0]]) -
                                            (f.x[tr[2]] - f.x[tr[0]]) * (f.y[tr[1]] - f.y[tr[0]])) * 1e-12;
            if (a > 0) pmax = std::max(pmax, f.power_w[t] / a);
        }
    j["range"] = {{"voltage", {vmin, vmax}}, {"current_density", {0.0, jmax}}, {"power_density", {0.0, pmax}}};
    j["ports"] = nlohmann::json::array();
    for (const auto& p : r.ports)
        j["ports"].push_back({{"name", p.name}, {"voltage_v", p.voltage_v}, {"current_a", p.current_a}, {"voltage_set", p.voltage_set}});
    j["power_w"] = r.power_ports_w;
    j["power_copper_w"] = r.power_copper_w;
    j["power_barrels_w"] = r.power_barrels_w;
    j["estimated_error"] = r.estimated_error;
    j["triangles"] = r.triangles;
    nlohmann::json holes = nlohmann::json::array();
    for (const auto& h : r.barrels)
        holes.push_back({{"drill", h.drill}, {"layers", {h.layer_a, h.layer_b}}, {"current_a", h.current_a}});
    j["barrels"] = holes;
    (void)board;
    return j;
}

// Thermal solution of one copper layer, little endian: uint32 V, uint32 T, float32
// positions [3V] (mm, display origin, z at the top of the copper), uint32 indices [3T],
// float32 temperature [V] (C), uint32 copper flag [T].
std::string encode_thermal(const Board& board, const ThermalLayerField& f) {
    const Point64 o = display_origin(board);
    const auto zs = stackup_z(board);
    double ztop = 0;
    for (std::size_t i = 0; i < board.stackup.size(); ++i)
        if (board.stackup[i].type == LayerType::Copper && board.stackup[i].copper_index == f.layer) ztop = zs[i].top_nm;
    std::string out;
    auto put32 = [&](std::uint32_t v) { out.append(reinterpret_cast<const char*>(&v), 4); };
    auto putf = [&](float v) { out.append(reinterpret_cast<const char*>(&v), 4); };
    put32(static_cast<std::uint32_t>(f.x.size()));
    put32(static_cast<std::uint32_t>(f.triangles.size()));
    for (std::size_t i = 0; i < f.x.size(); ++i) {
        putf(static_cast<float>(to_mm(f.x[i] - double(o.x))));
        putf(static_cast<float>(to_mm(f.y[i] - double(o.y))));
        putf(static_cast<float>(to_mm(ztop)));
    }
    for (const auto& t : f.triangles)
        for (auto v : t) put32(v);
    for (double t : f.temperature_c) putf(static_cast<float>(t));
    for (auto c : f.copper) put32(c);
    return out;
}

nlohmann::json thermal_json(const ThermalResult& r) {
    nlohmann::json j;
    j["layers"] = nlohmann::json::array();
    double lo = INFINITY, hi = -INFINITY;
    for (const auto& f : r.layers) {
        j["layers"].push_back(f.layer);
        for (double t : f.temperature_c) {
            lo = std::min(lo, t);
            hi = std::max(hi, t);
        }
    }
    if (r.transient)
        for (const auto& snap : r.transient->snapshot_c)
            for (const auto& layer : snap)
                for (double t : layer) {
                    lo = std::min(lo, t);
                    hi = std::max(hi, t);
                }
    j["range"] = {lo, hi};
    j["max_c"] = r.max_c;
    j["max_position_mm"] = {to_mm(r.max_x_nm), to_mm(r.max_y_nm)};
    j["max_layer"] = r.max_layer;
    j["estimated_error"] = r.estimated_error;
    j["heat_in_w"] = r.heat_sources_w + r.heat_joule_w;
    j["sources"] = nlohmann::json::array();
    for (const auto& s : r.sources) {
        nlohmann::json e{{"name", s.name}, {"power_w", s.power_w}, {"mean_c", s.mean_c}, {"max_c", s.max_c}, {"fixed", s.temperature_set}};
        if (s.junction_c) e["junction_c"] = *s.junction_c;
        j["sources"].push_back(e);
    }
    j["heatsinks"] = nlohmann::json::array();
    for (const auto& h : r.heatsinks)
        j["heatsinks"].push_back({{"name", h.name}, {"heat_w", h.heat_w}, {"temperature_c", h.temperature_c}, {"fixed", h.temperature_set}});
    if (r.transient) {
        const auto& t = *r.transient;
        nlohmann::json tr;
        tr["time_s"] = t.time_s;
        tr["max_c"] = t.max_c;
        tr["snapshot_time_s"] = t.snapshot_time_s;
        tr["series"] = nlohmann::json::array();
        for (std::size_t a = 0; a < r.sources.size(); ++a) {
            if (r.sources[a].temperature_set) continue;
            const bool junction = r.sources[a].junction_c.has_value();
            tr["series"].push_back({{"name", r.sources[a].name + (junction ? " junction" : "")},
                                    {"values", junction ? t.source_junction_c[a] : t.source_max_c[a]}});
        }
        for (std::size_t k = 0; k < r.heatsinks.size(); ++k)
            if (!r.heatsinks[k].temperature_set) tr["series"].push_back({{"name", r.heatsinks[k].name}, {"values", t.heatsink_c[k]}});
        j["transient"] = tr;
    }
    return j;
}

// Sheet current of a loop on one layer, little endian: uint32 V, uint32 T, float32 positions
// [3V] (mm, display origin, z at the top of the copper), uint32 indices [3T], float32 sheet
// current per ampere of loop current [T] (A/mm).
std::string encode_loop(const Board& board, const LoopSheetField& f) {
    const Point64 o = display_origin(board);
    const auto zs = stackup_z(board);
    double ztop = 0;
    for (std::size_t i = 0; i < board.stackup.size(); ++i)
        if (board.stackup[i].type == LayerType::Copper && board.stackup[i].copper_index == f.layer) ztop = zs[i].top_nm;
    std::string out;
    auto put32 = [&](std::uint32_t v) { out.append(reinterpret_cast<const char*>(&v), 4); };
    auto putf = [&](float v) { out.append(reinterpret_cast<const char*>(&v), 4); };
    put32(static_cast<std::uint32_t>(f.x.size()));
    put32(static_cast<std::uint32_t>(f.triangles.size()));
    for (std::size_t i = 0; i < f.x.size(); ++i) {
        putf(static_cast<float>(to_mm(f.x[i] - double(o.x))));
        putf(static_cast<float>(to_mm(f.y[i] - double(o.y))));
        putf(static_cast<float>(to_mm(ztop)));
    }
    for (const auto& t : f.triangles)
        for (auto v : t) put32(v);
    for (const auto& k : f.current_a_per_m) putf(static_cast<float>(std::hypot(k[0], k[1]) * 1e-3));
    return out;
}

nlohmann::json emi_json(const Board& board, const EmiResult& r) {
    const Point64 o = display_origin(board);
    const auto zs = stackup_z(board);
    auto ztop = [&](int layer) {
        for (std::size_t i = 0; i < board.stackup.size(); ++i)
            if (board.stackup[i].type == LayerType::Copper && board.stackup[i].copper_index == layer) return to_mm(zs[i].top_nm);
        return 0.0;
    };
    nlohmann::json j;
    j["loops"] = nlohmann::json::array();
    for (const auto& l : r.loops) {
        double kmax = 0;
        nlohmann::json layers = nlohmann::json::array();
        for (const auto& sh : l.sheets) {
            layers.push_back(sh.layer);
            for (const auto& k : sh.current_a_per_m) kmax = std::max(kmax, std::hypot(k[0], k[1]) * 1e-3);
        }
        const double m = std::sqrt(l.moment_m2[0] * l.moment_m2[0] + l.moment_m2[1] * l.moment_m2[1] + l.moment_m2[2] * l.moment_m2[2]);
        j["loops"].push_back({{"name", l.name}, {"inductance_nh", l.inductance_nh}, {"inductance_coarse_nh", l.inductance_coarse_nh},
                              {"moment_mm2", m * 1e6}, {"layers", layers}, {"max_current_a_per_mm", kmax},
                              {"max_field_dbuv_m", l.max_field_v_per_m > 0 ? 20 * std::log10(l.max_field_v_per_m / 1e-6) : -999.0},
                              {"max_field_frequency_hz", l.max_field_frequency_hz}});
    }
    j["findings"] = nlohmann::json::array();
    for (const auto& f : r.return_path) {
        nlohmann::json outline = nlohmann::json::array();
        if (!f.piece.empty())
            for (const auto& p : f.piece.front().outer) outline.push_back({to_mm(double(p.x - o.x)), to_mm(double(p.y - o.y))});
        j["findings"].push_back({{"net", f.net}, {"layer", f.layer}, {"reference_layer", f.reference_layer}, {"kind", f.kind},
                                 {"reference_nets", f.reference_nets}, {"position_mm", {to_mm(f.x_nm), to_mm(f.y_nm)}},
                                 {"span_mm", f.span_mm}, {"area_mm2", f.area_mm2}, {"z_mm", ztop(f.layer)}, {"outline", outline}});
    }
    j["coupling"] = nlohmann::json::array();
    for (const auto& c : r.coupling)
        j["coupling"].push_back({{"net", c.net}, {"plate", c.plate}, {"capacitance_pf", c.capacitance_pf},
                                 {"common_mode_current_ma", c.common_mode_current_ma}});
    return j;
}

}  // namespace

void serve(const Board& board, const ServeOptions& options, const DcResult* dc, const ThermalResult* thermal, const EmiResult* emi,
           const SiResult* si) {
    httplib::Server server;
    nlohmann::json shown = display_json(board);
    if (dc) shown["dc"] = dc_json(board, *dc);
    if (thermal) shown["thermal"] = thermal_json(*thermal);
    if (emi) shown["emi"] = emi_json(board, *emi);
    if (si) shown["si"] = nlohmann::json::parse(si_json(board, *si));
    const std::string board_text = shown.dump();
    std::map<std::string, std::string> cache;
    std::mutex cache_mutex;

    server.Get("/api/board", [&](const httplib::Request&, httplib::Response& res) {
        res.set_content(board_text, "application/json");
    });
    server.Get("/api/mesh", [&](const httplib::Request& req, httplib::Response& res) {
        const std::string key = req.get_param_value("key");
        try {
            std::string body;
            {
                std::lock_guard<std::mutex> lock(cache_mutex);
                auto it = cache.find(key);
                if (it != cache.end()) body = it->second;
            }
            if (body.empty()) {
                body = encode_mesh(board_mesh(board, key));
                std::lock_guard<std::mutex> lock(cache_mutex);
                cache[key] = body;
            }
            res.set_content(body, "application/octet-stream");
        } catch (const std::exception& e) {
            res.status = 404;
            res.set_content(e.what(), "text/plain");
        }
    });
    server.Get("/api/field", [&](const httplib::Request& req, httplib::Response& res) {
        const int layer = std::atoi(req.get_param_value("layer").c_str());
        if (dc)
            for (const auto& f : dc->layers)
                if (f.layer == layer) {
                    res.set_content(encode_field(board, f), "application/octet-stream");
                    return;
                }
        res.status = 404;
        res.set_content("no field on this layer", "text/plain");
    });
    server.Get("/api/thermal", [&](const httplib::Request& req, httplib::Response& res) {
        const int layer = std::atoi(req.get_param_value("layer").c_str());
        if (thermal)
            for (const auto& f : thermal->layers)
                if (f.layer == layer) {
                    res.set_content(encode_thermal(board, f), "application/octet-stream");
                    return;
                }
        res.status = 404;
        res.set_content("no thermal solution on this layer", "text/plain");
    });
    server.Get("/api/loop", [&](const httplib::Request& req, httplib::Response& res) {
        const long loop = std::atol(req.get_param_value("loop").c_str());
        const int layer = std::atoi(req.get_param_value("layer").c_str());
        if (emi && loop >= 0 && static_cast<std::size_t>(loop) < emi->loops.size())
            for (const auto& f : emi->loops[static_cast<std::size_t>(loop)].sheets)
                if (f.layer == layer) {
                    res.set_content(encode_loop(board, f), "application/octet-stream");
                    return;
                }
        res.status = 404;
        res.set_content("no loop current on this layer", "text/plain");
    });
    // Temperatures of a transient snapshot on a layer's thermal mesh: float32 [V].
    server.Get("/api/thermal_values", [&](const httplib::Request& req, httplib::Response& res) {
        const int layer = std::atoi(req.get_param_value("layer").c_str());
        const long snap = std::atol(req.get_param_value("snapshot").c_str());
        if (thermal && thermal->transient && snap >= 0 && static_cast<std::size_t>(snap) < thermal->transient->snapshot_c.size())
            for (std::size_t l = 0; l < thermal->layers.size(); ++l)
                if (thermal->layers[l].layer == layer) {
                    const auto& v = thermal->transient->snapshot_c[static_cast<std::size_t>(snap)][l];
                    std::string out(v.size() * 4, '\0');
                    for (std::size_t i = 0; i < v.size(); ++i) {
                        const float f = static_cast<float>(v[i]);
                        std::memcpy(&out[4 * i], &f, 4);
                    }
                    res.set_content(out, "application/octet-stream");
                    return;
                }
        res.status = 404;
        res.set_content("no such snapshot", "text/plain");
    });
    server.Get(R"(/(.*))", [&](const httplib::Request& req, httplib::Response& res) {
        std::string path = req.matches[1];
        if (path.empty()) path = "index.html";
        for (std::size_t i = 0; i < viewer_files_count; ++i) {
            if (path == viewer_files[i].path) {
                res.set_content(reinterpret_cast<const char*>(viewer_files[i].data), viewer_files[i].size, mime_type(path));
                return;
            }
        }
        res.status = 404;
        res.set_content("not found", "text/plain");
    });

    int port = options.port;
    if (port == 0) {
        port = server.bind_to_any_port(options.host);
    } else if (!server.bind_to_port(options.host, port)) {
        throw std::runtime_error("cannot listen on " + options.host + ":" + std::to_string(port));
    }
    const std::string url = "http://" + options.host + ":" + std::to_string(port) + "/";
    std::cout << "viewer at " << url << " (Ctrl+C to stop)" << std::endl;
    if (options.open_browser) open_url(url);
    if (!options.should_stop) {
        server.listen_after_bind();
        return;
    }
    std::thread listener([&server]() { server.listen_after_bind(); });
    server.wait_until_ready();
    while (server.is_running()) {
        if (options.should_stop()) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    server.stop();
    listener.join();
}

}  // namespace kairoq
