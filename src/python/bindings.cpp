// Python module kairoq._core. The kairoq package wraps it (src/python/kairoq).
#include "kairoq/cross_section.hpp"
#include "kairoq/dc.hpp"
#include "kairoq/emi.hpp"
#include "kairoq/deembed.hpp"
#include "kairoq/eye.hpp"
#include "kairoq/rf.hpp"
#include "kairoq/ibis.hpp"
#include "kairoq/mesh.hpp"
#include "kairoq/model_io.hpp"
#include "kairoq/network.hpp"
#include "kairoq/open.hpp"
#include "kairoq/openems.hpp"
#include "kairoq/pdn.hpp"
#include "kairoq/provenance.hpp"
#include "kairoq/server.hpp"
#include "kairoq/si.hpp"
#include "kairoq/thermal.hpp"
#include "kairoq/timedomain.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <nlohmann/json.hpp>

#include <Python.h>

#include <map>

namespace nb = nanobind;
using namespace kairoq;

namespace {

nb::bytes to_bytes(const std::string& s) { return nb::bytes(s.data(), s.size()); }

ImportOptions options(double arc_tolerance_um, const std::optional<std::string>& stackup,
                      const std::optional<std::array<double, 4>>& crop_mm = std::nullopt) {
    ImportOptions o;
    o.geometry.arc_tolerance_nm = arc_tolerance_um * 1000.0;
    if (stackup) o.stackup_yaml_path = *stackup;
    if (crop_mm) o.crop = Rect64(nm_from_mm((*crop_mm)[0]), nm_from_mm((*crop_mm)[1]), nm_from_mm((*crop_mm)[2]), nm_from_mm((*crop_mm)[3]));
    return o;
}

const CopperLayer& layer_of(const Board& b, int index) {
    for (const auto& l : b.copper)
        if (l.index == index) return l;
    throw nb::index_error(("no copper layer " + std::to_string(index)).c_str());
}

template <typename T>
nb::ndarray<nb::numpy, T> array(std::vector<T> v, std::initializer_list<std::size_t> shape) {
    auto* heap = new std::vector<T>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept { delete static_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, T>(heap->data(), shape, owner);
}

}  // namespace

NB_MODULE(_core, m) {
    m.doc() = "kairoq core library";

    nb::class_<Board>(m, "Board")
        .def("json", [](const Board& b) { return board_json(b, false).dump(); }, "Board data without geometry, as JSON.")
        .def("summary", &board_summary)
        .def("save", &save_board, nb::arg("path"), "Write the canonical model (directory, or .zip archive).")
        .def(
            "copper_wkb",
            [](const Board& b, int index) {
                nb::list out;
                for (const auto& nc : layer_of(b, index).nets) out.append(nb::make_tuple(nc.net, to_bytes(encode_wkb(nc.polygons))));
                return out;
            },
            nb::arg("layer"), "[(net id, WKB MultiPolygon in mm)] of a copper layer.")
        .def("outline_wkb", [](const Board& b) { return to_bytes(encode_wkb(b.outline)); })
        .def(
            "mask_wkb",
            [](const Board& b, const std::string& side) -> std::optional<nb::bytes> {
                for (const auto& m : b.masks)
                    if (to_string(m.side) == side) return to_bytes(encode_wkb(m.openings));
                return std::nullopt;
            },
            nb::arg("side"))
        .def(
            "pad_wkb", [](const Board& b, std::size_t i) { return to_bytes(encode_wkb(b.pads.at(i).shape)); }, nb::arg("index"))
        .def(
            "copper_area_mm2",
            [](const Board& b, int index, std::optional<int> net) {
                double a = 0;
                for (const auto& nc : layer_of(b, index).nets)
                    if (!net || nc.net == *net) a += area_nm2(nc.polygons);
                return a / 1e12;
            },
            nb::arg("layer"), nb::arg("net") = nb::none())
        .def("mesh_keys", &mesh_keys)
        .def(
            "mesh",
            [](const Board& b, const std::string& key) {
                Mesh mesh = board_mesh(b, key);
                const std::size_t nv = mesh.ids.size(), nt = mesh.indices.size() / 3;
                nb::dict d;
                d["positions"] = array(std::move(mesh.positions), {nv, 3});
                d["triangles"] = array(std::move(mesh.indices), {nt, 3});
                d["ids"] = array(std::move(mesh.ids), {nv});
                const Point64 o = display_origin(b);
                d["origin_mm"] = nb::make_tuple(to_mm(double(o.x)), to_mm(double(o.y)));
                return d;
            },
            nb::arg("key"), "Display mesh: positions (mm, relative to origin_mm), triangles, per-vertex ids.");

    nb::class_<DcResult>(m, "DcResult")
        .def_ro("power_copper_w", &DcResult::power_copper_w)
        .def_ro("power_barrels_w", &DcResult::power_barrels_w)
        .def_ro("power_ports_w", &DcResult::power_ports_w)
        .def_ro("unknowns", &DcResult::unknowns)
        .def_ro("triangles", &DcResult::triangles)
        .def_ro("islands", &DcResult::islands)
        .def_ro("iterations", &DcResult::iterations)
        .def_ro("residual", &DcResult::residual)
        .def_ro("seconds_mesh", &DcResult::seconds_mesh)
        .def_ro("seconds_assemble", &DcResult::seconds_assemble)
        .def_ro("seconds_solve", &DcResult::seconds_solve)
        .def_ro("warnings", &DcResult::warnings)
        .def_ro("estimated_error", &DcResult::estimated_error)
        .def_prop_ro("passes",
                     [](const DcResult& r) {
                         nb::list out;
                         for (const auto& p : r.passes) {
                             nb::dict d;
                             d["triangles"] = p.triangles;
                             d["estimated_error"] = p.estimated_error;
                             d["power_w"] = p.power_w;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("ports",
                     [](const DcResult& r) {
                         nb::list out;
                         for (const auto& p : r.ports) {
                             nb::dict d;
                             d["name"] = p.name;
                             d["voltage_v"] = p.voltage_v;
                             d["current_a"] = p.current_a;
                             d["voltage_set"] = p.voltage_set;
                             d["nodes"] = p.nodes;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("barrels",
                     [](const DcResult& r) {
                         nb::list out;
                         for (const auto& h : r.barrels) {
                             nb::dict d;
                             d["drill"] = h.drill;
                             d["net"] = h.net;
                             d["layer_a"] = h.layer_a;
                             d["layer_b"] = h.layer_b;
                             d["resistance_ohm"] = h.resistance_ohm;
                             d["current_a"] = h.current_a;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("hot_spots",
                     [](const DcResult& r) {
                         nb::list out;
                         for (const auto& h : r.hot_spots) {
                             nb::dict d;
                             d["layer"] = h.layer;
                             d["net"] = h.net;
                             d["position_mm"] = nb::make_tuple(to_mm(h.x_nm), to_mm(h.y_nm));
                             d["sheet_a_per_mm"] = h.sheet_a_per_mm;
                             d["density_a_per_mm2"] = h.density_a_per_mm2;
                             out.append(d);
                         }
                         return out;
                     })
        .def(
            "layer",
            [](const DcResult& r, int layer) {
                for (const auto& f : r.layers) {
                    if (f.layer != layer) continue;
                    const std::size_t nv = f.x.size(), nt = f.triangles.size();
                    std::vector<double> xy(2 * nv);
                    for (std::size_t i = 0; i < nv; ++i) {
                        xy[2 * i] = to_mm(f.x[i]);
                        xy[2 * i + 1] = to_mm(f.y[i]);
                    }
                    std::vector<std::uint32_t> tri(3 * nt);
                    for (std::size_t t = 0; t < nt; ++t)
                        for (int k = 0; k < 3; ++k) tri[3 * t + k] = f.triangles[t][k];
                    std::vector<double> j(2 * nt);
                    for (std::size_t t = 0; t < nt; ++t) {
                        j[2 * t] = f.jx[t] * 1e-3;
                        j[2 * t + 1] = f.jy[t] * 1e-3;
                    }
                    nb::dict d;
                    d["layer"] = f.layer;
                    d["thickness_mm"] = to_mm(f.thickness_nm);
                    d["conductivity_s_per_m"] = f.conductivity_s_per_m;
                    d["points_mm"] = array(std::move(xy), {nv, 2});
                    d["triangles"] = array(std::move(tri), {nt, 3});
                    d["voltage_v"] = array(std::vector<double>(f.voltage), {nv});
                    d["sheet_current_a_per_mm"] = array(std::move(j), {nt, 2});
                    d["power_w"] = array(std::vector<double>(f.power_w), {nt});
                    d["net"] = array(std::vector<int>(f.net), {nt});
                    return d;
                }
                throw nb::index_error(("no solution on layer " + std::to_string(layer)).c_str());
            },
            nb::arg("layer"),
            "Solution on a copper layer: points (mm), triangles, nodal voltage, per-triangle sheet current density "
            "(A/mm), Joule power (W) and net id.")
        .def_prop_ro("layers", [](const DcResult& r) {
            std::vector<int> out;
            for (const auto& f : r.layers) out.push_back(f.layer);
            return out;
        });

    nb::class_<ThermalResult>(m, "ThermalResult")
        .def_ro("max_c", &ThermalResult::max_c)
        .def_ro("max_layer", &ThermalResult::max_layer)
        .def_prop_ro("max_position_mm", [](const ThermalResult& r) { return nb::make_tuple(to_mm(r.max_x_nm), to_mm(r.max_y_nm)); })
        .def_ro("heat_sources_w", &ThermalResult::heat_sources_w)
        .def_ro("heat_joule_w", &ThermalResult::heat_joule_w)
        .def_ro("heat_top_w", &ThermalResult::heat_top_w)
        .def_ro("heat_bottom_w", &ThermalResult::heat_bottom_w)
        .def_ro("heat_fixed_w", &ThermalResult::heat_fixed_w)
        .def_ro("heat_sinks_w", &ThermalResult::heat_sinks_w)
        .def_prop_ro("heatsinks",
                     [](const ThermalResult& r) {
                         nb::list out;
                         for (const auto& h : r.heatsinks) {
                             nb::dict d;
                             d["name"] = h.name;
                             d["heat_w"] = h.heat_w;
                             d["temperature_c"] = h.temperature_c;
                             d["temperature_set"] = h.temperature_set;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("transient",
                     [](const ThermalResult& r) -> nb::object {
                         if (!r.transient) return nb::none();
                         const auto& t = *r.transient;
                         const std::size_t n = t.time_s.size();
                         auto vec = [&](const std::vector<double>& v) { return array(std::vector<double>(v), {v.size()}); };
                         auto per = [&](const std::vector<std::vector<double>>& v) {
                             std::vector<double> flat;
                             for (const auto& row : v) flat.insert(flat.end(), row.begin(), row.end());
                             return array(std::move(flat), {v.size(), n});
                         };
                         nb::dict d;
                         d["time_s"] = vec(t.time_s);
                         d["max_c"] = vec(t.max_c);
                         d["source_mean_c"] = per(t.source_mean_c);
                         d["source_max_c"] = per(t.source_max_c);
                         d["source_junction_c"] = per(t.source_junction_c);
                         d["heatsink_c"] = per(t.heatsink_c);
                         d["heatsink_w"] = per(t.heatsink_w);
                         d["heat_in_w"] = vec(t.heat_in_w);
                         d["heat_out_w"] = vec(t.heat_out_w);
                         d["heat_stored_w"] = vec(t.heat_stored_w);
                         d["snapshot_time_s"] = vec(t.snapshot_time_s);
                         d["steps"] = t.steps;
                         d["matrix_setups"] = t.matrix_setups;
                         d["solver_iterations"] = t.solver_iterations;
                         d["dc_solutions"] = t.dc_solutions;
                         return d;
                     })
        .def(
            "snapshot",
            [](const ThermalResult& r, std::size_t index, int layer) {
                if (!r.transient || index >= r.transient->snapshot_c.size()) throw nb::index_error("no such snapshot");
                for (std::size_t l = 0; l < r.layers.size(); ++l)
                    if (r.layers[l].layer == layer) {
                        const auto& v = r.transient->snapshot_c[index][l];
                        return array(std::vector<double>(v), {v.size()});
                    }
                throw nb::index_error(("no solution on layer " + std::to_string(layer)).c_str());
            },
            nb::arg("index"), nb::arg("layer"))
        .def_ro("triangles", &ThermalResult::triangles)
        .def_ro("unknowns", &ThermalResult::unknowns)
        .def_ro("iterations", &ThermalResult::iterations)
        .def_ro("estimated_error", &ThermalResult::estimated_error)
        .def_ro("coupling_iterations", &ThermalResult::coupling_iterations)
        .def_ro("coupling_change_k", &ThermalResult::coupling_change_k)
        .def_ro("seconds_mesh", &ThermalResult::seconds_mesh)
        .def_ro("seconds_assemble", &ThermalResult::seconds_assemble)
        .def_ro("seconds_solve", &ThermalResult::seconds_solve)
        .def_ro("seconds_dc", &ThermalResult::seconds_dc)
        .def_ro("warnings", &ThermalResult::warnings)
        .def_prop_ro("passes",
                     [](const ThermalResult& r) {
                         nb::list out;
                         for (const auto& p : r.passes) {
                             nb::dict d;
                             d["triangles"] = p.triangles;
                             d["estimated_error"] = p.estimated_error;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("sources",
                     [](const ThermalResult& r) {
                         nb::list out;
                         for (const auto& s : r.sources) {
                             nb::dict d;
                             d["name"] = s.name;
                             d["power_w"] = s.power_w;
                             d["temperature_set"] = s.temperature_set;
                             d["mean_c"] = s.mean_c;
                             d["max_c"] = s.max_c;
                             d["junction_c"] = s.junction_c ? nb::cast(*s.junction_c) : nb::none();
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("dc", [](const ThermalResult& r) -> std::optional<DcResult> { return r.dc; })
        .def_prop_ro("layers",
                     [](const ThermalResult& r) {
                         std::vector<int> out;
                         for (const auto& f : r.layers) out.push_back(f.layer);
                         return out;
                     })
        .def(
            "layer",
            [](const ThermalResult& r, int layer) {
                for (const auto& f : r.layers) {
                    if (f.layer != layer) continue;
                    const std::size_t nv = f.x.size(), nt = f.triangles.size();
                    std::vector<double> xy(2 * nv);
                    for (std::size_t i = 0; i < nv; ++i) {
                        xy[2 * i] = to_mm(f.x[i]);
                        xy[2 * i + 1] = to_mm(f.y[i]);
                    }
                    std::vector<std::uint32_t> tri(3 * nt);
                    for (std::size_t t = 0; t < nt; ++t)
                        for (int k = 0; k < 3; ++k) tri[3 * t + k] = f.triangles[t][k];
                    nb::dict d;
                    d["layer"] = f.layer;
                    d["points_mm"] = array(std::move(xy), {nv, 2});
                    d["triangles"] = array(std::move(tri), {nt, 3});
                    d["temperature_c"] = array(std::vector<double>(f.temperature_c), {nv});
                    d["copper"] = array(std::vector<std::uint8_t>(f.copper), {nt});
                    return d;
                }
                throw nb::index_error(("no solution on layer " + std::to_string(layer)).c_str());
            },
            nb::arg("layer"));

    m.def(
        "solve_thermal",
        [](const Board& b, const std::string& yaml, const std::string& base_dir) {
            ThermalProject p = parse_thermal_project(yaml, base_dir);
            if (!p.thermal) throw std::runtime_error("no thermal section");
            std::vector<std::string> warnings;
            const ThermalSettings s = resolve_thermal_project(b, p, warnings);
            nb::gil_scoped_release release;
            ThermalResult r = solve_thermal(b, s);
            r.warnings.insert(r.warnings.begin(), warnings.begin(), warnings.end());
            return r;
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def("thermal_summary", &thermal_summary);

    nb::class_<EmiResult>(m, "EmiResult")
        .def_ro("warnings", &EmiResult::warnings)
        .def_prop_ro("loops",
                     [](const EmiResult& r) {
                         nb::list out;
                         for (const auto& l : r.loops) {
                             nb::dict d;
                             d["name"] = l.name;
                             d["inductance_nh"] = l.inductance_nh;
                             d["inductance_coarse_nh"] = l.inductance_coarse_nh;
                             d["element_size_mm"] = l.element_size_mm;
                             d["memory_mb"] = l.memory_mb;
                             d["warnings"] = l.warnings;
                             d["element_inductance_nh"] = l.element_inductance_nh;
                             d["moment_m2"] = nb::make_tuple(l.moment_m2[0], l.moment_m2[1], l.moment_m2[2]);
                             d["window_margin_mm"] = l.window_margin_mm;
                             d["moment_change"] = l.moment_change;
                             d["triangles"] = l.triangles;
                             d["unknowns"] = l.unknowns;
                             d["filaments"] = l.filaments;
                             d["seconds"] = l.seconds;
                             d["max_field_v_per_m"] = l.max_field_v_per_m;
                             d["max_field_frequency_hz"] = l.max_field_frequency_hz;
                             nb::list imp;
                             for (const auto& z : l.impedance) {
                                 nb::list sc;
                                 for (const auto& c : z.stitch_current) sc.append(c);
                                 imp.append(nb::make_tuple(z.frequency_hz, z.resistance_ohm, z.inductance_nh, sc));
                             }
                             d["impedance"] = imp;
                             d["stitch_current"] = l.stitch_current;
                             nb::list fpa;
                             for (const auto& [f, e] : l.field_per_ampere) fpa.append(nb::make_tuple(f, e));
                             d["field_per_ampere"] = fpa;
                             {
                                 const std::size_t nc = l.currents.size();
                                 std::vector<double> pos(3 * nc), mom(3 * nc);
                                 for (std::size_t i = 0; i < nc; ++i)
                                     for (int c = 0; c < 3; ++c) {
                                         pos[3 * i + c] = l.currents[i].position_m[c];
                                         mom[3 * i + c] = l.currents[i].moment_a_m[c];
                                     }
                                 d["current_positions_m"] = array(std::move(pos), {nc, 3});
                                 d["current_moments_a_m"] = array(std::move(mom), {nc, 3});
                             }
                             std::vector<double> f, a, e;
                             for (const auto& h : l.harmonics) {
                                 f.push_back(h.frequency_hz);
                                 a.push_back(h.current_a);
                                 e.push_back(h.field_v_per_m);
                             }
                             const std::size_t nh = f.size();
                             d["harmonic_frequency_hz"] = array(std::move(f), {nh});
                             d["harmonic_current_a"] = array(std::move(a), {nh});
                             d["harmonic_field_v_per_m"] = array(std::move(e), {nh});
                             nb::list sheets;
                             for (const auto& sh : l.sheets) {
                                 nb::dict q;
                                 const std::size_t nv = sh.x.size(), nt = sh.triangles.size();
                                 std::vector<double> xy(2 * nv), k(2 * nt);
                                 for (std::size_t i = 0; i < nv; ++i) {
                                     xy[2 * i] = to_mm(sh.x[i]);
                                     xy[2 * i + 1] = to_mm(sh.y[i]);
                                 }
                                 std::vector<std::uint32_t> tri(3 * nt);
                                 for (std::size_t t = 0; t < nt; ++t) {
                                     for (int c = 0; c < 3; ++c) tri[3 * t + c] = sh.triangles[t][c];
                                     k[2 * t] = sh.current_a_per_m[t][0];
                                     k[2 * t + 1] = sh.current_a_per_m[t][1];
                                 }
                                 q["layer"] = sh.layer;
                                 q["points_mm"] = array(std::move(xy), {nv, 2});
                                 q["triangles"] = array(std::move(tri), {nt, 3});
                                 q["current_a_per_m"] = array(std::move(k), {nt, 2});
                                 sheets.append(q);
                             }
                             d["sheets"] = sheets;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("return_path",
                     [](const EmiResult& r) {
                         nb::list out;
                         for (const auto& f : r.return_path) {
                             nb::dict d;
                             d["net"] = f.net;
                             d["layer"] = f.layer;
                             d["reference_layer"] = f.reference_layer;
                             d["kind"] = f.kind;
                             d["reference_nets"] = f.reference_nets;
                             d["position_mm"] = nb::make_tuple(to_mm(f.x_nm), to_mm(f.y_nm));
                             d["area_mm2"] = f.area_mm2;
                             d["span_mm"] = f.span_mm;
                             out.append(d);
                         }
                         return out;
                     })
        .def_prop_ro("coupling", [](const EmiResult& r) {
            nb::list out;
            for (const auto& c : r.coupling) {
                nb::dict d;
                d["net"] = c.net;
                d["plate"] = c.plate;
                d["capacitance_pf"] = c.capacitance_pf;
                d["parallel_plate_pf"] = c.parallel_plate_pf;
                d["total_pf"] = c.total_pf;
                d["mesh_pf"] = c.mesh_pf;
                d["coarse_pf"] = c.coarse_pf;
                nb::list nets;
                for (const auto& [net, pf] : c.nets_pf) nets.append(nb::make_tuple(net, pf));
                d["nets_pf"] = nets;
                d["panels"] = c.panels;
                d["window_margin_mm"] = c.window_margin_mm;
                d["unshielded_area_mm2"] = c.unshielded_area_mm2;
                d["shielded_area_mm2"] = c.shielded_area_mm2;
                d["common_mode_current_ma"] = c.common_mode_current_ma;
                nb::list layers;
                for (const auto& [l, pf] : c.layers) layers.append(nb::make_tuple(l, pf));
                d["layers"] = layers;
                nb::list nbs;
                for (const auto& [net, l, pf] : c.neighbours) nbs.append(nb::make_tuple(net, l, pf));
                d["neighbours"] = nbs;
                out.append(d);
            }
            return out;
        });
    m.def(
        "analyse_emi",
        [](const Board& b, const std::string& yaml, const std::string& base_dir) {
            EmiProject p = parse_emi_project(yaml, base_dir);
            if (!p.emi) throw std::runtime_error("no emi section");
            nb::gil_scoped_release release;
            return analyse_emi(b, *p.emi);
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def("emi_summary", &emi_summary);
    nb::class_<SiResult>(m, "SiResult")
        .def("json", [](const SiResult& r, const Board& b) { return si_json(b, r); })
        .def("summary", [](const SiResult& r, const Board& b) { return si_summary(b, r); });
    m.def(
        "analyse_si",
        [](const Board& b, const std::string& yaml, const std::string& base_dir) {
            SiProject p = parse_si_project(yaml, base_dir);
            if (!p.si) throw std::runtime_error("no si section");
            nb::gil_scoped_release release;
            return analyse_si(b, *p.si);
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def(
        "analyse_channels",
        [](const Board& b, const std::string& yaml, const std::string& base_dir) {
            SiProject p = parse_si_project(yaml, base_dir);
            if (!p.si) throw std::runtime_error("no si section");
            const ChannelProject cp = parse_channel_project(yaml);
            if (cp.channels.empty()) throw std::runtime_error("no si.channels");
            std::vector<ChannelResult> r;
            {
                nb::gil_scoped_release release;
                r = analyse_channels(b, *p.si, cp.channels, cp.network, base_dir);
            }
            std::vector<std::string> ts, tm;
            for (const auto& c : r) {
                ts.push_back(touchstone(c, false));
                tm.push_back(c.mixed.empty() ? std::string() : touchstone(c, true));
            }
            return nb::make_tuple(channels_json(r), channels_summary(r), ts, tm);
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def(
        "simulate",
        [](const Board& b, const std::string& yaml, const std::string& base_dir) {
            SiProject p = parse_si_project(yaml, base_dir);
            if (!p.si) throw std::runtime_error("no si section");
            const ChannelProject cp = parse_channel_project(yaml);
            const auto sims = parse_simulations(yaml);
            if (sims.empty()) throw std::runtime_error("no si.simulations");
            std::vector<SimulationResult> out;
            {
                nb::gil_scoped_release release;
                std::map<std::string, ChannelResult> channels;
                for (const auto& sp : sims) {
                    if (!channels.count(sp.channel)) {
                        const ChannelSpec* c = nullptr;
                        for (const auto& x : cp.channels)
                            if (x.name == sp.channel) c = &x;
                        if (!c) throw std::runtime_error("simulation " + sp.name + ": no channel " + sp.channel);
                        channels[sp.channel] = analyse_channels(b, *p.si, {*c}, cp.network, base_dir).front();
                    }
                    out.push_back(simulate(channels[sp.channel], sp, base_dir));
                }
            }
            return std::make_pair(simulations_json(out), simulations_summary(out));
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def(
        "pdn",
        [](const Board& b, const std::string& yaml) {
            const PdnSettings settings = parse_pdn(yaml);
            PdnResult r;
            {
                nb::gil_scoped_release release;
                r = analyse_pdn(b, settings);
            }
            nlohmann::json j = nlohmann::json::parse(pdn_json(r));
            j["summary"] = pdn_summary(r, settings);
            return j.dump();
        },
        nb::arg("board"), nb::arg("yaml"));
    m.def(
        "eye",
        [](const Board& b, const std::string& yaml, const std::string& base_dir) {
            SiProject p = parse_si_project(yaml, base_dir);
            if (!p.si) throw std::runtime_error("no si section");
            const ChannelProject cp = parse_channel_project(yaml);
            const auto eyes = parse_eyes(yaml);
            if (eyes.empty()) throw std::runtime_error("no si.eyes");
            std::vector<EyeResult> out;
            std::vector<std::string> svg;
            {
                nb::gil_scoped_release release;
                std::map<std::string, ChannelResult> channels;
                for (const auto& e : eyes) {
                    if (!channels.count(e.channel)) {
                        const ChannelSpec* c = nullptr;
                        for (const auto& x : cp.channels)
                            if (x.name == e.channel) c = &x;
                        if (!c) throw std::runtime_error("eye " + e.name + ": no channel " + e.channel);
                        channels[e.channel] = analyse_channels(b, *p.si, {*c}, cp.network, base_dir).front();
                    }
                    out.push_back(eye(channels[e.channel], e));
                    svg.push_back(eye_svg(out.back()));
                }
            }
            return nb::make_tuple(eyes_json(out), eyes_summary(out), svg);
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def(
        "ibis_switching",
        [](const std::string& path, const std::string& model, bool rising, int corner) {
            const IbisFile f = read_ibis(path);
            const IbisModel* m = f.model(model);
            if (!m) throw std::runtime_error(path + " has no model " + model);
            std::vector<std::string> warnings = f.warnings;
            const IbisSwitching sw = switching_coefficients(*m, rising, corner, &warnings);
            return nb::make_tuple(sw.time_s, sw.ku, sw.kd, warnings, m->type, m->c_comp.at(corner), m->vpu(corner));
        },
        nb::arg("path"), nb::arg("model"), nb::arg("rising"), nb::arg("corner") = 0);
    m.def("prbs", [](const std::string& pattern, std::size_t bits) {
        Stimulus s;
        s.pattern = pattern;
        s.bits = bits;
        return s.sequence();
    });
    m.def("deembed", [](const std::string& yaml, const std::string& base_dir) {
        const DeembedResult r = deembed(yaml, base_dir);
        nb::list s;
        for (const auto& m : r.data.s) {
            nb::list row;
            for (const auto& c : m) row.append(nb::make_tuple(c.real(), c.imag()));
            s.append(row);
        }
        return nb::make_tuple(r.data.frequency_hz, s, r.data.reference_ohm, r.steps, r.output);
    });
    m.def("patch_antenna", [](double er, double h, double f, double w, double l, double feed) {
        return patch_json(patch_antenna(PatchSpec{er, h, f, w, l, feed}));
    });
    m.def("two_slot_directivity", &two_slot_directivity);
    m.def("parse_touchstone", [](const std::string& text, std::size_t ports) {
        const Touchstone t = parse_touchstone(text, ports);
        nb::list s;
        for (const auto& m : t.s) {
            nb::list row;
            for (const auto& c : m) row.append(nb::make_tuple(c.real(), c.imag()));
            s.append(row);
        }
        return nb::make_tuple(t.frequency_hz, s, t.reference_ohm);
    });
    m.def(
        "cross_section",
        [](const std::string& spec) {
            nb::gil_scoped_release release;
            return solve_cross_section_json(spec);
        },
        nb::arg("spec"));
    m.def(
        "export_openems",
        [](const Board& b, const std::string& yaml, const std::string& base_dir, const std::string& out_dir, const std::string& loop) {
            EmiProject p = parse_emi_project(yaml, base_dir);
            if (!p.emi || p.emi->loops.empty()) throw std::runtime_error("no emi loops");
            const LoopSpec* l = &p.emi->loops.front();
            if (!loop.empty()) {
                l = nullptr;
                for (const auto& x : p.emi->loops)
                    if (x.name == loop) l = &x;
                if (!l) throw std::runtime_error("no loop " + loop);
            }
            nlohmann::json model;
            {
                nb::gil_scoped_release release;
                model = out_dir.empty() ? openems_model(b, *l, *p.emi) : export_openems(b, *l, *p.emi, out_dir);
            }
            return std::make_pair(model.dump(), openems_summary(model));
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir"), nb::arg("out_dir"), nb::arg("loop") = "");
    m.def(
        "export_openems_channel",
        [](const Board& b, const std::string& yaml, const std::string& base_dir, const std::string& out_dir, const std::string& channel) {
            SiProject p = parse_si_project(yaml, base_dir);
            if (!p.si) throw std::runtime_error("no si section");
            const ChannelProject cp = parse_channel_project(yaml);
            const ChannelSpec* c = cp.channels.empty() ? nullptr : &cp.channels.front();
            if (!channel.empty()) {
                c = nullptr;
                for (const auto& x : cp.channels)
                    if (x.name == channel) c = &x;
            }
            if (!c) throw std::runtime_error("no channel " + channel);
            const OpenEmsChannelSettings s = parse_openems_channel(yaml);
            nlohmann::json model;
            {
                nb::gil_scoped_release release;
                model = out_dir.empty() ? openems_channel_model(b, *p.si, *c, s) : export_openems_channel(b, *p.si, *c, s, out_dir);
            }
            return std::make_pair(model.dump(), openems_channel_summary(model));
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir"), nb::arg("out_dir"), nb::arg("channel") = "");
    m.def("trapezoid_harmonics", [](double amplitude, double frequency, double duty, double rise_s, double fall_s, double f_max) {
        LoopWaveform w{amplitude, frequency, duty, rise_s, fall_s};
        return trapezoid_harmonics(w, f_max);
    });
    m.def(
        "loop_field_per_ampere",
        [](const std::vector<std::array<double, 3>>& positions, const std::vector<std::array<double, 3>>& moments, double frequency,
           const std::string& site, double distance, double board_height, double h_min, double h_max, bool far_only) {
            if (positions.size() != moments.size()) throw std::runtime_error("positions and moments differ in length");
            std::vector<CurrentElement> c(positions.size());
            std::array<double, 3> mo{0, 0, 0};
            for (std::size_t i = 0; i < c.size(); ++i) {
                c[i] = {positions[i], moments[i]};
                const auto& p = positions[i];
                const auto& d = moments[i];
                mo[0] += 0.5 * (p[1] * d[2] - p[2] * d[1]);
                mo[1] += 0.5 * (p[2] * d[0] - p[0] * d[2]);
                mo[2] += 0.5 * (p[0] * d[1] - p[1] * d[0]);
            }
            LoopSettings s;
            s.field_site = site;
            s.distance_m = distance;
            s.board_height_m = board_height;
            s.antenna_min_height_m = h_min;
            s.antenna_max_height_m = h_max;
            nb::gil_scoped_release release;
            return loop_field_per_ampere(c, mo, frequency, s, far_only);
        },
        nb::arg("positions"), nb::arg("moments"), nb::arg("frequency_hz"), nb::arg("site") = "free", nb::arg("distance_m") = 3.0,
        nb::arg("board_height_m") = 0.8, nb::arg("antenna_min_height_m") = 1.0, nb::arg("antenna_max_height_m") = 4.0,
        nb::arg("far_only") = false);

    m.def(
        "solve_dc",
        [](const Board& b, const std::string& settings_yaml, const std::string& base_dir) {
            const Project p = parse_project(settings_yaml, base_dir);
            nb::gil_scoped_release release;
            return solve_dc(b, p.dc);
        },
        nb::arg("board"), nb::arg("settings"), nb::arg("base_dir") = ".");
    m.def("dc_summary", &dc_summary);

    m.def(
        "open",
        [](const std::vector<std::string>& inputs, double arc_tolerance_um, std::optional<std::string> stackup,
           std::optional<std::array<double, 4>> crop_mm) {
            nb::gil_scoped_release release;
            return open_board(inputs, options(arc_tolerance_um, stackup, crop_mm));
        },
        nb::arg("inputs"), nb::arg("arc_tolerance_um") = 1.0, nb::arg("stackup") = nb::none(), nb::arg("crop_mm") = nb::none());

    m.def(
        "serve",
        [](const Board& b, const std::string& host, int port, bool open_browser, const DcResult* dc, const ThermalResult* th,
           const EmiResult* emi, const SiResult* si) {
            ServeOptions o;
            o.host = host;
            o.port = port;
            o.open_browser = open_browser;
            o.should_stop = []() {
                nb::gil_scoped_acquire acquire;
                return PyErr_CheckSignals() != 0;
            };
            {
                nb::gil_scoped_release release;
                serve(b, o, dc, th, emi, si);
            }
            if (PyErr_Occurred()) throw nb::python_error();
        },
        nb::arg("board"), nb::arg("host") = "127.0.0.1", nb::arg("port") = 8150, nb::arg("open_browser") = true,
        nb::arg("dc").none() = nb::none(), nb::arg("thermal").none() = nb::none(), nb::arg("emi").none() = nb::none(),
        nb::arg("si").none() = nb::none());

    m.attr("MODEL_VERSION") = kModelVersion;
    m.def("build_info", [] {
        const BuildInfo b = build_info();
        return nb::make_tuple(b.version, b.commit, b.compiler, b.platform);
    });
    m.def(
        "provenance",
        [](const std::vector<std::string>& inputs, const std::string& project, const Board* board) {
            return provenance_json(make_provenance({}, inputs, project, board)).dump();
        },
        nb::arg("inputs"), nb::arg("project"), nb::arg("board").none());
}
