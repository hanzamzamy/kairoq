// DC current solver: copper layers as conducting sheets (linear triangular finite
// elements), plated holes as barrel conductances between the layers they connect,
// ports as equipotential areas with a set voltage or an injected current.
// Method and validation: docs/solvers/dc.md.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/mesh2d.hpp"

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace kairoq {

// A port: an equipotential area of copper with either a voltage or a current.
struct PortSpec {
    std::string name;
    // Selection, one of:
    std::vector<std::string> pads;  // "REF.PIN": every pad of that pin, on every layer
    std::string component, net;     // every pad of the component on the net
    std::optional<Point64> at;      // a contact disc of `diameter_nm` on `layer`
    int layer = 0;
    double diameter_nm = 0.0;
    // Excitation, one of:
    std::optional<double> voltage_v;
    std::optional<double> current_a;  // flowing into the copper; a load is negative
};

struct DcSettings {
    double temperature_c = 20.0;
    double temperature_coefficient = 0.00393;  // of the copper resistivity, 1/K, at 20 C
    // Conductivity at 20 C for copper layers whose stackup gives none (IEC 60028
    // annealed copper); copper_conductivity overrides all layers when set.
    double default_conductivity_s_per_m = 5.8e7;
    std::optional<double> copper_conductivity_s_per_m;
    double plating_thickness_nm = 18000.0;
    double plating_conductivity_s_per_m = 5.8e7;  // at 20 C
    // Initial mesh: quality refinement only; the adaptive passes add the grading.
    MeshSettings mesh = [] {
        MeshSettings m;
        m.grading = 0.0;
        m.boundary_tolerance_nm = 5000.0;
        return m;
    }();
    // Adaptive refinement: passes stop when the estimated relative error of the
    // solution in the energy norm is below target_error. It bounds the relative error
    // of the total dissipation (and of a two-port resistance) by about its square;
    // quantities of a path that carries a small share of the power can be less accurate.
    double target_error = 0.1;
    int max_passes = 6;
    // Refinement stops before the mesh would exceed this many triangles (about 1 kB
    // of memory per triangle over the whole solution).
    std::size_t max_triangles = 5'000'000;
    double tolerance = 1e-10;  // relative residual of the linear solver
    int max_iterations = 2000;
    std::vector<PortSpec> ports;
    // Copper temperature (degrees C) by copper layer and position (nm); replaces
    // temperature_c where set (electrothermal coupling).
    std::function<double(int layer, double x_nm, double y_nm)> temperature_field;
    // Initial mesh vertices per copper island, from DcResult::mesh_seeds of an earlier
    // solution of the same board and ports; with max_passes = 1 the mesh is reused.
    std::vector<std::vector<std::array<double, 3>>> seeds;
};

// Solution on one copper layer: the mesh, nodal voltages and per-triangle results.
struct DcLayerField {
    int layer = 0;
    double thickness_nm = 0.0;
    double conductivity_s_per_m = 0.0;  // at the solution temperature
    std::vector<double> x, y;  // nm
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<double> voltage;   // per vertex, V
    std::vector<double> jx, jy;    // per triangle: sheet current density, A/m
    std::vector<double> power_w;   // per triangle: Joule heating, W
    std::vector<int> net;          // per triangle
};

struct DcPortResult {
    std::string name;
    double voltage_v = 0.0;
    double current_a = 0.0;  // into the copper
    bool voltage_set = false;
    std::size_t nodes = 0;  // mesh nodes tied together by the port
};

struct DcBarrelResult {
    std::size_t drill = 0;  // index into Board::drills
    int net = 0;
    int layer_a = 0, layer_b = 0;  // connected copper layers, layer_a < layer_b
    double resistance_ohm = 0.0;
    double current_a = 0.0;  // from layer_a to layer_b
};

struct DcHotSpot {
    int layer = 0;
    int net = 0;
    double x_nm = 0.0, y_nm = 0.0;
    double sheet_a_per_mm = 0.0;
    double density_a_per_mm2 = 0.0;
};

struct DcResult {
    std::vector<DcLayerField> layers;
    std::vector<DcPortResult> ports;
    std::vector<DcBarrelResult> barrels;
    std::vector<DcHotSpot> hot_spots;  // largest current densities, at most 20, apart by 1 mm
    double power_copper_w = 0.0;
    double power_barrels_w = 0.0;
    double power_ports_w = 0.0;  // sum of V I over the ports; equals the dissipation
    std::size_t unknowns = 0, triangles = 0, islands = 0;
    // Discretisation error estimate (gradient recovery) in the energy norm, relative,
    // and its history over the adaptive passes.
    double estimated_error = 0.0;
    struct Pass {
        std::size_t triangles;
        double estimated_error;
        double power_w;
    };
    std::vector<Pass> passes;
    int iterations = 0;
    double residual = 0.0;
    double seconds_mesh = 0.0, seconds_assemble = 0.0, seconds_solve = 0.0;
    std::vector<std::string> warnings;
    std::vector<std::vector<std::array<double, 3>>> mesh_seeds;  // see DcSettings::seeds
};

DcResult solve_dc(const Board& board, const DcSettings& settings);

// Project file (YAML): board inputs and DC settings. Paths are relative to the file.
struct Project {
    std::vector<std::string> board;
    DcSettings dc;
};
Project parse_project(std::string_view yaml, const std::string& base_dir);
Project load_project(const std::string& path);

// Summary text and JSON of a result.
std::string dc_summary(const Board& board, const DcResult& result);

}  // namespace kairoq
