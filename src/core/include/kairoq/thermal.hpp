// Thermal solver: layered conduction in copper and dielectric, convection and radiation
// from the faces, heat sinks through a thermal interface, heat from components and from the
// DC current (electrothermal coupling); steady state and transients. Method and validation:
// docs/solvers/thermal.md.
#pragma once

#include "kairoq/dc.hpp"

#include <array>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace kairoq {

struct SurfaceSpec {
    double convection_w_m2k = 10.0;  // natural convection is about 5 to 10
    double emissivity = 0.9;         // solder mask; 0 switches radiation off
};

// A heat source or a fixed temperature on an area selected as a DC port is.
struct HeatSpec {
    PortSpec area;  // name and selection; a component without a net selects all its pads
    std::optional<double> power_w;
    // Power over time for transients: (time s, power W) points, linear in between and
    // constant beyond the ends. A steady solution uses power_w, or else the power at t = 0.
    std::vector<std::array<double, 2>> power_schedule;
    std::optional<double> temperature_c;     // fixed temperature (a heat sink)
    std::optional<double> junction_k_per_w;  // junction-to-case resistance: the power
                                             // enters a junction node behind it
    double heat_capacity_j_per_k = 0.0;      // of the junction node, or else of the area
};

// A heat sink pressed onto a face of the board through a thermal interface material.
struct HeatsinkSpec {
    std::string name;
    struct Region {
        bool bottom = true;  // face
        PolygonSet area;     // empty: the whole face
    };
    std::vector<Region> regions;
    double tim_m2k_per_w = 0.0;  // resistance per area of the interface (thickness / k)
    // The sink: a fixed temperature, or a body with a resistance to ambient (as a heat sink
    // datasheet gives it) and a heat capacity.
    std::optional<double> temperature_c;
    std::optional<double> resistance_k_per_w;
    double heat_capacity_j_per_k = 0.0;
};

struct TransientSettings {
    double duration_s = 0.0;
    double step_s = 0.0;      // first step after the start and after each change of a schedule
    double max_step_s = 0.0;  // 0: duration / 50
    int steps_per_doubling = 8;
    bool start_steady = false;  // start from the steady solution at t = 0 instead of ambient
    std::vector<double> snapshot_times_s;  // fields kept at these times; empty: about 20
    double dc_update_k = 1.0;  // the DC solution is renewed when the temperature has
                               // changed by this much since the last one
    // Each step solves for the change of the temperature, to this tolerance relative to the
    // residual of the previous temperature (0: the whole state to the solver tolerance).
    double step_tolerance = 1e-4;
};

struct ThermalSettings {
    double ambient_c = 25.0;
    SurfaceSpec top, bottom;
    double copper_k = 385.0;          // W/(m K)
    double dielectric_k_xy = 0.8;     // FR-4, in-plane
    double dielectric_k_z = 0.3;      // FR-4, through the thickness
    std::map<std::string, std::array<double, 2>> layer_k;  // stackup layer name -> {k_xy, k_z}
    double mask_k = 0.25;
    // The space between the copper of the outer layers: air (nothing), or solder mask that
    // fills it to the copper's thickness.
    bool outer_gaps_mask = false;
    // Volumetric heat capacities, J/(m^3 K), for transients.
    double copper_heat_capacity = 3.45e6;
    double dielectric_heat_capacity = 2.2e6;
    double mask_heat_capacity = 1.8e6;
    std::map<std::string, double> layer_heat_capacity;  // stackup layer name -> J/(m^3 K)
    // Dielectrics are resolved by conduction sheets: the first slice at each face of a
    // dielectric has dielectric_first_slice, the next ones double up to dielectric_slice
    // (the layered model converges to a three-dimensional solution as they decrease;
    // docs/validation/thermal.md).
    double dielectric_first_slice_nm = 25000.0;
    double dielectric_slice_nm = 200000.0;
    double plating_thickness_nm = 18000.0;
    std::vector<HeatSpec> sources;
    std::vector<HeatsinkSpec> heatsinks;
    MeshSettings mesh = [] {
        MeshSettings m;
        m.grading = 0.0;
        m.boundary_tolerance_nm = 10000.0;
        m.max_edge_nm = 2.0e6;
        return m;
    }();
    double target_error = 0.1;  // relative, energy norm, as for the DC solver
    int max_passes = 6;
    std::size_t max_triangles = 5'000'000;
    double tolerance = 1e-10;
    int max_iterations = 2000;
    // Electrothermal coupling: the DC problem whose Joule heat is solved for, with the
    // copper conductivity at the computed temperature, until the temperature settles.
    std::optional<DcSettings> dc;
    double coupling_tolerance_k = 0.02;
    int max_coupling_iterations = 30;
    bool coupling_acceleration = true;  // Anderson acceleration of the coupling (depth 3)
    std::optional<TransientSettings> transient;
};

struct ThermalLayerField {
    int layer = 0;
    std::vector<double> x, y;  // nm
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<double> temperature_c;  // per vertex
    std::vector<std::uint8_t> copper;   // per triangle
};

struct ThermalSourceResult {
    std::string name;
    double power_w = 0.0;        // set power, or heat taken out by a fixed temperature
    bool temperature_set = false;
    double mean_c = 0.0, max_c = 0.0;  // over the area
    std::optional<double> junction_c;
};

struct ThermalHeatsinkResult {
    std::string name;
    double heat_w = 0.0;         // into the sink through the interface
    double temperature_c = 0.0;  // of the sink
    bool temperature_set = false;
};

// Time series of a transient: one entry per step, the start included.
struct ThermalTransientResult {
    std::vector<double> time_s;
    std::vector<double> max_c;  // over the copper layers
    std::vector<std::vector<double>> source_mean_c, source_max_c, source_junction_c;  // [source][step]
    std::vector<std::vector<double>> heatsink_c, heatsink_w;                          // [heatsink][step]
    std::vector<double> heat_in_w;      // sources and Joule heat
    std::vector<double> heat_out_w;     // faces, heat sinks and fixed temperatures
    std::vector<double> heat_stored_w;  // rate of change of the stored heat
    // Fields at selected times, on the mesh of `ThermalResult::layers`.
    std::vector<double> snapshot_time_s;
    std::vector<std::vector<std::vector<double>>> snapshot_c;  // [snapshot][layer][vertex]
    int steps = 0, matrix_setups = 0, dc_solutions = 0;
    long solver_iterations = 0;    // conjugate-gradient iterations over all steps
    double seconds_setup = 0.0;    // multigrid set-ups (included in seconds_solve)
    double seconds_steady = 0.0;   // steady solutions for the mesh and the start
};

struct ThermalResult {
    std::vector<ThermalLayerField> layers;
    std::vector<ThermalSourceResult> sources;
    std::vector<ThermalHeatsinkResult> heatsinks;
    double max_c = 0.0;
    int max_layer = 0;
    double max_x_nm = 0.0, max_y_nm = 0.0;
    double heat_sources_w = 0.0, heat_joule_w = 0.0;  // in
    double heat_top_w = 0.0, heat_bottom_w = 0.0, heat_fixed_w = 0.0, heat_sinks_w = 0.0;  // out
    std::size_t triangles = 0, unknowns = 0;
    int iterations = 0;
    double residual = 0.0;
    double estimated_error = 0.0;
    std::vector<DcResult::Pass> passes;
    int coupling_iterations = 0;
    double coupling_change_k = 0.0;  // last change of the temperature between iterations
    std::optional<DcResult> dc;
    std::optional<ThermalTransientResult> transient;  // the fields above are then the final state
    double seconds_mesh = 0.0, seconds_assemble = 0.0, seconds_solve = 0.0, seconds_dc = 0.0;
    std::vector<std::string> warnings;
};

ThermalResult solve_thermal(const Board& board, const ThermalSettings& settings);

// Component power from a bill of materials or a similar table (CSV): one source per
// component. Columns are found by name; references may be grouped in one cell ("R1, R2").
struct ComponentPowerTable {
    std::string path;
    std::string reference_column;      // empty: Reference, Designator, Ref or RefDes
    std::string power_column;          // empty: Power, Power_W, Power (W) or P
    std::string junction_column;       // optional: junction-to-case resistance, K/W
    std::string heat_capacity_column;  // optional: J/K
    bool ignore_unknown = false;       // references not on the board are skipped, not errors
};
std::vector<HeatSpec> read_component_power(const Board& board, const ComponentPowerTable& table,
                                           std::vector<std::string>& warnings);

// Project file section `thermal` (with `dc` for the coupling). Component power tables are
// read when the board is known (`resolve_thermal_project`).
struct ThermalProject {
    std::vector<std::string> board;
    std::optional<ThermalSettings> thermal;
    DcSettings dc;
    std::vector<ComponentPowerTable> component_power;
};
ThermalProject load_thermal_project(const std::string& path);
ThermalProject parse_thermal_project(std::string_view yaml, const std::string& base_dir);
// Settings with the component power tables read into sources.
ThermalSettings resolve_thermal_project(const Board& board, const ThermalProject& project, std::vector<std::string>& warnings);

std::string thermal_summary(const Board& board, const ThermalResult& result);

}  // namespace kairoq
