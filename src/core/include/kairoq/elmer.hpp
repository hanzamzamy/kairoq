// Export of a DC or a thermal problem to Elmer, an open-source finite element package, as a
// 3D reference for the 2.5D solvers. DC: the copper of the ports' nets as volumes of their
// thickness and the plated holes as columns, on one triangulation of the board shared by
// all layers; Elmer computes the conductances between the equipotential ports, and the
// ports' voltages and currents follow. Thermal: the whole board, copper and dielectric,
// through its thickness on such a triangulation. Description: docs/solvers/elmer.md.
#pragma once

#include "kairoq/dc.hpp"
#include "kairoq/thermal.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

struct ElmerSettings {
    double max_edge_nm = 0.5e6;   // largest triangle edge of the shared triangulation
    double min_size_nm = 0.05e6;  // at copper edges, holes and ports
    double grading = 0.3;         // growth of the size away from them
    int copper_divisions = 1;     // elements through the thickness of a copper layer
    std::string solver = "ElmerSolver";
};

struct ElmerExport {
    std::size_t nodes = 0, elements = 0, boundary_elements = 0, triangles = 0;
    std::vector<std::string> ports;
    std::vector<std::string> warnings;
};

// Writes mesh/ (Elmer's mesh files) and port_<k>.sif, one case per port at 1 V with the
// others at 0 V, into the directory.
ElmerExport export_elmer_dc(const Board& board, const DcSettings& dc, const ElmerSettings& settings, const std::string& dir);

struct ElmerDcResult {
    std::vector<DcPortResult> ports;
    std::vector<std::vector<double>> conductance_s;  // between the ports, I = G V
    double seconds = 0;
    std::vector<std::string> warnings;
};

// Runs the solver on each case of an exported directory and solves the ports' conditions.
ElmerDcResult run_elmer_dc(const DcSettings& dc, const std::string& dir, const std::string& solver);

// Writes mesh/, thermal.sif and probes.txt (the nodes of each heat source and copper layer)
// for the steady temperature with the thermal settings' sources and surfaces. Radiation is
// left out; heat sinks, fixed temperatures, junction resistances and Joule heat are not
// exported.
ElmerExport export_elmer_thermal(const Board& board, const ThermalSettings& thermal, const ElmerSettings& settings,
                                 const std::string& dir);

struct ElmerThermalResult {
    struct Source {
        std::string name;
        double mean_c = 0.0, max_c = 0.0;  // over the source's copper
    };
    std::vector<Source> sources;
    std::vector<double> layer_max_c;  // per copper layer, from index 1
    double seconds = 0;
};

// Runs the solver on an exported thermal case and reads the temperatures of the probes.
ElmerThermalResult run_elmer_thermal(const std::string& dir, const std::string& solver);

// Project file: `dc.elmer` or `thermal.elmer`.
ElmerSettings parse_elmer_settings(std::string_view yaml, const std::string& section = "dc");

std::string elmer_dc_summary(const ElmerExport& exported, const ElmerDcResult* result, const DcResult* native);
std::string elmer_thermal_summary(const ElmerExport& exported, const ElmerThermalResult* result, const ThermalResult* native);

}  // namespace kairoq
