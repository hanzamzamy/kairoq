// Export of a loop to openEMS, an open-source FDTD solver, for a full-wave check of the
// loop estimates: the copper, plated holes, dielectric and elements of the loop's model,
// with the first element as a lumped port. Description: docs/solvers/openems.md.
#pragma once

#include "kairoq/emi.hpp"
#include "kairoq/network.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace kairoq {

// The model (model.json): geometry in metres, simulation settings and, if asked for, the
// loop's estimate.
nlohmann::json openems_model(const Board& board, const LoopSpec& loop, const EmiSettings& settings);

// Write model.json and simulate.py, the script that builds and runs the simulation, into
// the directory (created if missing). Returns the model.
nlohmann::json export_openems(const Board& board, const LoopSpec& loop, const EmiSettings& settings, const std::string& dir);

std::string openems_summary(const nlohmann::json& model);

// Export of a channel's structure (docs/solvers/rf.md): all copper, plated holes and
// dielectric within a margin around the copper of the ports' nets, with a lumped port from
// each port's pad to the copper beneath it. The script excites each port in turn and
// writes the S-parameters as a Touchstone file; with far-field frequencies, the radiation
// of port 1 (directivity, efficiency, gain).
struct OpenEmsChannelSettings {
    double f_max_hz = 20.0e9;          // the excitation covers 0 to f_max
    std::size_t points = 201;          // frequencies of the result, linear from f_max / points
    double margin_nm = 3.0e6;          // around the ports' nets' copper
    double min_cell_nm = 0.0;          // 0: a quarter of the thinnest dielectric, at most 0.1 mm
    double board_cell_nm = 0.0;        // 0: a twentieth of the wavelength in the densest dielectric at f_max, at most 4 min_cell
    double air_margin_nm = 0.0;        // 0: a quarter wavelength at f_max
    std::vector<double> far_field_hz;  // radiation of port 1 at these frequencies
    double end_criterion_db = -40.0;
    std::size_t max_timesteps = 500000;
};

nlohmann::json openems_channel_model(const Board& board, const SiSettings& si, const ChannelSpec& channel,
                                     const OpenEmsChannelSettings& settings);
nlohmann::json export_openems_channel(const Board& board, const SiSettings& si, const ChannelSpec& channel,
                                      const OpenEmsChannelSettings& settings, const std::string& dir);
std::string openems_channel_summary(const nlohmann::json& model);

// Project file: `si.openems`.
OpenEmsChannelSettings parse_openems_channel(std::string_view yaml);

}  // namespace kairoq
