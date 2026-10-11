// Closed-form estimates of RF structures: the rectangular microstrip patch antenna by the
// transmission-line (two-slot) model. Method: docs/solvers/rf.md.
#pragma once

#include <string>
#include <vector>

namespace kairoq {

struct PatchSpec {
    double epsilon_r = 4.3;
    double height_m = 1.6e-3;     // substrate
    double frequency_hz = 0.0;    // design frequency: the width and length where not given
    double width_m = 0.0;         // radiating edges
    double length_m = 0.0;        // resonant dimension
    double feed_ohm = 50.0;       // for the inset of the feed
};

struct PatchResult {
    double width_m = 0, length_m = 0;
    double epsilon_eff = 0, delta_l_m = 0;  // fringing extension of each edge
    double resonant_hz = 0;
    double g1_s = 0, g12_s = 0;             // slot conductance, mutual conductance
    double edge_resistance_ohm = 0;         // input resistance at a radiating edge
    double inset_m = 0;                     // feed point from the edge for feed_ohm
    double directivity_dbi = 0;             // two slots over an infinite ground plane
    double bandwidth = 0;                   // fraction, VSWR below 2 (thin substrates)
    std::vector<std::string> warnings;
};

PatchResult patch_antenna(const PatchSpec& spec);

// Directivity (linear) of two parallel slots of length w at spacing s over an infinite
// ground plane, both in wavelengths.
double two_slot_directivity(double w_wavelengths, double s_wavelengths);

std::string patch_json(const PatchResult& result);

}  // namespace kairoq
