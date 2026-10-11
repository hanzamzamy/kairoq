// Impedance of a power distribution network from DC to a few GHz: the pairs of power and
// ground planes as cavities (their resonances), the decoupling capacitors and the
// regulator as branches at the vias that join them to the planes, seen from the loads.
// Method: docs/solvers/pdn.md.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/cross_section.hpp"

#include <complex>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

struct PdnCapacitor {
    std::string component;
    double c_f = 0, esr_ohm = 0, esl_h = 0;
    std::optional<double> mount_h;  // mounting inductance; estimated from the vias when not given
};

struct PdnSettings {
    std::string power_net, ground_net;
    double f_min_hz = 1e5, f_max_hz = 1e9;
    std::size_t points = 300;  // logarithmic
    std::optional<double> target_ohm;
    std::vector<std::string> loads;            // components; Z(f) is reported at each
    std::vector<PdnCapacitor> capacitors;
    // Every other part of two pads, one on each net, whose reference starts with this
    // prefix, with these values.
    std::optional<PdnCapacitor> other_capacitors;
    std::string capacitor_prefix = "C";
    struct Regulator {
        std::string component;
        double r_ohm = 1e-3, l_h = 1e-9;
    };
    std::optional<Regulator> regulator;
    double plane_fraction = 0.1;  // a net's copper on a layer this share of the board is a plane
    double via_reach_nm = 2e6;     // a pad joins the planes through vias of its net this close
    DielectricModel dielectric{4.3, 0.02, 1.0e9, true};  // where the stackup gives none
};

struct PdnPair {
    int power_layer = 0, ground_layer = 0;
    double height_m = 0, area_m2 = 0;
    std::size_t modes = 0, nodes = 0, lanczos_steps = 0;
    double seconds = 0;  // modes and static part
};

struct PdnBranch {
    std::string component;
    std::string kind;            // load, capacitor, regulator
    double x_mm = 0, y_mm = 0;   // the via (or pin) joining it to the planes
    double mount_h = 0;          // mounting inductance used (branches other than loads)
};

struct PdnResult {
    std::vector<PdnPair> pairs;
    std::vector<PdnBranch> branches;
    std::vector<std::string> loads;
    std::vector<double> frequency_hz;
    std::vector<std::vector<std::complex<double>>> z_ohm;  // per load, per frequency
    double seconds = 0;
    std::vector<std::string> warnings;
};

PdnResult analyse_pdn(const Board& board, const PdnSettings& settings);

// Project file: the `pdn` section.
PdnSettings parse_pdn(std::string_view yaml);

std::string pdn_summary(const PdnResult& result, const PdnSettings& settings);
std::string pdn_json(const PdnResult& result);

}  // namespace kairoq
