// Time-domain response of a channel with its drivers and receivers: linear sources and
// loads, or IBIS buffers, at the channel's ports; the waveforms at the ports over time.
// Method: docs/solvers/timedomain.md.
#pragma once

#include "kairoq/network.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

// The bits a driver sends: "edge" (low, then high from the delay on), "pulse" (one bit
// high), "prbsN" (N = 7, 9, 15, 23, 31; a pseudo-random sequence), or a string of 0 and 1.
struct Stimulus {
    std::string pattern = "edge";
    double bit_rate = 1e9;    // bit/s
    std::size_t bits = 64;    // of a PRBS
    double delay_s = 100e-12; // before the first bit
    std::vector<int> sequence() const;  // the bits, in order
};

struct DriverSpec {
    std::size_t port = 0;      // channel port, 0-based
    Stimulus stimulus;
    // Linear: a source of the low and high levels behind a resistance; transitions are
    // linear over 1.25 times the 10-90 % time.
    double v_low = 0, v_high = 1, r_ohm = 50, rise_s = 50e-12, fall_s = -1;
    // IBIS: a model of a file, at a corner (0 typ, 1 min, 2 max); its package by default.
    std::string ibis_file, ibis_model;
    int corner = 0;
    bool package = true;
};

struct ReceiverSpec {
    std::size_t port = 0;
    double r_ohm = 0, v_term = 0, c_f = 0;  // r_ohm 0: no resistor
    std::string ibis_file, ibis_model;      // instead: an IBIS input model
    int corner = 0;
    bool package = true;
};

struct SimulationSpec {
    std::string name, channel;
    std::vector<DriverSpec> drivers;
    std::vector<ReceiverSpec> receivers;  // other ports are open
    double duration_s = 0;                // 0: the stimulus and 2 ns, at most half the sweep's period
    std::size_t oversample = 4;           // time step 1 / (2 oversample f_max)
};

struct PortMetrics {
    std::string port;
    double v_min = 0, v_max = 0, v_final = 0;
    double overshoot = 0, undershoot = 0;  // beyond the final and initial levels, V; single edges only
};

struct SimulationResult {
    std::string name, channel;
    std::vector<std::string> ports;
    std::vector<double> time_s;
    std::vector<std::vector<double>> voltage;  // per port
    std::vector<PortMetrics> metrics;
    double dt = 0;
    double seconds = 0;
    std::vector<std::string> warnings;
};

// Project file: `si.simulations`, each of a channel of `si.channels`.
std::vector<SimulationSpec> parse_simulations(std::string_view yaml);

SimulationResult simulate(const ChannelResult& channel, const SimulationSpec& spec, const std::string& base_dir = ".");

std::string simulations_summary(const std::vector<SimulationResult>& results);
std::string simulations_json(const std::vector<SimulationResult>& results);

}  // namespace kairoq
