// Networks of a board's nets: S-parameters between pads, with the tracks as
// transmission lines (single, or coupled where a pair runs together), vias, series parts
// and loads; mixed-mode parameters of pairs; time-domain reflectometry; Touchstone files.
// Method: docs/solvers/networks.md.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/si.hpp"

#include <complex>
#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

using cplx = std::complex<double>;

// S-parameters of a part, read from a Touchstone file.
struct Touchstone {
    std::vector<double> frequency_hz;
    std::vector<std::vector<cplx>> s;  // per frequency, N x N row major
    std::size_t ports = 0;
    double reference_ohm = 50.0;
    // S at a frequency, interpolated linearly in magnitude and unwrapped phase; the end
    // values outside the range.
    std::vector<cplx> at(double frequency_hz) const;
};
Touchstone read_touchstone(const std::string& path);
Touchstone parse_touchstone(const std::string& text, std::size_t ports);

// A part between two pads of the channel (a series capacitor, resistor, ferrite...).
struct ChannelPart {
    std::string from_pad, to_pad;  // REF.PIN
    double r_ohm = 0.0, l_h = 0.0, c_f = 0.0;  // in series; c_f 0: no capacitor
    std::string touchstone;        // a two-port file instead (path), port 1 at from_pad
};
// A load from a pad to the reference.
struct ChannelLoad {
    std::string pad;
    double r_ohm = 0.0;  // 0: no resistor
    double c_f = 0.0;    // in parallel
};

struct ChannelSpec {
    std::string name;
    std::vector<std::string> ports;           // pads (REF.PIN), one single-ended port each
    std::vector<std::pair<std::size_t, std::size_t>> differential;  // port pairs (0-based) of the mixed-mode result
    double impedance_ohm = 50.0;              // reference of each single-ended port
    std::vector<ChannelPart> parts;
    std::vector<ChannelLoad> loads;
    std::vector<std::size_t> drive;           // ports (0-based) to step for waveforms
};

struct NetworkSettings {
    double f_max_hz = 20.0e9;
    std::size_t points = 400;     // linear from f_max / points to f_max
    bool vias = true;             // vias in the channel (false: ideal connections)
    // "physical": barrel capacitances to the reference planes and the planes' impedance;
    // "coaxial": coaxial lines between the layers, of the antipad's diameter.
    std::string via_model = "physical";
    double rise_time_s = 35e-12;  // of the TDR step (10 to 90 %), at least 0.9 / f_max
    double tdr_duration_s = 0.0;  // 0: twice the summed delay of the channel's nets, at least 1 ns
};

struct TdrTrace {
    std::string port;           // port name, or "d1" for a differential port
    std::vector<double> time_s;
    std::vector<double> impedance_ohm;
};

// Waveforms at the ports for a step into one: fractions of the incident step.
struct StepResponse {
    std::string driven;
    std::vector<std::string> ports;
    std::vector<double> time_s;
    std::vector<std::vector<double>> values;  // per port
};

struct ChannelResult {
    std::string name;
    std::vector<std::string> ports;
    double impedance_ohm = 50.0;  // reference of the single-ended ports
    std::vector<std::string> nets;
    std::vector<double> frequency_hz;
    std::vector<std::vector<cplx>> s;      // per frequency, P x P
    std::vector<std::vector<cplx>> mixed;  // per frequency, differential then common modes (2K x 2K)
    std::vector<TdrTrace> tdr;
    std::vector<StepResponse> steps;
    std::size_t nodes = 0, sections = 0, coupled_sections = 0, vias = 0;
    double seconds = 0;
    std::vector<std::string> warnings;
};

// Project file: `si.channels` and `si.network`.
struct ChannelProject {
    std::vector<ChannelSpec> channels;
    NetworkSettings network;
};
ChannelProject parse_channel_project(std::string_view yaml);

std::vector<ChannelResult> analyse_channels(const Board& board, const SiSettings& si, const std::vector<ChannelSpec>& channels,
                                            const NetworkSettings& settings, const std::string& base_dir = ".");

// A Touchstone file (version 1, real and imaginary parts, Hz) of the single-ended or the
// mixed-mode parameters.
std::string touchstone(const ChannelResult& r, bool mixed);
std::string channels_summary(const std::vector<ChannelResult>& results);
std::string channels_json(const std::vector<ChannelResult>& results);

// Two-port ABCD matrix of a uniform line (characteristic impedance zc, propagation gamma,
// length in metres), and the S-parameters of a two-port from ABCD with reference z0.
std::array<cplx, 4> line_abcd(cplx zc, cplx gamma, double length_m);
std::array<cplx, 4> s_from_abcd(const std::array<cplx, 4>& abcd, double z0);

}  // namespace kairoq
