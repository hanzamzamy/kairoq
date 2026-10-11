// Eye diagrams of a channel with linear drivers and receivers: the response to a single
// bit, the eye of a bit pattern by superposition, and the statistical eye with noise,
// jitter and equalisation (FFE, CTLE, DFE) at a bit error ratio. Method:
// docs/solvers/eye.md.
#pragma once

#include "kairoq/network.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

struct Ctle {
    double dc_gain_db = 0.0;
    double f_zero_hz = 0.0;   // 0: no CTLE
    double f_pole1_hz = 0.0, f_pole2_hz = 0.0;
};

struct EyeSpec {
    std::string name, channel;
    std::size_t tx_p = 0, tx_n = SIZE_MAX;  // driver ports (0-based); tx_n: differential
    std::size_t rx_p = 0, rx_n = SIZE_MAX;  // receiver ports
    double bit_rate = 5e9;
    double swing_v = 1.0;      // open-circuit peak-to-peak of the source (differential, or single-ended)
    double r_source_ohm = 50;  // per leg
    double rise_s = 30e-12;    // 10-90 % of the source's edges
    double r_term_ohm = 50, c_term_f = 0;  // receiver, per leg to the reference
    std::string pattern = "prbs7";
    std::size_t bits = 4096;
    std::size_t samples_per_ui = 32;
    // Statistical eye.
    double noise_rms_v = 0, rj_rms_s = 0, dj_pp_s = 0;
    double ber = 1e-12;
    // Equalisation.
    std::vector<double> ffe;   // transmitter taps in order, with ffe_main the cursor tap
    std::size_t ffe_main = 0;
    Ctle ctle;
    std::size_t dfe_taps = 0;
    double dfe_limit = 1.0;    // largest tap, as a fraction of the main cursor
};

struct EyeResult {
    std::string name, channel;
    double ui = 0;
    // The response to one bit (source high for one UI from 0, relative to low): time, V.
    std::vector<double> pulse_t, pulse_v;
    std::size_t main_index = 0;           // sample of the pulse's peak
    std::vector<double> cursors;          // at the best phase, from main - pre to the end
    std::size_t pre_cursors = 0;
    std::vector<double> dfe;              // tap values, V
    // Eye of the pattern by superposition: density over two UI.
    std::vector<double> eye_t;            // s, 2 UI
    std::vector<double> eye_v;            // V, bin centres
    std::vector<std::vector<double>> density;  // [time][voltage], fraction of the bits
    // Statistical eye at the BER: per phase over two UI centred on the response's peak,
    // the inner edges of the eye.
    std::vector<double> phase_t, upper_v, lower_v;  // lower of the 1 level, upper of the 0 level
    double eye_height_v = 0, eye_width_ui = 0, best_phase_s = 0;
    double peak_distortion_v = 0;         // worst-case opening at the best phase, without noise and jitter
    double threshold_v = 0;
    double seconds = 0;
    std::vector<std::string> warnings;
};

EyeResult eye(const ChannelResult& channel, const EyeSpec& spec);

// Project file: `si.eyes`.
std::vector<EyeSpec> parse_eyes(std::string_view yaml);
std::string eyes_summary(const std::vector<EyeResult>& results);
std::string eyes_json(const std::vector<EyeResult>& results);
std::string eye_svg(const EyeResult& result);

}  // namespace kairoq
