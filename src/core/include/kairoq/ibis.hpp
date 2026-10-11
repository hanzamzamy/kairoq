// IBIS buffer models (I/O Buffer Information Specification, versions 3 to 7): the
// keywords of a model that describe its behaviour at a pin, and its switching
// coefficients. Method: docs/solvers/timedomain.md.
#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace kairoq {

// A typical, minimum and maximum value; minimum and maximum default to typical.
struct IbisValue {
    double typ = 0, min = 0, max = 0;
    double at(int corner) const { return corner == 1 ? min : corner == 2 ? max : typ; }  // 0 typ, 1 min, 2 max
};

// A table of voltage against current (I-V) or against time (V-t), with its corners.
struct IbisTable {
    std::vector<double> x;                    // voltage (V) or time (s)
    std::array<std::vector<double>, 3> y;     // typ, min, max
    bool empty() const { return x.empty(); }
    double at(double v, int corner) const;    // linear interpolation, extended linearly beyond the ends
};

struct IbisWaveform {
    double r_fixture = 50, v_fixture = 0, c_fixture = 0, l_fixture = 0, r_dut = 0, l_dut = 0, c_dut = 0;
    IbisTable v;  // pin voltage against time
};

struct IbisModel {
    std::string name, type;   // Output, I/O, Input, 3-state, Open_drain, ...
    IbisValue c_comp;         // F
    IbisValue voltage_range;  // V
    std::optional<IbisValue> pullup_ref, pulldown_ref, power_clamp_ref, gnd_clamp_ref;
    IbisValue vinl, vinh;
    IbisTable pullup, pulldown, power_clamp, gnd_clamp;  // current into the pin, A
    IbisValue ramp_rise, ramp_fall;  // dV/dt of 20 to 80 % into r_load, V/s
    double ramp_r_load = 50;
    std::vector<IbisWaveform> rising, falling;
    // References of the tables at a corner.
    double vpu(int corner) const;  // pullup and power clamp
    double vpd(int corner) const;  // pulldown and ground clamp
    // Current into the pin with the output stage weighted by ku and kd, at pin voltage v.
    double current(double v, double ku, double kd, int corner) const;
    double clamp_current(double v, int corner) const;
};

struct IbisPackage {
    IbisValue r, l, c;  // R_pkg, L_pkg, C_pkg
};

struct IbisPin {
    std::string pin, signal, model;
    std::optional<double> r, l, c;  // per pin; else the package's
};

struct IbisFile {
    std::string version, component;
    IbisPackage package;
    std::vector<IbisPin> pins;
    std::vector<IbisModel> models;
    std::vector<std::string> warnings;
    const IbisModel* model(const std::string& name) const;
};

IbisFile parse_ibis(const std::string& text);
IbisFile read_ibis(const std::string& path);

// Switching coefficients of the output stage over time for a rising or falling edge:
// ku (pullup) and kd (pulldown) such that the model reproduces its two V-t waveforms
// into their fixtures; from the ramp when the waveforms are missing.
struct IbisSwitching {
    std::vector<double> time_s, ku, kd;
};
IbisSwitching switching_coefficients(const IbisModel& m, bool rising, int corner, std::vector<std::string>* warnings = nullptr);

}  // namespace kairoq
