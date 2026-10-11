// Rectangular microstrip patch by the transmission-line model (docs/solvers/rf.md), after
// Balanis, Antenna Theory, chapter 14, and Jackson and Alexopoulos (1991) for the bandwidth.
#include "kairoq/rf.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <stdexcept>

namespace kairoq {

namespace {

constexpr double kC = 299792458.0;
constexpr double kPiR = 3.14159265358979323846;

// Composite Simpson rule of f over [a, b] with n (even) intervals.
template <class F>
double simpson(F f, double a, double b, int n) {
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
    return s * h / 3;
}

// sin(a cos t) / cos t, with its limit a at cos t = 0.
double slot(double a, double t) {
    const double c = std::cos(t);
    return std::abs(c) < 1e-12 ? a : std::sin(a * c) / c;
}

}  // namespace

double two_slot_directivity(double w, double s) {
    // Ground plane z = 0, slots along y at x = +-s/2: U = sinc^2(k w/2 sin t sin p)
    // (1 - sin^2 t sin^2 p) cos^2(k s/2 sin t cos p) over the upper half-space, 1 at broadside.
    const double k = 2 * kPiR;
    auto u = [&](double t, double p) {
        const double st = std::sin(t);
        const double x = 0.5 * k * w * st * std::sin(p);
        const double sinc = std::abs(x) < 1e-12 ? 1.0 : std::sin(x) / x;
        const double a = std::cos(0.5 * k * s * st * std::cos(p));
        return sinc * sinc * (1 - st * st * std::sin(p) * std::sin(p)) * a * a;
    };
    const double total = simpson([&](double t) { return std::sin(t) * simpson([&](double p) { return u(t, p); }, 0, 2 * kPiR, 360); }, 0,
                                 0.5 * kPiR, 180);
    return 4 * kPiR / total;
}

PatchResult patch_antenna(const PatchSpec& sp) {
    if (!(sp.epsilon_r >= 1) || !(sp.height_m > 0)) throw std::runtime_error("patch: epsilon_r at least 1 and height_m positive");
    if (!(sp.frequency_hz > 0) && !(sp.width_m > 0 && sp.length_m > 0))
        throw std::runtime_error("patch: a design frequency, or the width and length");
    PatchResult r;
    const double er = sp.epsilon_r, h = sp.height_m;
    r.width_m = sp.width_m > 0 ? sp.width_m : kC / (2 * sp.frequency_hz) * std::sqrt(2 / (er + 1));
    const double w = r.width_m;
    r.epsilon_eff = (er + 1) / 2 + (er - 1) / 2 / std::sqrt(1 + 12 * h / w);
    const double ee = r.epsilon_eff;
    r.delta_l_m = 0.412 * h * (ee + 0.3) * (w / h + 0.264) / ((ee - 0.258) * (w / h + 0.8));
    r.length_m = sp.length_m > 0 ? sp.length_m : kC / (2 * sp.frequency_hz * std::sqrt(ee)) - 2 * r.delta_l_m;
    const double l = r.length_m;
    if (!(l > 0)) throw std::runtime_error("patch: no positive length");
    r.resonant_hz = kC / (2 * (l + 2 * r.delta_l_m) * std::sqrt(ee));
    const double k0 = 2 * kPiR * r.resonant_hz / kC, lambda = kC / r.resonant_hz;
    // Slot conductance and mutual conductance of the two radiating edges.
    const double a = 0.5 * k0 * w;
    const double i1 = simpson([&](double t) { const double v = slot(a, t); return v * v * std::pow(std::sin(t), 3); }, 0, kPiR, 4000);
    const double i12 = simpson([&](double t) {
        const double v = slot(a, t);
        return v * v * std::cyl_bessel_j(0.0, k0 * l * std::sin(t)) * std::pow(std::sin(t), 3);
    }, 0, kPiR, 4000);
    r.g1_s = i1 / (120 * kPiR * kPiR);
    r.g12_s = i12 / (120 * kPiR * kPiR);
    r.edge_resistance_ohm = 1 / (2 * (r.g1_s + r.g12_s));
    if (sp.feed_ohm > 0 && sp.feed_ohm <= r.edge_resistance_ohm)
        r.inset_m = l / kPiR * std::acos(std::sqrt(sp.feed_ohm / r.edge_resistance_ohm));
    else if (sp.feed_ohm > 0)
        r.warnings.push_back("the feed impedance exceeds the edge resistance: no inset matches it");
    r.directivity_dbi = 10 * std::log10(two_slot_directivity(w / lambda, (l + 2 * r.delta_l_m) / lambda));
    r.bandwidth = 3.771 * (er - 1) / (er * er) * (w / l) * (h / lambda);
    if (h > 0.05 * lambda / std::sqrt(er))
        r.warnings.push_back("the substrate is thicker than about 0.05 wavelengths in it: the model's bandwidth and resistance lose accuracy");
    if (er == 1.0) r.warnings.push_back("air substrate: the bandwidth formula gives zero");
    return r;
}

std::string patch_json(const PatchResult& r) {
    nlohmann::json j = {{"width_m", r.width_m},
                        {"length_m", r.length_m},
                        {"epsilon_eff", r.epsilon_eff},
                        {"delta_l_m", r.delta_l_m},
                        {"resonant_hz", r.resonant_hz},
                        {"g1_s", r.g1_s},
                        {"g12_s", r.g12_s},
                        {"edge_resistance_ohm", r.edge_resistance_ohm},
                        {"inset_m", r.inset_m},
                        {"directivity_dbi", r.directivity_dbi},
                        {"bandwidth", r.bandwidth},
                        {"warnings", r.warnings}};
    return j.dump();
}

}  // namespace kairoq
