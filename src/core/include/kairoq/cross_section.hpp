// Transmission-line cross-sections: the per-unit-length resistance, inductance,
// conductance and capacitance of conductors in a stack of dielectrics, by
// two-dimensional finite elements in the quasi-TEM approximation. Method:
// docs/solvers/si.md.
#pragma once

#include "kairoq/geometry.hpp"

#include <complex>
#include <string>
#include <vector>

namespace kairoq {

// Relative permittivity over frequency, as epsilon' - j epsilon''. Wideband: the
// Djordjevic-Sarkar model (a causal, nearly constant loss tangent) through the values at
// frequency_hz; otherwise constant.
struct DielectricModel {
    double epsilon_r = 1.0;
    double loss_tangent = 0.0;
    double frequency_hz = 1.0e9;
    bool wideband = true;
    double omega_low = 1.0e4, omega_high = 1.0e12;  // rad/s, the model's band
    std::complex<double> at(double frequency_hz) const;
};

// Conductor material: conductivity and surface roughness.
struct ConductorModel {
    double conductivity_s_per_m = 5.8e7;
    enum class Roughness { None, Hammerstad, Huray };
    Roughness roughness = Roughness::None;
    double rms_m = 0.0;           // Hammerstad: RMS height
    double huray_radius_m = 0.0;  // Huray: sphere radius
    double huray_ratio = 0.0;     // Huray: sphere surface per flat area
    // Increase of the high-frequency surface impedance over a smooth surface.
    double roughness_factor(double frequency_hz) const;
};

// The cut: x across, y up (nm).
struct CsConductor {
    std::string name;
    PolygonSet shape;
    bool reference = false;  // part of the return; all reference conductors are one
    ConductorModel metal;
};
struct CsDielectric {
    PolygonSet shape;  // later entries replace earlier ones where they overlap; elsewhere vacuum
    DielectricModel material;
};

struct CsSettings {
    double target_error = 0.02;  // estimated relative energy error; capacitance error about its square
    std::size_t max_triangles = 300000;
    int max_passes = 14;
    double margin = 1.0;  // grounded box beyond the geometry, in units of its larger extent
    std::vector<double> frequencies_hz{1.0e9};
    // Weight of each signal conductor's solution in the error estimate (empty: all 1);
    // conductors whose own quantities are not wanted can be solved less accurately.
    std::vector<double> weights;
    // Solve the field again at frequencies spread over the range (two per decade) instead
    // of weighting the energies of one solution.
    bool exact_dispersion = false;
    // Body of revolution: x is the radius from an axis at x = 0, which bounds the domain
    // without a condition; capacitances are those of the whole body (F, not F/m). Only the
    // capacitances are meaningful.
    bool axisymmetric = false;
};

struct CsFrequency {
    double frequency_hz = 0.0;
    std::vector<double> r, l, g, c;  // N x N, row major: ohm/m, H/m, S/m, F/m
};

// The frequency-independent results of a solved cross-section, from which its
// per-unit-length parameters follow at any frequency: the external inductance, each
// conductor's surface integral of the high-frequency currents with its effective
// thickness and material, and the energy of the field in each dielectric region.
struct CsModel {
    std::size_t n = 0;
    std::vector<double> l_ext;  // N x N, H/m
    struct Body {
        std::vector<double> m;  // N x N, 1/m
        double teff_m = 0.0;
        ConductorModel metal;
    };
    std::vector<Body> bodies;
    struct Region {
        std::vector<double> energy;  // N x N, F/m per unit of relative permittivity
        DielectricModel material;
        bool vacuum = false;
    };
    std::vector<Region> regions;
    std::vector<double> c_vacuum;  // N x N, F/m
    // Microstrip dispersion: set for a cut whose signal conductors lie on an outer layer
    // over one reference plane. at() then moves the capacitances towards the substrate's
    // as the frequency rises and scales the impedance, as the Kirschning-Jansen formulas
    // do for a single microstrip of the self conductor's width and height.
    struct Microstrip {
        bool enabled = false;
        std::size_t self = 0;  // signal conductor whose width and height set the formulas
        double width_m = 0, height_m = 0;
        DielectricModel substrate;
    } microstrip;
    // R, L, G, C at a frequency. The dielectric energies are those of the field at the
    // solved frequency, weighted with each region's permittivity at f: exact at the solved
    // frequency, and in error by the square of the relative change of the permittivities'
    // ratios elsewhere.
    CsFrequency at(double frequency_hz) const;
};

struct CsResult {
    std::vector<std::string> conductors;  // the signal conductors, in order
    std::vector<double> c0, l_ext;        // vacuum capacitance and external inductance (N x N)
    std::vector<CsFrequency> frequencies;
    CsModel model;
    double estimated_error = 0.0;
    std::size_t triangles = 0;
    int passes = 0;
    double seconds = 0.0;
    std::vector<std::string> warnings;
};

CsResult solve_cross_section(const std::vector<CsConductor>& conductors, const std::vector<CsDielectric>& dielectrics,
                             const CsSettings& settings);

// Quantities of one line (conductor i, the others at rest) and of a pair (conductors i
// and j) at one frequency.
struct LineQuantities {
    std::complex<double> z0;  // ohm
    double epsilon_eff = 0.0;
    double delay_s_per_m = 0.0;
    double loss_db_per_m = 0.0;
    double conductor_loss_db_per_m = 0.0, dielectric_loss_db_per_m = 0.0;
};
LineQuantities line_quantities(const CsFrequency& f, std::size_t n, std::size_t i);

struct PairQuantities {
    double z_odd = 0, z_even = 0, z_differential = 0, z_common = 0;  // ohm, lossless at the frequency
    double epsilon_eff_odd = 0, epsilon_eff_even = 0;
};
PairQuantities pair_quantities(const CsFrequency& f, std::size_t n, std::size_t i, std::size_t j);

// A cross-section given as JSON (docs/solvers/si.md), solved; the result as JSON.
std::string solve_cross_section_json(const std::string& spec);

// Upper end of the dispersion formulas' range, f h in GHz mm (Jansen and Kirschning).
constexpr double kMicrostripFnMax = 20.0;

// Kirschning and Jansen (1982): effective permittivity of a microstrip of width/height u
// at fn = f h (GHz mm) from its static value e0 on a substrate of permittivity er; Jansen
// and Kirschning (1983): the ratio of its characteristic impedance to the static one.
void kirschning_jansen(double u, double fn, double er, double e0, double& e_f, double& z_ratio);

}  // namespace kairoq
