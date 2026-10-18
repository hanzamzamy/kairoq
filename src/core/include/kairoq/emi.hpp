// First-order electromagnetic interference estimates: inductance and magnetic moment of
// switching loops (partial inductance, quasi-static), radiated field of their harmonics,
// return paths of traces over their reference planes, and the capacitance of switch nodes
// to a metal plate. Method and validation: docs/solvers/emi.md.
#pragma once

#include "kairoq/board.hpp"

#include <array>
#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace kairoq {

// One element of a loop: a component (or any part) carrying the loop current from one of
// its pads to another, as a straight conductor between the pad centres.
struct LoopElement {
    std::string from_pad, to_pad;      // "REF.PIN"
    double radius_nm = 200000.0;       // of the equivalent round conductor
    double extra_inductance_nh = 0.0;  // inside the package (equivalent series inductance)
};

// Current of a loop over one period: trapezoidal pulses.
struct LoopWaveform {
    double amplitude_a = 0.0;  // 0: no radiated estimate
    double frequency_hz = 0.0;
    double duty = 0.5;   // pulse width at half amplitude / period
    double rise_s = 0.0;
    double fall_s = -1.0;  // negative: as rise_s
};

// A part that offers the return current a path between two pads outside the loop's own
// chain: a stitching capacitor across a split or slot, or a strap. A straight conductor
// between the pad centres with a series capacitance (0: none), resistance and inductance.
struct LoopStitch {
    std::string from_pad, to_pad;
    double radius_nm = 200000.0;
    double capacitance_f = 0.0;
    double resistance_ohm = 0.0;
    double extra_inductance_nh = 0.0;
};

struct LoopSpec {
    std::string name;
    std::vector<LoopElement> elements;  // in the order the current flows; closed
    std::vector<LoopStitch> stitches;
    LoopWaveform waveform;
};

struct LoopSettings {
    double window_margin_nm = 10.0e6;     // copper this far around the loop is included
    double passive_min_area_mm2 = 10.0;   // smaller islands of other nets are left out
    double min_size_nm = 0.0;             // element size at the loop (0: from the pads)
    std::size_t max_triangles = 4000;
    int depth_layers = 2;                 // copper layers included beyond the loop's own (towards the inside)
    double via_radius_nm = 3.0e6;         // plated holes this close to the loop are modelled
    bool estimate_error = true;           // a second solution with about half the triangles
    double distance_m = 3.0;              // for the radiated field
    double max_frequency_hz = 1.0e9;
    // The window grows (its margin doubled, at most max_window_doublings times) until the
    // moment changes by less than moment_tolerance or the window stops growing; 0: one window.
    double moment_tolerance = 0.05;
    int max_window_doublings = 3;
    bool loop_basis = true;  // solve in a basis of loops (false: branches with the current law)
    // The dense matrix of partial inductances is kept up to this size (MiB); above it, in the
    // basis of loops, its rows are formed by blocks and multiplied out (about twice the
    // kernel evaluations, a fraction of the memory).
    double dense_matrix_mb = 2048.0;
    // Resistance and inductance of the loop at these frequencies (empty: the high-frequency
    // inductance only): the sheets' surface impedance and the holes' plating in the solve.
    std::vector<double> frequencies_hz;
    double plating_nm = 25.0e3;  // wall of plated holes
    // Radiated field: "free" (the field of the solved currents at distance_m in free space,
    // in the worst direction), "ground" (a measurement site: the board horizontal at
    // board_height_m over a reflecting ground, the antenna at distance_m, its height scanned,
    // both polarisations, the board turned), or "small_loop" (the far field of the moment).
    std::string field_site = "free";
    double board_height_m = 0.8;
    double antenna_min_height_m = 1.0, antenna_max_height_m = 4.0;
};

// A current element of the solved loop: its position (m) and current moment I dl (A m),
// per ampere of loop current.
struct CurrentElement {
    std::array<double, 3> position_m{0, 0, 0}, moment_a_m{0, 0, 0};
};

struct LoopImpedance {
    double frequency_hz = 0, resistance_ohm = 0, inductance_nh = 0;
    std::vector<std::complex<double>> stitch_current;  // per stitch, per ampere of loop current
};

struct LoopHarmonic {
    double frequency_hz = 0, current_a = 0, field_v_per_m = 0;
};

struct LoopSheetField {
    int layer = 0;
    std::vector<double> x, y;  // nm, vertices
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<std::array<double, 2>> current_a_per_m;  // sheet current at the centroid, per ampere of loop current
};

struct LoopResult {
    std::string name;
    double inductance_nh = 0.0;
    // At the previous refinement pass, with about half the triangles (0: not computed). The
    // inductance converges from above; the error is at least the difference, and at small
    // budgets several times it (docs/solvers/emi.md, "Accuracy and cost").
    double inductance_coarse_nh = 0.0;
    double element_size_mm = 0.0;  // of the first mesh at the loop, before refinement
    double element_inductance_nh = 0.0; // added inside the packages
    std::array<double, 3> moment_m2{0, 0, 0};  // magnetic moment per ampere (A m^2 / A)
    std::size_t triangles = 0, unknowns = 0, filaments = 0;
    std::vector<LoopHarmonic> harmonics;  // with a waveform
    std::vector<LoopImpedance> impedance;  // at LoopSettings::frequencies_hz
    std::vector<CurrentElement> currents;  // of the last solution
    // Current through each stitch per ampere of loop current in the high-frequency limit,
    // where a stitch is its inductance alone.
    std::vector<double> stitch_current;
    std::vector<std::pair<double, double>> field_per_ampere;  // frequency (Hz), largest field (V/m per A)
    std::string field_site;  // as LoopSettings::field_site, with the distance below
    double field_distance_m = 0.0;
    double max_field_v_per_m = 0.0, max_field_frequency_hz = 0.0;
    std::vector<LoopSheetField> sheets;
    double window_margin_mm = 0.0, window_area_mm2 = 0.0;
    double moment_change = 0.0;  // relative change of the moment at the last growth of the window
    double seconds = 0.0;
    std::vector<std::string> warnings;
};

// Return paths: traces of the checked nets over gaps in the copper of their reference
// layers (the nearest plane layers above and below).
struct ReturnPathSettings {
    std::vector<std::string> nets;            // empty: every net that is not a reference
    std::vector<std::string> reference_nets;  // empty: nets that fill a quarter of a layer
    double plane_fraction = 0.25;
    double min_area_mm2 = 0.002;              // smaller pieces are ignored
    double pad_clearance_nm = 300000.0;       // around the net's own pads and holes
};

struct ReturnPathFinding {
    std::string net;
    int layer = 0, reference_layer = 0;
    std::string kind;                     // "crossing", "split", "edge"
    std::vector<std::string> reference_nets;  // on either side
    double x_nm = 0, y_nm = 0;            // centroid of the piece over the gap
    double area_mm2 = 0, span_mm = 0;     // span: along the trace, area / trace width
    PolygonSet piece;
};

// Capacitance of switch nodes to a metal plate on a face of the board.
struct PlateSpec {
    std::string name;
    bool bottom = true;
    PolygonSet area;            // empty: the whole face
    double gap_nm = 0.0;        // insulation between face and plate (pad, film, air)
    double gap_epsilon_r = 1.0;
};
struct SwitchNodeSpec {
    std::string net;
    std::string plate;
    double dv_dt_v_per_ns = 0.0;  // for the common-mode current estimate
};
struct CouplingResult {
    std::string net, plate;
    double capacitance_pf = 0.0;          // to the plate, with fringing (method of moments)
    double parallel_plate_pf = 0.0;       // the unshielded area's parallel-plate capacitance
    double total_pf = 0.0;                // to the plate and all copper in the window
    // To the plate on the mesh used and on one with elements twice as large: the capacitance
    // is extrapolated from the two (0: no coarser mesh; the capacitance is then not extrapolated).
    double mesh_pf = 0.0, coarse_pf = 0.0;
    std::vector<std::pair<std::string, double>> nets_pf;  // to each net in the window
    std::size_t panels = 0;
    double window_margin_mm = 0.0;
    std::vector<std::string> warnings;
    std::vector<std::pair<int, double>> layers;  // layer, capacitance pF
    double unshielded_area_mm2 = 0.0, shielded_area_mm2 = 0.0;
    double common_mode_current_ma = 0.0;
    // Overlap with the copper of other nets on the neighbouring layers: net, layer, pF.
    std::vector<std::tuple<std::string, int, double>> neighbours;
};

// Export of a loop to openEMS (src/core/include/kairoq/openems.hpp).
struct OpenEmsSettings {
    double f_max_hz = 1.0e9;           // the excitation covers 0 to f_max
    double min_cell_nm = 100000.0;     // mesh cell at copper edges, sheets and elements
    double board_cell_nm = 1.0e6;      // largest cell within the window
    double air_margin_nm = 0.0;        // air around the window; 0: a quarter wavelength at f_max
    double port_resistance_ohm = 50.0;
    std::vector<double> far_field_hz;  // empty: f_max / 10, / 4 and / 2
    double end_criterion_db = -40.0;   // energy decay at which the simulation stops
    std::size_t max_timesteps = 1000000;
    bool estimate = true;              // include the loop's estimate for comparison
};

struct EmiSettings {
    std::vector<LoopSpec> loops;
    LoopSettings loop;
    std::optional<ReturnPathSettings> return_path;
    std::vector<PlateSpec> plates;
    std::vector<SwitchNodeSpec> switch_nodes;
    double dielectric_epsilon_r = 4.3;  // where the stackup gives none
    double mask_epsilon_r = 3.8;
    // Switch-node capacitance by the method of moments: the copper within this margin of
    // the node (0: three times its height over the plate, 2 to 30 mm), in at most this many
    // triangles.
    struct Coupling {
        double window_margin_nm = 0.0;
        std::size_t max_panels = 4000;
        double min_size_nm = 0.0;  // at the node's edges (0: a third of its distance to the nearest other conductor)
    } coupling;
    OpenEmsSettings openems;
};

struct EmiResult {
    std::vector<LoopResult> loops;
    std::vector<ReturnPathFinding> return_path;
    std::vector<CouplingResult> coupling;
    std::vector<std::string> warnings;
};

LoopResult solve_loop(const Board& board, const LoopSpec& loop, const LoopSettings& settings);
std::vector<ReturnPathFinding> check_return_paths(const Board& board, const ReturnPathSettings& settings,
                                                  std::vector<std::string>& warnings);
CouplingResult switch_node_coupling(const Board& board, const SwitchNodeSpec& node, const PlateSpec& plate,
                                    const EmiSettings& settings);
EmiResult analyse_emi(const Board& board, const EmiSettings& settings);

// Harmonics of a trapezoidal current, one-sided amplitudes up to f_max (exact Fourier
// coefficients of the piecewise-linear waveform).
std::vector<std::pair<double, double>> trapezoid_harmonics(const LoopWaveform& w, double f_max_hz);
// Largest far field (V/m) of a small loop of moment m (A m^2) at distance r, frequency f.
double small_loop_field(double moment_a_m2, double frequency_hz, double distance_m);
// Largest field (V/m) per ampere of loop current at frequency f, at the settings' site,
// from the current elements (or the moment, for "small_loop"). far_only: the far-field term
// alone, scaled to distance_m, as a near-to-far-field transformation gives it (free site).
double loop_field_per_ampere(const std::vector<CurrentElement>& currents, const std::array<double, 3>& moment_m2, double frequency_hz,
                             const LoopSettings& settings, bool far_only = false);

// Project file section `emi`.
struct EmiProject {
    std::vector<std::string> board;
    std::optional<EmiSettings> emi;
};
EmiProject load_emi_project(const std::string& path);
EmiProject parse_emi_project(std::string_view yaml, const std::string& base_dir);
std::string emi_summary(const Board& board, const EmiResult& result);

}  // namespace kairoq
