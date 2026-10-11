// Signal integrity of a board's tracks: the cross-section of each track, cut across it at
// intervals, solved for its per-unit-length parameters (kairoq/cross_section.hpp); the
// impedance profile, delay and loss of each net, the differential impedance of pairs,
// and findings where they leave their target or lose their reference. Method:
// docs/solvers/si.md.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/cross_section.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace kairoq {

struct SiTarget {
    std::vector<std::string> nets;  // patterns (* and ?); empty: every analysed net
    double ohm = 0.0;               // single-ended, or differential for pairs
    bool differential = false;
    double tolerance = 0.1;         // relative
};

struct SiPair {
    std::string p, n;
};

struct SiSettings {
    std::vector<std::string> nets;            // patterns; empty: every net with tracks but the references
    std::vector<std::string> net_names;       // exact names, besides the patterns
    std::vector<std::string> reference_nets;  // empty: nets that fill plane_fraction of a layer
    double plane_fraction = 0.25;
    std::vector<SiPair> pairs;
    bool pair_by_name = true;                 // pair nets whose names differ in P/N or +/- at the end
    std::vector<SiTarget> targets;
    double reference_hz = 1.0e9;              // frequency of the reported quantities
    std::map<std::string, DielectricModel> layer_materials;  // by stackup layer name
    DielectricModel default_dielectric{4.3, 0.02, 1.0e9, true};
    double mask_epsilon_r = 3.8, mask_loss_tangent = 0.025;
    ConductorModel copper;
    double etch_nm = 0.0;                     // each side of a conductor's far face is this much narrower
    std::map<std::string, bool> etch_bottom;  // per copper or stackup layer name: the narrow face is the bottom one
    double step_nm = 1.0e6;                   // largest distance between cuts along a track
    double lateral_nm = 25.0e6;               // largest half width of a cut (it reaches twenty times the distance to the shields)
    double neighbour_nm = 1.0e6;              // other nets' copper this close is part of the cut
    double large_copper_nm = 1.0e6;           // other nets' copper this wide counts as reference
    double target_error = 0.08;               // of each cross-section (impedance error about 0.2 %)
    double neighbour_heights = 3.0;           // quiet lines within this many shield distances
    double rise_time_s = 100e-12;             // of the aggressor step, for crosstalk
    double crosstalk_threshold = 0.05;        // near- or far-end crosstalk reported above this
    double min_finding_nm = 0.5e6;            // shorter runs of a finding are not reported
    double reference_merge = 0.1;             // reference pieces join across gaps below this fraction of their distance to the track (0: never)
    bool dispersion = true;                   // microstrip dispersion (Kirschning and Jansen)
    double finding_rise_s = 100e-12, finding_rise_fraction = 0.2;  // runs of less delay are electrically short
};

struct SiSample {
    int net = 0, layer = 0;
    double x_nm = 0, y_nm = 0;
    double distance_nm = 0;   // along the net from its start, -1 where unreached
    double z0 = 0;            // ohm, single-ended
    double z_differential = 0, z_common = 0;  // with the pair partner in the cut, else 0
    double epsilon_eff = 0, delay_ps_per_mm = 0, loss_db_per_m = 0;
    std::vector<int> reference_layers;
    std::size_t cross_section = 0;  // index into SiResult::cross_sections
    bool no_reference = false;
    std::size_t track = 0;          // index into Board::tracks
    double along_nm = 0;            // position along the track from its start
    bool coupled = false;           // the pair partner is conductor 1 of the cross-section
    double partner_x_nm = 0, partner_y_nm = 0;  // middle of the partner's copper in the cut
    // Other signal conductors of the cut (the partner and quiet neighbours): net, index
    // among the cross-section's conductors, middle of their copper.
    struct Other {
        int net = 0;
        std::size_t conductor = 0;
        double x_nm = 0, y_nm = 0;
    };
    std::vector<Other> others;
};

// Coupling of a victim net to an aggressor running beside it: the backward and forward
// coupling coefficients of their cross-sections (weak coupling), summed over the length
// they run together, and the near- and far-end crosstalk of a step of the given rise time.
struct SiCoupling {
    std::string aggressor, victim;
    int layer = 0;
    double length_nm = 0;          // coupled
    double delay_ps = 0;           // of the coupled length
    double kb_max = 0, kb_mean = 0;  // backward coefficient (1/4)(Lm/Ls + Cm/Cs)
    double kf_ps_per_mm = 0;       // forward coefficient (1/2)(Cm/Cs - Lm/Ls) sqrt(Ls Cs), mean
    double next = 0, fext = 0;     // fractions of the aggressor's step at the victim's ends
};

struct SiFinding {
    std::string net, kind;  // "impedance", "differential impedance", "no reference", "reference change", "uncoupled", "crosstalk"
    int layer = 0;
    double x_nm = 0, y_nm = 0, length_nm = 0;
    double value = 0, target = 0;
    std::string text;
    double delay_ps = 0;            // delay of the run
    bool electrically_short = false;  // delay below finding_rise_fraction of the rise time
};

struct SiNet {
    std::string name;
    std::string partner;          // pair partner, if any
    double length_nm = 0, delay_ps = 0;
    double z_min = 0, z_max = 0, z_mean = 0;  // single-ended, length-weighted mean
    double zd_min = 0, zd_max = 0, zd_mean = 0;  // differential, where coupled
    double coupled_nm = 0;        // length with the partner in the cut
    double loss_db = 0;           // at the reference frequency
    std::vector<int> layers;
    std::optional<SiTarget> target;
};

struct SiCrossSection {
    std::vector<std::string> conductors;  // roles: "self", "partner", "neighbour <net>"
    CsResult result;
};

struct SiResult {
    std::vector<SiNet> nets;
    std::vector<SiSample> samples;
    std::vector<SiCoupling> coupling;  // by decreasing crosstalk
    std::vector<SiFinding> findings;
    std::vector<SiCrossSection> cross_sections;
    std::size_t cuts = 0;
    double finding_rise_ps = 0, finding_rise_fraction = 0;  // of the electrically short findings
    double seconds = 0;
    std::vector<std::string> warnings;
};

SiResult analyse_si(const Board& board, const SiSettings& settings);

// Project file section `si`.
struct SiProject {
    std::vector<std::string> board;
    std::optional<SiSettings> si;
};
SiProject load_si_project(const std::string& path);
SiProject parse_si_project(std::string_view yaml, const std::string& base_dir);
std::string si_summary(const Board& board, const SiResult& result);
// The result as JSON (millimetres): nets, samples, the analysed tracks with their samples
// and the top of their copper, findings, cross-sections.
std::string si_json(const Board& board, const SiResult& result);

}  // namespace kairoq
