// The canonical board model. Every importer produces it; docs/canonical-model.md
// describes it and its file format.
#pragma once

#include "kairoq/geometry.hpp"

#include <optional>
#include <string>
#include <vector>

namespace kairoq {

constexpr int kModelVersion = 2;  // 2: tracks

enum class LayerType { Copper, Dielectric, SolderMask, Silkscreen, SolderPaste, Other };
enum class Side { Top, Inner, Bottom };
enum class DrillFunction { Via, Component, Mechanical, Unknown };

const char* to_string(LayerType t);
const char* to_string(Side s);
const char* to_string(DrillFunction f);
LayerType layer_type_from_string(const std::string& s);
Side side_from_string(const std::string& s);
DrillFunction drill_function_from_string(const std::string& s);

// One physical layer, ordered from the top of the board to the bottom.
struct StackupLayer {
    std::string name;
    LayerType type = LayerType::Other;
    double thickness_nm = 0.0;
    std::string material;
    std::optional<double> epsilon_r;
    std::optional<double> loss_tangent;
    std::optional<double> conductivity_s_per_m;  // copper layers
    int copper_index = 0;  // 1..N for copper layers, 0 otherwise
};

// Merged copper of one net on one layer.
struct NetCopper {
    int net = 0;
    PolygonSet polygons;
};

struct CopperLayer {
    int index = 0;  // 1 = top
    std::string name;
    Side side = Side::Inner;
    std::vector<NetCopper> nets;  // ascending net id
};

// Solder mask openings of one side.
struct MaskLayer {
    Side side = Side::Top;
    PolygonSet openings;
};

struct Pad {
    std::string component;
    std::string pin;
    std::string function;  // ComponentPad, SMDPad, ...
    int net = 0;
    int layer = 0;  // copper layer index
    Point64 at{0, 0};
    PolygonSet shape;
};

struct Drill {
    Point64 at{0, 0};
    std::optional<Point64> slot_end;
    double diameter_nm = 0.0;
    bool plated = true;
    DrillFunction function = DrillFunction::Unknown;
    int from_layer = 0;  // copper layer indices, from <= to
    int to_layer = 0;
    int net = 0;
};

// A drawn track: a straight or circular stroke of a round tool, as the design data
// gives it. Its copper is also part of the layer's merged copper.
struct Track {
    int layer = 0;  // copper layer index
    int net = 0;
    double width_nm = 0.0;
    Point64 start{0, 0}, end{0, 0};
    std::optional<Point64> centre;  // arcs: centre, direction from start to end
    bool ccw = false;
};

struct Component {
    std::string refdes;
    std::string value;
    std::string package;
    Point64 at{0, 0};
    double rotation_deg = 0.0;
    Side side = Side::Top;
    std::string position_source;  // "placement", "ipc2581", "pads"
};

struct SourceInfo {
    std::string importer;   // "fabrication" or "ipc2581"
    std::string generator;  // generating software, when the files name it
    std::vector<std::string> files;
    double arc_tolerance_nm = 0.0;
};

struct Board {
    std::string name;
    SourceInfo source;
    PolygonSet outline;
    std::string stackup_source;  // "gbrjob", "ipc2581", "yaml", "default"
    std::vector<StackupLayer> stackup;
    std::vector<CopperLayer> copper;  // index order, top first
    std::vector<MaskLayer> masks;
    std::vector<std::string> nets;  // id -> name; id 0 is "no net" with an empty name
    std::string net_source;         // "x2", "ipc2581", "ipc356", "geometric"
    std::vector<Component> components;
    std::vector<Pad> pads;
    std::vector<Drill> drills;
    std::vector<Track> tracks;
    std::vector<std::string> warnings;

    int net_id(const std::string& name);  // adds the net when missing
    int find_net(const std::string& name) const;  // -1 when missing
    const StackupLayer* stackup_for_copper(int index) const;
    double thickness_nm() const;
};

// Vertical position of each stackup layer: z = 0 at the bottom surface of the
// bottom layer, z increasing upwards.
struct LayerZ {
    double bottom_nm = 0.0;
    double top_nm = 0.0;
};
std::vector<LayerZ> stackup_z(const Board& board);

// The part of a board inside a rectangle (nm): outline, copper, mask openings and pads
// clipped to it; holes, tracks and components kept when they lie wholly inside (a component
// when its origin does). The cut edges are open: no copper or heat crosses them.
Board crop_board(const Board& board, const Rect64& area);

}  // namespace kairoq
