// Readers for the auxiliary fabrication files, and the importer that combines a set of
// fabrication files (Gerber X2, drill, job file, netlist, placement) into a Board.
#pragma once

#include "kairoq/board.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

// ---------------------------------------------------------------- Gerber job file

struct JobFileEntry {
    std::string path;
    std::vector<std::string> function;  // .FileFunction fields
    std::string polarity;
};

struct JobFile {
    std::string generator;
    std::string project;
    int layer_count = 0;
    double board_thickness_nm = 0.0;
    std::vector<JobFileEntry> files;
    std::vector<StackupLayer> stackup;  // top to bottom; empty without MaterialStackup
};

JobFile parse_job_file(std::string_view json);

// ---------------------------------------------------------------- IPC-D-356 netlist

struct NetlistPoint {
    std::string net;
    std::string refdes;
    std::string pin;
    Point64 at{0, 0};
    double hole_nm = 0.0;  // 0 for surface pads
    bool plated = false;
    int access = 0;  // 0: both sides, otherwise the copper layer number
    bool midpoint = false;
};

struct Netlist {
    std::vector<NetlistPoint> points;
    std::vector<std::string> warnings;
};

Netlist parse_ipc356(std::string_view text);
bool looks_like_ipc356(std::string_view text);

// ---------------------------------------------------------------- placement (pick and place)

struct PlacementRow {
    std::string refdes, value, package;
    Point64 at{0, 0};
    double rotation_deg = 0.0;
    Side side = Side::Top;
};

// CSV with a header row. Columns are recognised by name: Ref/Designator,
// Val/Value/Comment, Package/Footprint, PosX/Mid X/X, PosY/Mid Y/Y, Rot/Rotation,
// Side/Layer. Lengths are millimetres unless a column name or value says mils or in.
std::vector<PlacementRow> parse_placement_csv(std::string_view text);
bool looks_like_placement_csv(std::string_view text);

// ---------------------------------------------------------------- stackup YAML

// layers: [{name, type, thickness_mm, material, epsilon_r, loss_tangent}], top to bottom.
std::vector<StackupLayer> parse_stackup_yaml(std::string_view text);

// ---------------------------------------------------------------- importer

struct ImportOptions {
    GeometryOptions geometry;
    std::string stackup_yaml_path;  // overrides the job file stackup when set
    // Board edge loops whose ends are closer than this are joined (nm).
    double outline_join_nm = 2000.0;
    // The board cropped to this rectangle (nm) after reading (crop_board), when set.
    std::optional<Rect64> crop;
};

// Import a set of fabrication files. Each input is a file or a directory (its files,
// not recursive). File roles come from X2 .FileFunction attributes, the job file, the
// file contents and, last, the file name extension.
Board import_fabrication(const std::vector<std::string>& inputs, const ImportOptions& options = {});

}  // namespace kairoq
