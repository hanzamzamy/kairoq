// Drill data: Excellon (IPC-NC-349 with the KiCad/X2 attribute comments) and Gerber
// drill files (X2 .FileFunction Plated/NonPlated).
#pragma once

#include "kairoq/gerber.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

struct DrillTool {
    int number = 0;
    double diameter_nm = 0.0;
    std::string function;  // X2 drill function: ViaDrill, ComponentDrill, MechanicalDrill, ...
};

struct DrillHit {
    int tool = 0;  // index into DrillFile::tools
    Point64 at{0, 0};
    std::optional<Point64> slot_end;  // routed or G85 slot
};

struct DrillFile {
    // From .FileFunction: Plated,<from>,<to>,<PTH|Blind|Buried> or NonPlated,...
    std::optional<bool> plated;
    int from_layer = 0;  // copper layer numbers, 0 when not given
    int to_layer = 0;
    std::string span_type;  // PTH, NPTH, Blind, Buried
    std::vector<DrillTool> tools;
    std::vector<DrillHit> hits;
    std::vector<std::string> warnings;
};

DrillFile parse_excellon(std::string_view text);
bool looks_like_excellon(std::string_view text);

// Drill file from a Gerber image: flashes are holes, strokes are slots.
DrillFile drill_from_gerber(const GerberImage& image);

// Apply a .FileFunction (Plated,1,4,PTH,...) to a drill file.
void apply_drill_file_function(DrillFile& drill, const std::vector<std::string>& fields);

}  // namespace kairoq
