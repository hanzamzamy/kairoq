// Steps shared by the importers once the layers are rendered.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/fabrication.hpp"
#include "kairoq/feature.hpp"

#include <map>
#include <string>
#include <vector>

namespace kairoq::detail {

// Evenly divided FR-4 stackup with 35 um copper, for inputs without stackup data.
std::vector<StackupLayer> default_stackup(int copper_layers, double board_thickness_nm);

// Store per-net copper on a layer, interning net names.
void set_layer_copper(Board& board, CopperLayer& layer, const std::map<std::string, PolygonSet>& by_net);

// Net of each plated drill from the copper at its centre on the layers it spans.
void assign_drill_nets(Board& board);

// Nets from geometric connectivity: copper islands joined by plated drills, named
// from the netlist points where available. Replaces the copper of every layer.
void geometric_nets(Board& board, const std::vector<PolygonSet>& copper, const Netlist* netlist);

// Places an IPC-D-356 netlist on the board: generators write it in another origin, or with
// y upwards. The translation (and mirror of y) under which most of its holes land on
// plated drills is applied; returns false when none places a majority of them.
bool register_netlist(Netlist& netlist, const Board& board, std::string& note);

// Compare pad nets with the netlist; mismatches become warnings.
void cross_check_netlist(Board& board, const Netlist& netlist);

// Components located at the mean of their pads, for inputs without placement data.
void components_from_pads(Board& board);

// Tracks from the strokes of each copper layer (index 0 = layer 1), with the net of the
// copper at their middle. Call once the layers' copper and nets are final.
void collect_tracks(Board& board, const std::vector<const std::vector<Feature>*>& layers);

// Warn about copper of different nets overlapping on a layer.
void check_net_overlaps(Board& board);

std::string normalise_net_name(const std::string& name);

}  // namespace kairoq::detail
