// Capacitance of a net's copper to a plate on a face of the board and to the other nets
// around it, by the method of moments in the stratified medium of the stack (internal).
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/emi.hpp"

#include <string>
#include <utility>
#include <vector>

namespace kairoq::detail {

struct NodeCapacitance {
    bool solved = false;      // false: no copper of the node over the plate, or too many panels
    double to_plate_f = 0.0;  // the plate's share of the node's charge at 1 V
    double total_f = 0.0;     // the node's charge at 1 V, all other conductors at 0 V
    // To the plate on the mesh used and on one with elements twice as large, before the
    // extrapolation from the two (0: no second mesh).
    double unextrapolated_to_plate_f = 0.0, coarse_to_plate_f = 0.0;
    std::vector<std::pair<int, double>> to_nets;  // net, F
    std::size_t panels = 0;
    double window_margin_m = 0.0, seconds = 0.0;
    std::vector<std::string> warnings;
};

NodeCapacitance node_capacitance(const Board& board, int net, const PlateSpec& plate, const EmiSettings& settings);

}  // namespace kairoq::detail
