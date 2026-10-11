// Open any supported input as a Board.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/fabrication.hpp"

#include <string>
#include <vector>

namespace kairoq {

// A single canonical model (directory or zip) is loaded; a single IPC-2581 file is
// imported; anything else is imported as a set of fabrication files.
Board open_board(const std::vector<std::string>& inputs, const ImportOptions& options = {});

// Human-readable summary: stackup, layers, nets, drills, warnings.
std::string board_summary(const Board& board);

}  // namespace kairoq
