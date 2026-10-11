// Reading and writing the canonical board model (docs/canonical-model.md).
#pragma once

#include "kairoq/board.hpp"

#include <nlohmann/json_fwd.hpp>

#include <string>

namespace kairoq {

// Write to a directory, or to a zip archive when the path ends in ".zip".
void save_board(const Board& board, const std::string& path);
Board load_board(const std::string& path);

// True if the path holds a canonical model (a directory with board.json, or a zip
// archive containing one).
bool is_board_model(const std::string& path);

// board.json content. Geometry is referenced by file and byte range; with
// geometry = false the references are left out (summaries).
nlohmann::json board_json(const Board& board, bool geometry);

// OGC Well-Known Binary MultiPolygon (little endian, coordinates in mm).
std::string encode_wkb(const PolygonSet& polygons);
PolygonSet decode_wkb(std::string_view data);

}  // namespace kairoq
