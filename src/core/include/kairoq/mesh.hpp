// Triangle meshes of extruded polygons, for display.
#pragma once

#include "kairoq/board.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace kairoq {

struct Mesh {
    std::vector<float> positions;   // x, y, z in mm, relative to the mesh origin
    std::vector<std::uint32_t> indices;
    std::vector<std::uint32_t> ids;  // per vertex: net id, or 0
};

// Triangulation of a polygon with holes, as indices into its rings concatenated
// (outer first).
std::vector<std::uint32_t> triangulate(const Polygon& polygon);

// Prism over each polygon between z0 and z1 (nm), with top, bottom and side walls.
// Coordinates are shifted by -origin so that float precision suffices.
void extrude(Mesh& mesh, const PolygonSet& polygons, double z0_nm, double z1_nm, Point64 origin, std::uint32_t id);

// Centre of the board outline's bounding box, used as the display origin.
Point64 display_origin(const Board& board);

// Display meshes of a board, keyed "copper/<index>", "dielectric/<stackup position>",
// "mask/<side>", "drills".
std::vector<std::string> mesh_keys(const Board& board);
Mesh board_mesh(const Board& board, const std::string& key);

// Binary encoding served to the viewer: uint32 vertex count, uint32 index count,
// float32 positions, uint32 indices, uint32 ids (all little endian).
std::string encode_mesh(const Mesh& mesh);

}  // namespace kairoq
