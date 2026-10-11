// The model of a loop as the loop solver builds it: the window, the copper islands in it,
// the plated holes and the elements. The openEMS export writes the same geometry.
#pragma once

#include "kairoq/emi.hpp"

#include <map>
#include <string>
#include <vector>

namespace kairoq::detail {

struct LoopGeometry {
    PolygonSet window;                           // nm, within the board outline
    std::map<int, std::vector<Polygon>> copper;  // layer -> islands, simplified away from the loop
    std::vector<double> z_m;                     // mid-plane of each copper layer (index 1..N), m
    struct Via {
        Point64 at;
        double radius_nm;
        int from_layer, to_layer;  // outermost layers it joins in the model
    };
    std::vector<Via> vias;
    struct Element {
        std::string from_pad, to_pad;
        Point64 from, to;  // pad centres, nm
        int layer;         // the side it lies on: 1 or N
        double z_m, radius_m;  // axis height, radius of the round conductor
        bool stitch = false;    // a stitch, after the loop's elements
        double capacitance_f = 0;  // a stitch's series capacitance (0: none)
    };
    std::vector<Element> elements;
};

LoopGeometry loop_geometry(const Board& board, const LoopSpec& loop, const LoopSettings& settings);

}  // namespace kairoq::detail
