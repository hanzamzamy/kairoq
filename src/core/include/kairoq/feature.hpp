// Graphical objects of one layer, as produced by the importers, and their rendering
// into merged copper per net.
#pragma once

#include "kairoq/geometry.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace kairoq {

enum class FeatureKind { Flash, Stroke, Region };

// One graphical object. The shape is its filled area; strokes also keep their centre
// line, from which board outlines are reconstructed.
struct Feature {
    FeatureKind kind = FeatureKind::Flash;
    bool dark = true;  // false: clear polarity, removes material drawn before it
    Paths64 shape;

    bool has_net = false;
    std::string net;        // X2 .N or IPC-2581 Set net; empty for "no net"
    bool pad = false;       // X2 .P present, or IPC-2581 PinRef; the pin name may be empty
    std::string component;  // X2 .P / .C, IPC-2581 PinRef componentRef
    std::string pin;
    std::string pin_function;
    std::string function;   // aperture function (X2 .AperFunction) or IPC-2581 usage

    Point64 at{0, 0};   // flash position or stroke start
    Point64 end{0, 0};  // stroke end
    std::optional<Point64> arc_centre;
    bool arc_ccw = false;
    double stroke_width_nm = 0.0;  // circular aperture diameter of a stroke

    bool is_pad() const { return pad; }
};

// Copper merged per net name, with clear-polarity objects applied in order.
std::map<std::string, PolygonSet> render_by_net(const std::vector<Feature>& features);

// All dark material merged, nets ignored.
PolygonSet render_all(const std::vector<Feature>& features);

// Closed loops formed by the stroke centre lines (and region contours), filled with
// the even-odd rule: a board outline with its cut-outs. Loose ends are reported.
struct OutlineResult {
    PolygonSet polygons;
    int open_chains = 0;
};
OutlineResult outline_from_features(const std::vector<Feature>& features, double tolerance_nm,
                                    double join_distance_nm);

}  // namespace kairoq
