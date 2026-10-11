// Planar triangle meshes for the finite element solvers.
#pragma once

#include "kairoq/geometry.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

namespace kairoq {

struct MeshSettings {
    // Outlines are simplified so that no point moves by more than this (nm).
    double boundary_tolerance_nm = 2000.0;
    // Element size field: h(x) = min(max_edge, min_size + grading * d(x)), with d the
    // distance to the nearest boundary (copper edges, holes, ports). Elements are small
    // where the current crowds and grow away from the edges.
    double max_edge_nm = 1.0e6;
    double min_size_nm = 50000.0;
    double grading = 0.3;
    // Smallest triangle angle the refinement aims for (degrees, at most about 30).
    double min_angle_deg = 20.0;
    // Edges this short are not split further (nm): limits refinement at sharp input corners.
    double min_edge_nm = 500.0;
    // Budget of inserted vertices per region.
    std::size_t max_vertices = 4'000'000;
    // Size at the outer ring of the region and at constraint edges (nm); 0: the minimum
    // size. A coarse size suits a far boundary and long interfaces between materials,
    // which need no refinement of their own.
    double outer_size_nm = 0.0;
    double constraint_size_nm = 0.0;
    // Size at a boundary segment from the position of its middle (nm); when set, it
    // replaces the minimum size and outer_size_nm on the region's outline.
    std::function<double(double, double)> boundary_size;
};

struct TriMesh {
    std::vector<double> x, y;  // nm
    std::vector<std::array<std::uint32_t, 3>> triangles;  // counter-clockwise
    std::size_t unrefined = 0;  // refinements the mesher could not perform
};

// Quality mesh of a polygon with holes. The boundaries of `constraints` (typically
// ports) inside the region become mesh edges, so that a constraint polygon is an
// exact union of triangles.
// `seeds` are additional vertices {x, y, spacing} (nm), for example from adaptive
// refinement. Seeds outside the region, or closer to its boundary than 0.3 times their
// spacing, are ignored.
TriMesh mesh_region(const Polygon& region, const PolygonSet& constraints, const MeshSettings& settings,
                    const std::vector<std::array<double, 3>>& seeds = {});

// The region after outline simplification, as meshed by mesh_region.
PolygonSet simplified_region(const Polygon& region, const MeshSettings& settings);

double triangle_area(const TriMesh& mesh, std::size_t t);  // nm^2, positive
double min_angle_deg(const TriMesh& mesh, std::size_t t);

}  // namespace kairoq
