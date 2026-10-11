// Kernels of the partial-inductance (PEEC) solver: potentials of triangles, partial
// inductances of straight conductors, and dense symmetric positive definite algebra
// (internal). Lengths in metres.
#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace kairoq::detail {

constexpr double kMu0Over4Pi = 1e-7;  // H/m

struct Vec2 {
    double x = 0, y = 0;
};

// Integrals over a triangle in a horizontal plane (vertices a, b, c, counter-clockwise) of
// 1/R and of (rho - rho0)/R, with R the distance to the point (rho0, height d above or
// below the plane). Closed forms after Wilton et al. (IEEE Trans. AP-32, 1984) and Graglia
// (IEEE Trans. AP-41, 1993).
struct TrianglePotential {
    double scalar = 0;  // integral of 1/R (m)
    Vec2 vector;        // integral of (rho - rho0)/R (m^2)
};
TrianglePotential triangle_potential(const std::array<Vec2, 3>& v, Vec2 rho0, double d);

// Quadrature points on a triangle: barycentric coordinates and weights summing to 1.
struct TriangleRule {
    std::vector<std::array<double, 3>> bary;
    std::vector<double> weight;
};
const TriangleRule& triangle_rule_7();   // degree 5 (Dunavant)
const TriangleRule& triangle_rule_28();  // the 7-point rule on the 4 sub-triangles

// Partial mutual inductance (H) of two parallel straight filaments along a common axis
// direction, at perpendicular distance rho, spanning [a1, b1] and [a2, b2] along it.
// With rho the radius and equal spans this is the partial self-inductance of a tube.
double parallel_filaments(double a1, double b1, double a2, double b2, double rho);

// Partial mutual inductance (H) of two straight segments in space, by Gauss quadrature
// subdivided by their distance (not for overlapping segments).
double segment_mutual(const std::array<double, 3>& p0, const std::array<double, 3>& p1, const std::array<double, 3>& q0,
                      const std::array<double, 3>& q1);

// Dense symmetric positive definite matrix, row-major n x n, factorised in place
// (lower triangle, A = C C^T). Returns false when not positive definite.
bool cholesky(std::vector<double>& a, std::size_t n);
// Solves C y = b (forward) or C^T x = y (backward) for several right-hand sides, stored as
// columns of a row-major n x m matrix.
void forward_solve(const std::vector<double>& c, std::size_t n, std::vector<double>& b, std::size_t m);
void backward_solve(const std::vector<double>& c, std::size_t n, std::vector<double>& b, std::size_t m);

}  // namespace kairoq::detail
