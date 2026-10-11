// Parts of the physics-based via model (docs/solvers/networks.md): the capacitances of a
// barrel and its pads to the reference planes around it, from an axisymmetric
// electrostatic solution, and the impedance of a pair of parallel planes seen from a
// barrel between them (a radial line).
#pragma once

#include <complex>
#include <vector>

namespace kairoq {

using cplx = std::complex<double>;

// Hankel functions of the second kind of orders 0 and 1, H = J - jY, for Re z > 0 and
// Im z <= 0 (outgoing waves in a lossy medium with time dependence exp(j omega t)).
void hankel2(cplx z, cplx& h0, cplx& h1);

// A via's surroundings, along its axis (z up, nm).
struct ViaCut {
    double drill_nm = 0;
    struct Copper {
        double z_bottom_nm = 0, z_top_nm = 0;
        double pad_nm = 0;      // pad diameter on this layer, 0 for none
        double antipad_nm = 0;  // a reference plane on this layer: the diameter of its hole; 0 for none
        bool solid = false;     // a reference plane without a hole (beyond the barrel's ends)
    };
    std::vector<Copper> copper;  // top to bottom
    struct Dielectric {
        double z_bottom_nm = 0, z_top_nm = 0, epsilon_r = 1.0;
    };
    std::vector<Dielectric> dielectrics;
    std::size_t from = 0, to = 0;  // the barrel's end layers (indices into copper)
};


// The barrel as a ladder: nodes along it (each copper layer it crosses, and slices of
// longer gaps), segments between them, the Maxwell capacitance matrix of the nodes (the
// barrel cut midway between nodes, with the pads of each node's layer, all planes
// grounded), and the inductance of the segments that the pair's higher modes and, outside
// the pairs, the wire and its image in the nearest plane give. The pairs' fundamental
// (TEM) mode is left to the caller, who couples it between vias.
struct ViaLadder {
    struct Node {
        double z_nm = 0;
        int layer = -1;  // index into ViaCut::copper, -1 for a slice
    };
    struct Segment {
        std::size_t upper = 0, lower = 0;  // node indices
        int plane_upper = -1, plane_lower = -1;  // copper indices of the pair it lies in; -1 outside
        double weight = 0;  // its share of the pair's height (TEM coupling)
    };
    std::vector<Node> nodes;
    std::vector<Segment> segments;
    std::vector<double> capacitance;  // nodes x nodes, F
    std::vector<double> inductance;   // segments x segments, H
};
ViaLadder via_ladder(const ViaCut& cut, double max_segment_nm, double target_error = 0.02);

// Infinite plane pair of separation h (m), wavenumber k (1/m, complex with loss): the
// impedance between the planes at a barrel of radius a (m), and the transfer impedance to
// a point at distance rho (m).
cplx radial_self(double omega, double h, cplx k, double a);
cplx radial_transfer(double omega, double h, cplx k, double rho);

}  // namespace kairoq
