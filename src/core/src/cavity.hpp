// A pair of parallel planes as a two-dimensional cavity: the voltage between the planes
// over the region where both exist, with an open (magnetic) wall at its outline, excited
// by the barrels that cross it (ports). Its modes come from finite elements on the
// region; the impedance between the ports at any frequency is their sum, with the static
// part of the modes not computed. Description: docs/solvers/pdn.md.
#pragma once

#include "kairoq/geometry.hpp"

#include <complex>
#include <cstddef>
#include <vector>

namespace kairoq {

struct CavityPort {
    double x_nm = 0, y_nm = 0, radius_nm = 0;
    bool shorted = false;  // a barrel joining the two planes: no voltage across the pair there
};

class Cavity {
public:
    // Modes up to the wavenumber k_max (1/m, in the dielectric) are computed; the mesh
    // resolves them. Ports must lie inside the region and must not overlap.
    Cavity(const PolygonSet& region, const std::vector<CavityPort>& ports, double k_max);

    // Impedance matrix between the open ports (row major, in the order given, shorted ports
    // left out) for planes h apart (m) filled with relative permittivity eps_r (complex,
    // eps' - j eps'').
    std::vector<std::complex<double>> impedance(double omega, double h_m, std::complex<double> eps_r) const;
    // The same with lossy plates: z_surface the surface impedance of each plate (ohm), which
    // adds 2 z_surface to the pair's series impedance j omega mu h per square.
    std::vector<std::complex<double>> impedance(double omega, double h_m, std::complex<double> eps_r,
                                                std::complex<double> z_surface) const;

    std::size_t modes() const { return lambda_.size(); }
    std::size_t nodes() const { return nodes_; }
    std::size_t lanczos_steps() const { return steps_; }
    const std::vector<double>& eigenvalues() const { return lambda_; }  // 1/m^2, ascending
    double area_m2() const { return area_; }

private:
    std::size_t open_ = 0, nodes_ = 0, steps_ = 0;
    double area_ = 0;
    std::vector<double> lambda_;     // kept modes
    std::vector<double> mean_;       // modes x open ports: each mode's mean over the port circle
    std::vector<double> remainder_;  // open x open: static sum of the modes not kept (m^2... per 1/m^2)
};

}  // namespace kairoq
