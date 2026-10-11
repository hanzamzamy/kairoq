// Electrostatic Green's function of a stratified medium: horizontal slabs of dielectric
// between a perfectly conducting plane below and a half-space above, or between two
// half-spaces. The potential of a unit point charge, from its Hankel transform, with the
// quasi-static direct term extracted in closed form (internal). Lengths in metres.
#pragma once

#include <cstddef>
#include <vector>

namespace kairoq::detail {

struct Stratified {
    // Slabs from the bottom up: thickness (m) and relative permittivity. The bottom of the
    // first slab is at z = 0.
    std::vector<double> thickness, epsilon_r;
    bool ground = true;           // a conducting plane at z = 0; else a half-space below
    double epsilon_below = 1.0;   // the half-space below, without the ground
    double epsilon_above = 1.0;   // the half-space above the last slab
};

// Hankel transform of the potential (V) at height z of a unit charge (C) at height zs, times
// exp(k |z - zs|) so that it stays finite as k grows.
double spectral_scaled(const Stratified& s, double k, double z, double zs);

// Potential of a unit point charge at height zs, observed at horizontal distance rho and
// height z: G = C / (2 pi R) + sum_m c_m / (2 pi R_m) + remainder(rho), with
// C = lim k->inf k phi(k) exp(k |z - zs|) (F^-1) and R the distance. In a plane (z = zs) the
// images c_m at depths 2 m d, d the distance to the nearest interface, are fitted to the
// transform so that the remainder, tabulated against rho and interpolated, varies on scales
// larger than d: an interface closer than the elements would otherwise put a near-singular
// term into it.
class LayeredGreen {
public:
    LayeredGreen(const Stratified& s, double z, double zs, double rho_max);
    double direct_coefficient() const { return c_; }  // C / (2 pi), the coefficient of 1/R
    struct Image {
        double coefficient = 0, depth = 0;  // c_m / (2 pi) and the added vertical distance (m)
    };
    const std::vector<Image>& images() const { return images_; }
    double remainder(double rho) const;               // V per C
    double operator()(double rho) const;              // V per C, the whole potential
    double dz() const { return dz_; }

private:
    double c_ = 0, dz_ = 0, rho_min_ = 0;
    std::vector<Image> images_;
    // The remainder times sqrt(rho^2 + a^2), a the stack's height, which tends to constants at
    // both ends; cubic interpolation in log rho.
    double a_ = 0;
    std::vector<double> log_rho_, rem_;
    double rem0_ = 0;  // at rho = 0
};

}  // namespace kairoq::detail
