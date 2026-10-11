#include "layered_green.hpp"

#include <doctest/doctest.h>

#include <cmath>

using namespace kairoq::detail;

namespace {
constexpr double kEps0 = 8.8541878128e-12;
constexpr double kPi = 3.14159265358979323846;
}  // namespace

TEST_CASE("Layered Green's function: a uniform medium") {
    Stratified s;
    s.ground = false;
    s.thickness = {0.4e-3, 1.0e-3, 0.2e-3};
    s.epsilon_r = {4.0, 4.0, 4.0};
    s.epsilon_below = s.epsilon_above = 4.0;
    for (double z : {0.4e-3, 0.9e-3}) {
        const LayeredGreen g(s, z, 0.4e-3, 0.02);
        const double dz = z - 0.4e-3;
        for (double rho : {0.0, 1e-5, 3e-4, 2e-3, 1.5e-2}) {
            const double R = std::hypot(rho, dz);
            if (R == 0) continue;
            CHECK(g(rho) == doctest::Approx(1.0 / (4 * kPi * kEps0 * 4.0 * R)).epsilon(1e-6));
        }
    }
}

TEST_CASE("Layered Green's function: a charge on a slab over a ground plane, by images") {
    // G = 1 / (2 pi eps0 (er + 1)) [1/rho - (1 + K) sum_n (-K)^(n-1) / sqrt(rho^2 + (2 n h)^2)],
    // K = (er - 1) / (er + 1), on the slab's surface.
    const double h = 1.5e-3, er = 4.4, K = (er - 1) / (er + 1);
    Stratified s;
    s.thickness = {h};
    s.epsilon_r = {er};
    const LayeredGreen g(s, h, h, 0.05);
    for (double rho : {2e-5, 2e-4, 1e-3, 5e-3, 3e-2}) {
        double sum = 1.0 / rho;
        for (int n = 1; n < 4000; ++n) sum -= (1 + K) * std::pow(-K, n - 1) / std::hypot(rho, 2 * n * h);
        const double ref = sum / (2 * kPi * kEps0 * (er + 1));
        CHECK(g(rho) == doctest::Approx(ref).epsilon(2e-4));
    }
}

TEST_CASE("Layered Green's function: reciprocity between two heights") {
    Stratified s;
    s.thickness = {0.2e-3, 0.035e-3, 1.0e-3, 0.035e-3, 0.02e-3};
    s.epsilon_r = {1.0, 4.3, 4.3, 4.3, 3.8};
    const LayeredGreen a(s, 0.2175e-3, 1.2525e-3, 0.03), b(s, 1.2525e-3, 0.2175e-3, 0.03);
    for (double rho : {0.0, 1e-4, 1e-3, 1e-2})
        CHECK(a(rho) == doctest::Approx(b(rho)).epsilon(1e-6));
}
