#include "cavity.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>

using namespace kairoq;

namespace {

constexpr double kPiT = 3.14159265358979323846;
constexpr double kMu0T = 1.25663706212e-6;
constexpr double kC0T = 299792458.0;

PolygonSet rectangle(double w_mm, double h_mm) {
    const std::int64_t w = std::llround(w_mm * 1e6), h = std::llround(h_mm * 1e6);
    const std::int64_t o = 0;
    const Path64 p{Point64(o, o), Point64(w, o), Point64(w, h), Point64(o, h)};
    return PolygonSet{Polygon{p, {}}};
}

// Transfer impedance of a rectangular cavity with open walls between two points, by the
// double sum of its modes (point ports, which the transfer term allows).
std::complex<double> rectangle_transfer(double a, double b, double h, std::complex<double> er, double f, double x1, double y1, double x2,
                                        double y2) {
    const std::complex<double> k2 = std::pow(2 * kPiT * f / kC0T, 2) * er;
    std::complex<double> s = 0;
    for (int m = 0; m < 3000; ++m)
        for (int n = 0; n < 800; ++n) {
            const double km = m * kPiT / a, kn = n * kPiT / b;
            const double e = (m ? 2.0 : 1.0) * (n ? 2.0 : 1.0);
            s += e * std::cos(km * x1) * std::cos(kn * y1) * std::cos(km * x2) * std::cos(kn * y2) / (km * km + kn * kn - k2);
        }
    return std::complex<double>(0, 2 * kPiT * f * kMu0T * h / (a * b)) * s;
}

}  // namespace

TEST_CASE("Rectangular plane pair: modes, plate capacitance and transfer impedance") {
    const double er = 4.5, h = 0.2e-3, f_max = 10e9;
    const double kmax = 2 * kPiT * f_max / kC0T * std::sqrt(er);
    const std::vector<CavityPort> ports{{12e6, 3e6, 0.15e6, false}, {27e6, 6e6, 0.15e6, false}};
    const Cavity c(rectangle(40, 10), ports, kmax);
    // Modes against (m pi / a)^2 + (n pi / b)^2: the ports' holes perturb them slightly.
    std::vector<double> exact;
    for (int m = 0; m < 30; ++m)
        for (int n = 0; n < 10; ++n) exact.push_back(std::pow(m * kPiT / 0.04, 2) + std::pow(n * kPiT / 0.01, 2));
    std::sort(exact.begin(), exact.end());
    REQUIRE(c.modes() > 10);
    for (std::size_t k = 1; k < 20; ++k) CHECK(c.eigenvalues()[k] == doctest::Approx(exact[k]).epsilon(k < 10 ? 2e-3 : 5e-3));
    // Low frequency: the plate capacitance eps A / h.
    const double w = 2 * kPiT * 1e6;
    const auto z = c.impedance(w, h, er);
    const double cplate = 8.8541878128e-12 * er * 0.04 * 0.01 / h;
    CHECK(std::abs(1.0 / (std::complex<double>(0, w) * z[1])) == doctest::Approx(cplate).epsilon(5e-3));
    // Transfer impedance between the ports against the double sum, with a loss tangent of
    // 0.02, across the band: within 4 % (and 0.01 ohm), most of it the modes' frequency
    // error, which the resonances' Q of about 50 magnifies near them.
    for (double f = 0.5e9; f < 10e9; f += 0.5e9) {
        const std::complex<double> lossy(er, -er * 0.02);
        const auto zf = c.impedance(2 * kPiT * f, h, lossy);
        const auto ref = rectangle_transfer(0.04, 0.01, h, lossy, f, 0.012, 0.003, 0.027, 0.006);
        CHECK(std::abs(zf[1] - ref) < 0.04 * std::abs(ref) + 0.01);
    }
}
