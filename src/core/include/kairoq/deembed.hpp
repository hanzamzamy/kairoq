// De-embedding of S-parameters: change of the reference impedance, lines added at or
// removed from ports, and fixtures removed by the half of a symmetric 2x-thru. Method:
// docs/solvers/rf.md.
#pragma once

#include "kairoq/network.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

// A uniform quasi-TEM line: impedance, effective permittivity, attenuation at 1 GHz
// scaling with the square root of frequency (conductor loss), length (negative: removed).
struct LineFixture {
    double z_ohm = 50.0, epsilon_eff = 1.0, loss_db_per_m_1ghz = 0.0, length_m = 0.0;
};

// S with every port's reference changed to z_ohm (real references).
Touchstone renormalise(const Touchstone& t, double z_ohm);

// S with a two-port in front of a port (0-based): its port 2 joins the network, its port 1
// is the new port. fixture_s per frequency of t, row major (S11 S12 S21 S22), at t's reference.
Touchstone connect_fixture(const Touchstone& t, std::size_t port, const std::vector<std::array<cplx, 4>>& fixture_s);

Touchstone add_line(const Touchstone& t, std::size_t port, const LineFixture& line);

// ABCD matrices of the half of a symmetric 2x-thru at the given frequencies: the square
// root of the thru's ABCD matrix, the branch followed from low frequency.
std::vector<std::array<cplx, 4>> half_thru_abcd(const Touchstone& thru, const std::vector<double>& frequency_hz);

// Removes the half of a 2x-thru from a port: its port 1 at the port, its port 2 the new
// reference plane.
Touchstone deembed_thru(const Touchstone& t, std::size_t port, const Touchstone& thru);

// Touchstone version 1 text (real and imaginary parts, Hz).
std::string touchstone_text(const Touchstone& t, const std::string& comment);

// A de-embedding project: `input` (a Touchstone file), `operations` in order
// (`renormalise_ohm`, `line`, `thru`), and `output` (a file to write, optional).
struct DeembedResult {
    Touchstone data;
    std::string output;  // path written, if any
    std::vector<std::string> steps;
};
DeembedResult deembed(std::string_view yaml, const std::string& base_dir = ".");

}  // namespace kairoq
