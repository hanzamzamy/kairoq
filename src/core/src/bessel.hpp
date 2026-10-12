// Bessel functions of order zero for real arguments (docs/methods/networks.md, section 4.9).
// The C++17 special functions (std::cyl_bessel_j and others) are not in every standard
// library; Apple's libc++ has none of them.
#pragma once

#include <cmath>

namespace kairoq::detail {

// J0(x), within 2e-14: the power series below 8 (its largest term is about 100); from 8 to
// 25 the trapezoidal rule with 48 intervals on J0(x) = (1/pi) int_0^pi cos(x sin t) dt,
// whose error for this periodic integrand is 2 (J_96(x) + J_192(x) + ...), negligible;
// beyond, Hankel's asymptotic expansion (Abramowitz and Stegun 9.2.5, 9.2.9, 9.2.10),
// summed to its smallest term, about exp(-2x).
inline double bessel_j0(double x) {
    constexpr double pi = 3.14159265358979323846;
    x = std::abs(x);
    if (x < 8.0) {
        const double q = -0.25 * x * x;
        double term = 1.0, sum = 1.0;
        for (int k = 1; k < 80; ++k) {
            term *= q / (double(k) * double(k));
            sum += term;
            if (std::abs(term) < 1e-18) break;
        }
        return sum;
    }
    if (x < 25.0) {
        constexpr int n = 48;
        double sum = 1.0;  // the two ends, each cos(0) with weight one half
        for (int k = 1; k < n; ++k) sum += std::cos(x * std::sin(pi * k / n));
        return sum / n;
    }
    // J0 = sqrt(2 / (pi x)) (P cos chi - Q sin chi), chi = x - pi/4, with the terms
    // b_k = 1^2 3^2 ... (2k-1)^2 / (k! (8x)^k): P = b_0 - b_2 + b_4 ..., Q = -b_1 + b_3 ...
    double b = 1.0, p = 0.0, q = 0.0;
    for (int k = 0; k < 200; ++k) {
        const double sign = (k / 2) % 2 ? -1.0 : 1.0;
        if (k % 2 == 0) p += sign * b;
        else q -= sign * b;
        const double next = b * double(2 * k + 1) * double(2 * k + 1) / (double(k + 1) * 8.0 * x);
        if (next > b || next < 1e-18) break;
        b = next;
    }
    const double chi = x - 0.25 * pi;
    return std::sqrt(2.0 / (pi * x)) * (p * std::cos(chi) - q * std::sin(chi));
}

// I0(x) by its power series, whose terms are all positive.
inline double bessel_i0(double x) {
    const double q = 0.25 * x * x;
    double term = 1.0, sum = 1.0;
    for (int k = 1; k < 1000; ++k) {
        term *= q / (double(k) * double(k));
        sum += term;
        if (term < 1e-17 * sum) break;
    }
    return sum;
}

// K0(x), x > 0, by the trapezoidal rule with step 0.1 on K0(x) = int_0^inf exp(-x cosh t) dt,
// which converges exponentially for this integrand, analytic in a strip about the real
// axis; within 3e-15 relative for x up to 20.
inline double bessel_k0(double x) {
    constexpr double h = 0.1;
    double sum = 0.5 * std::exp(-x);
    for (int k = 1; k < 10000; ++k) {
        const double v = std::exp(-x * std::cosh(k * h));
        sum += v;
        if (v < 1e-18 * sum) break;
    }
    return sum * h;
}

}  // namespace kairoq::detail
