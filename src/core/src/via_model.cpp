// Physics-based via model parts (via_model.hpp).
#include "via_model.hpp"

#include "bessel.hpp"
#include "kairoq/cross_section.hpp"
#include "kairoq/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace kairoq {

namespace {

constexpr double kMu0 = 1.25663706212e-6;
constexpr double kEuler = 0.57721566490153286;

// Power series of J0, J1, Y0, Y1 (Abramowitz and Stegun 9.1.10, 9.1.13, 9.1.11).
void series(cplx z, cplx& j0, cplx& j1, cplx& y0, cplx& y1) {
    const cplx q = -0.25 * z * z;  // -(z/2)^2
    cplx t0 = 1.0, t1 = 0.5 * z;   // terms of J0 and J1
    j0 = t0;
    j1 = t1;
    double h = 0;                  // harmonic number H_k
    cplx s0 = 0.0;                 // sum for Y0
    cplx s1 = -1.0 * t1 * (-kEuler + (1.0 - kEuler));  // k = 0 term of the Y1 sum: (psi(1) + psi(2)) (z/2)
    for (int k = 1; k < 200; ++k) {
        t0 *= q / double(k * k);
        t1 *= q / double(k * (k + 1));
        h += 1.0 / k;
        j0 += t0;
        j1 += t1;
        s0 -= h * t0;  // sum_{k>=1} (-1)^{k+1} H_k (z^2/4)^k / (k!)^2
        const double psi = 2 * (h - kEuler) + 1.0 / (k + 1);  // psi(k+1) + psi(k+2)
        s1 -= psi * t1;
        if (std::abs(t0) < 1e-17 * std::abs(j0) && std::abs(t1) < 1e-17 * std::abs(j1) && k > 2) break;
    }
    const cplx lg = std::log(0.5 * z);
    y0 = 2.0 / kPi * ((lg + kEuler) * j0 + s0);
    y1 = -2.0 / (kPi * z) + 2.0 / kPi * lg * j1 + s1 / kPi;
}

}  // namespace

void hankel2(cplx z, cplx& h0, cplx& h1) {
    const cplx j(0.0, 1.0);
    if (std::abs(z) < 12.0) {
        cplx j0, j1, y0, y1;
        series(z, j0, j1, y0, y1);
        h0 = j0 - j * y0;
        h1 = j1 - j * y1;
        return;
    }
    // Asymptotic expansion (Abramowitz and Stegun 9.2.8), summed to its smallest term.
    cplx s0 = 0.0, s1 = 0.0, a0 = 1.0, a1 = 1.0, mj = 1.0;
    for (int k = 0; k < 60; ++k) {
        s0 += mj * a0;
        s1 += mj * a1;
        const double m = 2 * k + 1;
        const cplx n0 = a0 * (0.0 - m * m) / ((k + 1) * 8.0) / z, n1 = a1 * (4.0 - m * m) / ((k + 1) * 8.0) / z;
        if (std::abs(n0) > std::abs(a0) || std::abs(n1) > std::abs(a1) || std::abs(n0) + std::abs(n1) < 1e-18) break;
        a0 = n0;
        a1 = n1;
        mj *= -j;
    }
    const cplx f = std::sqrt(2.0 / (kPi * z));
    h0 = f * std::exp(-j * (z - 0.25 * kPi)) * s0;
    h1 = f * std::exp(-j * (z - 0.75 * kPi)) * s1;
}

cplx radial_self(double omega, double h, cplx k, double a) {
    cplx h0, h1;
    hankel2(k * a, h0, h1);
    return cplx(0.0, omega * kMu0 * h) * h0 / (2 * kPi * k * a * h1);
}

cplx radial_transfer(double omega, double h, cplx k, double rho) {
    cplx h0, h1;
    hankel2(k * rho, h0, h1);
    return 0.25 * omega * kMu0 * h * h0;
}

namespace {

// K0(x) I0(x), the ring-to-ring kernel of a higher mode; asymptotic beyond x = 20.
double k0i0(double x) {
    if (x <= 20.0) return detail::bessel_k0(x) * detail::bessel_i0(x);
    const double r = 1.0 / (8 * x * x);
    return 1.0 / (2 * x) * (1 + r + 27 * r * r / 2);
}

// Thin-wire partial self inductance of a straight wire, and the mutual inductance of two
// collinear segments of lengths l1, l2 with a gap s (Neumann's formula), H.
double wire_self(double l, double a) {
    return kMu0 / (2 * kPi) * (l * std::asinh(l / a) - std::sqrt(l * l + a * a) + a);
}
double collinear_mutual(double l1, double l2, double s) {
    auto f = [](double x) { return x > 0 ? x * std::log(x) : 0.0; };
    return kMu0 / (4 * kPi) * (f(l1 + l2 + s) - f(l1 + s) - f(l2 + s) + f(s));
}

}  // namespace

ViaLadder via_ladder(const ViaCut& cut, double max_segment_nm, double target_error) {
    const auto& cu = cut.copper;
    if (cut.from >= cu.size() || cut.to >= cu.size() || cut.from > cut.to) throw std::runtime_error("via: bad barrel span");
    const double a = 0.5 * cut.drill_nm;
    ViaLadder lad;
    // Nodes: the copper layers of the span, and slices of longer gaps.
    for (std::size_t l = cut.from; l <= cut.to; ++l) {
        const double z = 0.5 * (cu[l].z_bottom_nm + cu[l].z_top_nm);
        if (!lad.nodes.empty()) {
            const double z0 = lad.nodes.back().z_nm;
            const int n = static_cast<int>(std::ceil((z0 - z) / max_segment_nm - 1e-9));
            for (int k = 1; k < n; ++k) lad.nodes.push_back({z0 - (z0 - z) * k / n, -1});
        }
        lad.nodes.push_back({z, static_cast<int>(l)});
    }
    const std::size_t nn = lad.nodes.size(), ns = nn - 1;
    // Planes: copper layers with a hole or solid copper.
    std::vector<int> planes;
    for (std::size_t l = 0; l < cu.size(); ++l)
        if (cu[l].antipad_nm > 0 || cu[l].solid) planes.push_back(static_cast<int>(l));
    // Segments and the pair each lies in.
    for (std::size_t k = 0; k < ns; ++k) {
        ViaLadder::Segment sg;
        sg.upper = k;
        sg.lower = k + 1;
        const double zm = 0.5 * (lad.nodes[k].z_nm + lad.nodes[k + 1].z_nm);
        for (std::size_t p = 1; p < planes.size(); ++p) {
            const auto& up = cu[static_cast<std::size_t>(planes[p - 1])];
            const auto& lo = cu[static_cast<std::size_t>(planes[p])];
            if (zm < up.z_top_nm && zm > lo.z_bottom_nm) {
                sg.plane_upper = planes[p - 1];
                sg.plane_lower = planes[p];
                const double top = up.z_bottom_nm, bottom = lo.z_top_nm;
                const double inside = std::max(0.0, std::min(top, lad.nodes[k].z_nm) - std::max(bottom, lad.nodes[k + 1].z_nm));
                sg.weight = inside / (top - bottom);
            }
        }
        lad.segments.push_back(sg);
    }
    // Inductance of the segments: higher modes of their pair (same pair only), or the wire
    // with its image outside the pairs.
    lad.inductance.assign(ns * ns, 0.0);
    for (std::size_t i = 0; i < ns; ++i) {
        const auto& si = lad.segments[i];
        const double zi0 = lad.nodes[si.lower].z_nm, zi1 = lad.nodes[si.upper].z_nm;
        if (si.plane_upper < 0) {
            // Outside: the nearest plane's face as the image plane.
            const double l = (zi1 - zi0) * 1e-9;
            double face = NAN, best = INFINITY;
            for (int p : planes)
                for (double zf : {cu[static_cast<std::size_t>(p)].z_top_nm, cu[static_cast<std::size_t>(p)].z_bottom_nm}) {
                    const double d = std::min(std::abs(zf - zi0), std::abs(zf - zi1));
                    if (d < best) {
                        best = d;
                        face = zf;
                    }
                }
            double v = wire_self(l, a * 1e-9);
            if (!std::isnan(face)) v += collinear_mutual(l, l, 2 * best * 1e-9);
            lad.inductance[i * ns + i] = v;
            continue;
        }
        const double top = cu[static_cast<std::size_t>(si.plane_upper)].z_bottom_nm, bottom = cu[static_cast<std::size_t>(si.plane_lower)].z_top_nm;
        const double h = (top - bottom) * 1e-9;
        for (std::size_t j = 0; j < ns; ++j) {
            const auto& sj = lad.segments[j];
            if (sj.plane_upper != si.plane_upper || sj.plane_lower != si.plane_lower) continue;
            const double zj0 = lad.nodes[sj.lower].z_nm, zj1 = lad.nodes[sj.upper].z_nm;
            double v = 0;
            for (int n = 1; n <= 4000; ++n) {
                const double g = n * kPi / h;
                auto c = [&](double z0, double z1) {
                    const double u0 = (std::clamp(z0, bottom, top) - bottom) * 1e-9, u1 = (std::clamp(z1, bottom, top) - bottom) * 1e-9;
                    return (std::sin(g * u1) - std::sin(g * u0)) / g;
                };
                v += 2.0 / h * c(zi0, zi1) * c(zj0, zj1) * k0i0(g * a * 1e-9) / (2 * kPi);
            }
            lad.inductance[i * ns + j] = kMu0 * v;
        }
    }
    // Node capacitances: the barrel cut midway between nodes, with each node's pads.
    double r_hole = a, span = 0, last = NAN;
    for (const auto& c : cu) {
        r_hole = std::max({r_hole, 0.5 * c.antipad_nm, 0.5 * c.pad_nm});
        if (c.antipad_nm > 0 || c.solid) {
            if (!std::isnan(last)) span = std::max(span, last - c.z_top_nm);
            last = c.z_bottom_nm;
        }
    }
    if (span == 0) span = cu.front().z_top_nm - cu.back().z_bottom_nm;
    const double R = r_hole + 3 * span;
    auto rect = [](double x0, double y0, double x1, double y1) {
        const Path64 p{Point64(std::llround(x0), std::llround(y0)), Point64(std::llround(x1), std::llround(y0)),
                       Point64(std::llround(x1), std::llround(y1)), Point64(std::llround(x0), std::llround(y1))};
        return PolygonSet{Polygon{p, {}}};
    };
    std::vector<CsConductor> conductors;
    for (std::size_t k = 0; k < nn; ++k) {
        const double z1 = k == 0 ? cu[cut.from].z_top_nm : 0.5 * (lad.nodes[k - 1].z_nm + lad.nodes[k].z_nm);
        const double z0 = k + 1 == nn ? cu[cut.to].z_bottom_nm : 0.5 * (lad.nodes[k].z_nm + lad.nodes[k + 1].z_nm);
        CsConductor c;
        c.name = "node" + std::to_string(k);
        c.shape = rect(0, z0, a, z1);
        if (lad.nodes[k].layer >= 0) {
            const auto& l = cu[static_cast<std::size_t>(lad.nodes[k].layer)];
            if (l.pad_nm > 2 * a) c.shape = unite(c.shape, rect(0, std::max(l.z_bottom_nm, z0), 0.5 * l.pad_nm, std::min(l.z_top_nm, z1)));
        }
        conductors.push_back(c);
    }
    for (int p : planes) {
        const auto& l = cu[static_cast<std::size_t>(p)];
        CsConductor c;
        c.name = "plane" + std::to_string(p);
        c.reference = true;
        c.shape = rect(l.solid ? 0.0 : 0.5 * l.antipad_nm, l.z_bottom_nm, R, l.z_top_nm);
        conductors.push_back(c);
    }
    std::vector<CsDielectric> dielectrics;
    for (const auto& d : cut.dielectrics) {
        CsDielectric k;
        k.shape = rect(0, d.z_bottom_nm, 2 * R, d.z_top_nm);
        k.material.epsilon_r = d.epsilon_r;
        k.material.wideband = false;
        dielectrics.push_back(k);
    }
    // The pieces touch, so the field of one piece against its neighbour is singular where
    // they meet and their mutual capacitances grow without bound as the mesh refines. Each
    // node takes the charge of its piece with the whole barrel at one potential, the sum of
    // its row, which converges: a capacitance to the planes, the mutual ones dropped. The
    // passes are few, as the estimate follows the singular solutions.
    CsSettings s;
    s.axisymmetric = true;
    s.target_error = target_error;
    s.margin = 0.2;
    s.max_passes = 4;
    s.max_triangles = 100000;
    const CsResult r = solve_cross_section(conductors, dielectrics, s);
    const auto& c = r.frequencies.front().c;  // with the dielectrics (c0 is the vacuum's)
    lad.capacitance.assign(nn * nn, 0.0);
    for (std::size_t i = 0; i < nn; ++i)
        for (std::size_t j = 0; j < nn; ++j) lad.capacitance[i * nn + i] += c[i * nn + j];
    return lad;
}

}  // namespace kairoq
