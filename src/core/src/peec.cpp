#include "peec.hpp"

#include "fem_util.hpp"

#include <algorithm>
#include <cmath>

namespace kairoq::detail {

TrianglePotential triangle_potential(const std::array<Vec2, 3>& v, Vec2 rho0, double d) {
    TrianglePotential out;
    const double ad = std::abs(d);
    for (int i = 0; i < 3; ++i) {
        const Vec2 p = v[i], q = v[(i + 1) % 3];
        const double ex = q.x - p.x, ey = q.y - p.y;
        const double len = std::hypot(ex, ey);
        if (len <= 0) continue;
        const double lx = ex / len, ly = ey / len;  // tangent
        const double ux = ly, uy = -lx;             // outward normal of a counter-clockwise triangle
        const double px = p.x - rho0.x, py = p.y - rho0.y;
        const double qx = q.x - rho0.x, qy = q.y - rho0.y;
        const double sm = px * lx + py * ly, sp = qx * lx + qy * ly;
        const double t0 = px * ux + py * uy;
        const double r0sq = t0 * t0 + d * d;
        const double rm = std::sqrt(r0sq + sm * sm), rp = std::sqrt(r0sq + sp * sp);
        double f2 = 0;
        const double scale = len * len;
        if (r0sq > 1e-24 * scale) {
            if (sm > 0) f2 = std::log((rp + sp) / (rm + sm));
            else if (sp < 0) f2 = std::log((rm - sm) / (rp - sp));
            else f2 = std::log((rp + sp) * (rm - sm) / r0sq);
        }
        double beta = 0;
        if (std::abs(t0) > 0 && ad > 0)
            beta = std::atan(t0 * sp / (r0sq + ad * rp)) - std::atan(t0 * sm / (r0sq + ad * rm));
        out.scalar += t0 * f2 - ad * beta;
        const double g = 0.5 * (r0sq * f2 + sp * rp - sm * rm);
        out.vector.x += ux * g;
        out.vector.y += uy * g;
    }
    return out;
}

const TriangleRule& triangle_rule_7() {
    static const TriangleRule rule = [] {
        TriangleRule r;
        auto add3 = [&](double a, double b, double w) {
            r.bary.push_back({a, b, b});
            r.bary.push_back({b, a, b});
            r.bary.push_back({b, b, a});
            for (int k = 0; k < 3; ++k) r.weight.push_back(w);
        };
        r.bary.push_back({1.0 / 3, 1.0 / 3, 1.0 / 3});
        r.weight.push_back(0.225);
        add3(0.059715871789770, 0.470142064105115, 0.132394152788506);
        add3(0.797426985353087, 0.101286507323456, 0.125939180544827);
        return r;
    }();
    return rule;
}

const TriangleRule& triangle_rule_28() {
    static const TriangleRule rule = [] {
        // Sub-triangles of the midpoint subdivision, in barycentric coordinates.
        const std::array<std::array<std::array<double, 3>, 3>, 4> sub{{
            {{{1, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}}},
            {{{0, 1, 0}, {0, 0.5, 0.5}, {0.5, 0.5, 0}}},
            {{{0, 0, 1}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}},
            {{{0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}}},
        }};
        TriangleRule r;
        const TriangleRule& base = triangle_rule_7();
        for (const auto& s : sub)
            for (std::size_t q = 0; q < base.bary.size(); ++q) {
                std::array<double, 3> b{0, 0, 0};
                for (int k = 0; k < 3; ++k)
                    for (int c = 0; c < 3; ++c) b[c] += base.bary[q][k] * s[k][c];
                r.bary.push_back(b);
                r.weight.push_back(base.weight[q] / 4.0);
            }
        return r;
    }();
    return rule;
}

double parallel_filaments(double a1, double b1, double a2, double b2, double rho) {
    auto G = [rho](double u) { return u * std::asinh(u / rho) - std::sqrt(u * u + rho * rho); };
    return kMu0Over4Pi * (G(b1 - a2) - G(a1 - a2) - G(b1 - b2) + G(a1 - b2));
}

namespace {

double point_segment_distance(const std::array<double, 3>& p, const std::array<double, 3>& a, const std::array<double, 3>& b) {
    std::array<double, 3> ab{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, ap{p[0] - a[0], p[1] - a[1], p[2] - a[2]};
    const double l2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    double t = l2 > 0 ? (ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / l2 : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    const double dx = ap[0] - t * ab[0], dy = ap[1] - t * ab[1], dz = ap[2] - t * ab[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

}  // namespace

double segment_mutual(const std::array<double, 3>& p0, const std::array<double, 3>& p1, const std::array<double, 3>& q0,
                      const std::array<double, 3>& q1) {
    static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
    static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    const std::array<double, 3> u{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]}, v{q1[0] - q0[0], q1[1] - q0[1], q1[2] - q0[2]};
    const double lu = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]), lv = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (lu <= 0 || lv <= 0) return 0.0;
    const double dot = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (lu * lv);
    if (std::abs(dot) < 1e-12) return 0.0;
    const double dist = std::max(1e-12, std::min({point_segment_distance(p0, q0, q1), point_segment_distance(p1, q0, q1),
                                                  point_segment_distance(q0, p0, p1), point_segment_distance(q1, p0, p1)}));
    const int nu = std::clamp(static_cast<int>(std::ceil(2.0 * lu / dist)), 1, 64);
    const int nv = std::clamp(static_cast<int>(std::ceil(2.0 * lv / dist)), 1, 64);
    double sum = 0;
    for (int a = 0; a < nu; ++a)
        for (int i = 0; i < 8; ++i) {
            const double s = (a + 0.5 * (gx[i] + 1)) / nu;
            const double wi = gw[i] * 0.5 / nu;
            const double x = p0[0] + s * u[0], y = p0[1] + s * u[1], z = p0[2] + s * u[2];
            for (int b = 0; b < nv; ++b)
                for (int j = 0; j < 8; ++j) {
                    const double t = (b + 0.5 * (gx[j] + 1)) / nv;
                    const double wj = gw[j] * 0.5 / nv;
                    const double dx = x - (q0[0] + t * v[0]), dy = y - (q0[1] + t * v[1]), dz = z - (q0[2] + t * v[2]);
                    sum += wi * wj / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
        }
    return kMu0Over4Pi * dot * lu * lv * sum;
}

bool cholesky(std::vector<double>& a, std::size_t n) {
    constexpr std::size_t nb = 96;
    for (std::size_t k0 = 0; k0 < n; k0 += nb) {
        const std::size_t k1 = std::min(n, k0 + nb);
        // Diagonal block.
        for (std::size_t j = k0; j < k1; ++j) {
            double s = a[j * n + j];
            for (std::size_t k = k0; k < j; ++k) s -= a[j * n + k] * a[j * n + k];
            if (!(s > 0)) return false;
            const double d = std::sqrt(s);
            a[j * n + j] = d;
            for (std::size_t i = j + 1; i < k1; ++i) {
                double t = a[i * n + j];
                for (std::size_t k = k0; k < j; ++k) t -= a[i * n + k] * a[j * n + k];
                a[i * n + j] = t / d;
            }
        }
        // Panel below the diagonal block.
        parallel_for(n > k1 ? (n - k1 + 63) / 64 : 0, [&](std::size_t chunk) {
            const std::size_t i0 = k1 + chunk * 64, i1 = std::min(n, i0 + 64);
            for (std::size_t i = i0; i < i1; ++i)
                for (std::size_t j = k0; j < k1; ++j) {
                    double t = a[i * n + j];
                    for (std::size_t k = k0; k < j; ++k) t -= a[i * n + k] * a[j * n + k];
                    a[i * n + j] = t / a[j * n + j];
                }
        });
        // Trailing update of the lower triangle.
        parallel_for(n > k1 ? (n - k1 + 31) / 32 : 0, [&](std::size_t chunk) {
            const std::size_t i0 = k1 + chunk * 32, i1 = std::min(n, i0 + 32);
            for (std::size_t i = i0; i < i1; ++i) {
                const double* ri = &a[i * n + k0];
                for (std::size_t j = k1; j <= i; ++j) {
                    const double* rj = &a[j * n + k0];
                    double t = 0;
                    for (std::size_t k = 0; k < k1 - k0; ++k) t += ri[k] * rj[k];
                    a[i * n + j] -= t;
                }
            }
        });
    }
    return true;
}

void forward_solve(const std::vector<double>& c, std::size_t n, std::vector<double>& b, std::size_t m) {
    const std::size_t chunk = 32;
    parallel_for((m + chunk - 1) / chunk, [&](std::size_t q) {
        const std::size_t j0 = q * chunk, j1 = std::min(m, j0 + chunk);
        for (std::size_t i = 0; i < n; ++i) {
            double* bi = &b[i * m];
            for (std::size_t k = 0; k < i; ++k) {
                const double cik = c[i * n + k];
                if (cik == 0) continue;
                const double* bk = &b[k * m];
                for (std::size_t j = j0; j < j1; ++j) bi[j] -= cik * bk[j];
            }
            const double d = c[i * n + i];
            for (std::size_t j = j0; j < j1; ++j) bi[j] /= d;
        }
    });
}

void backward_solve(const std::vector<double>& c, std::size_t n, std::vector<double>& b, std::size_t m) {
    const std::size_t chunk = 32;
    parallel_for((m + chunk - 1) / chunk, [&](std::size_t q) {
        const std::size_t j0 = q * chunk, j1 = std::min(m, j0 + chunk);
        for (std::size_t ii = n; ii-- > 0;) {
            double* bi = &b[ii * m];
            const double d = c[ii * n + ii];
            for (std::size_t j = j0; j < j1; ++j) bi[j] /= d;
            for (std::size_t k = 0; k < ii; ++k) {
                const double cik = c[ii * n + k];
                if (cik == 0) continue;
                double* bk = &b[k * m];
                for (std::size_t j = j0; j < j1; ++j) bk[j] -= cik * bi[j];
            }
        }
    });
}

}  // namespace kairoq::detail
