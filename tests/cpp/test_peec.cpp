#include "peec.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <random>

using namespace kairoq::detail;

namespace {

// Brute force: the 7-point rule on 4^depth sub-triangles.
void brute(const std::array<Vec2, 3>& v, Vec2 p, double d, int depth, double& s, Vec2& vec) {
    if (depth == 0) {
        const auto& r = triangle_rule_7();
        const double area = 0.5 * std::abs((v[1].x - v[0].x) * (v[2].y - v[0].y) - (v[2].x - v[0].x) * (v[1].y - v[0].y));
        for (std::size_t q = 0; q < r.bary.size(); ++q) {
            const double x = r.bary[q][0] * v[0].x + r.bary[q][1] * v[1].x + r.bary[q][2] * v[2].x;
            const double y = r.bary[q][0] * v[0].y + r.bary[q][1] * v[1].y + r.bary[q][2] * v[2].y;
            const double R = std::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + d * d);
            s += area * r.weight[q] / R;
            vec.x += area * r.weight[q] * (x - p.x) / R;
            vec.y += area * r.weight[q] * (y - p.y) / R;
        }
        return;
    }
    const Vec2 a{0.5 * (v[0].x + v[1].x), 0.5 * (v[0].y + v[1].y)}, b{0.5 * (v[1].x + v[2].x), 0.5 * (v[1].y + v[2].y)},
        c{0.5 * (v[2].x + v[0].x), 0.5 * (v[2].y + v[0].y)};
    brute({v[0], a, c}, p, d, depth - 1, s, vec);
    brute({a, v[1], b}, p, d, depth - 1, s, vec);
    brute({c, b, v[2]}, p, d, depth - 1, s, vec);
    brute({a, b, c}, p, d, depth - 1, s, vec);
}

}  // namespace

TEST_CASE("triangle potentials against brute-force quadrature") {
    const std::array<Vec2, 3> tri{Vec2{0.0, 0.0}, Vec2{2e-3, 0.2e-3}, Vec2{0.6e-3, 1.5e-3}};
    struct Case {
        Vec2 p;
        double d;
    };
    for (const Case c : {Case{{0.9e-3, 0.5e-3}, 0.3e-3}, Case{{3e-3, -1e-3}, 0.0}, Case{{-1e-3, 2e-3}, 0.7e-3},
                         Case{{0.9e-3, 0.5e-3}, 0.02e-3}, Case{{1.0e-3, 0.1e-3}, 0.0}}) {
        const auto t = triangle_potential(tri, c.p, c.d);
        double s = 0;
        Vec2 v;
        brute(tri, c.p, c.d, 6, s, v);
        CHECK(t.scalar == doctest::Approx(s).epsilon(2e-4));
        CHECK(t.vector.x == doctest::Approx(v.x).epsilon(2e-4));
        CHECK(t.vector.y == doctest::Approx(v.y).epsilon(2e-4));
    }
    // Far away the triangle is a point: area / R.
    const auto far = triangle_potential(tri, Vec2{0.0, 0.0}, 1.0);
    CHECK(far.scalar == doctest::Approx(1.45e-6).epsilon(1e-5));
}

TEST_CASE("partial inductances of straight conductors") {
    const double l = 10e-3, a = 0.1e-3, d = 2e-3;
    // Tube of radius a: mu0/2pi (l asinh(l/a) - sqrt(l^2 + a^2) + a).
    CHECK(parallel_filaments(0, l, 0, l, a) ==
          doctest::Approx(2e-7 * (l * std::asinh(l / a) - std::sqrt(l * l + a * a) + a)).epsilon(1e-12));
    // Parallel offset segments: closed form against quadrature.
    const double m = parallel_filaments(0, l, 3e-3, 3e-3 + 0.5 * l, d);
    const double q = segment_mutual({0, 0, 0}, {l, 0, 0}, {3e-3, d, 0}, {3e-3 + 0.5 * l, d, 0});
    CHECK(q == doctest::Approx(m).epsilon(1e-6));
    // Antiparallel segments change sign; orthogonal ones do not couple.
    CHECK(segment_mutual({0, 0, 0}, {l, 0, 0}, {l, d, 0}, {0, d, 0}) == doctest::Approx(-parallel_filaments(0, l, 0, l, d)).epsilon(1e-6));
    CHECK(segment_mutual({0, 0, 0}, {l, 0, 0}, {0, d, 0}, {0, d + l, 0}) == 0.0);
}

TEST_CASE("dense Cholesky solves") {
    const std::size_t n = 300, m = 7;
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> u(-1, 1);
    std::vector<double> g(n * n), a(n * n, 0.0);
    for (auto& x : g) x = u(rng);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < n; ++k) a[i * n + j] += g[i * n + k] * g[j * n + k];
            if (i == j) a[i * n + j] += n;
        }
    std::vector<double> x(n * m), b(n * m, 0.0);
    for (auto& v : x) v = u(rng);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < m; ++j)
            for (std::size_t k = 0; k < n; ++k) b[i * m + j] += a[i * n + k] * x[k * m + j];
    std::vector<double> c(a);
    REQUIRE(cholesky(c, n));
    forward_solve(c, n, b, m);
    backward_solve(c, n, b, m);
    double err = 0;
    for (std::size_t i = 0; i < n * m; ++i) err = std::max(err, std::abs(b[i] - x[i]));
    CHECK(err < 1e-10);
}
