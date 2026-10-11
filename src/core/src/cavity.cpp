// Plane-pair cavity (cavity.hpp).
//
// The voltage V between the planes satisfies div grad V + k^2 V = 0 on the region, with
// dV/dn = 0 at its outline and, at a port circle of radius a carrying current I into the
// pair, dV/dn = j omega mu h I / (2 pi a) (the radial line's condition). With linear
// elements, (K - k^2 M) V = j omega mu h sum_p I_p f_p, f_p the circle's load vector divided
// by its length, so that f_p^T V is the mean voltage on the circle; shorted ports hold V = 0.
// Z = j omega mu h F^T (K - k^2 M)^-1 F, expanded in the M-orthonormal modes K phi = lambda
// M phi: sum_m (f^T phi_m)(f^T phi_m)^T / (lambda_m - k^2). Modes up to (2 k_max)^2 come
// from shift-invert Lanczos; the rest enter by their static sum, F^T (K + s M)^-1 F less the
// kept modes' part at the same shift.
#include "cavity.hpp"

#include "fem_util.hpp"
#include "kairoq/mesh2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <random>
#include <stdexcept>

namespace kairoq {

namespace {

using namespace detail;

constexpr double kMu0 = 1.25663706212e-6;
constexpr double kC0 = 299792458.0;

// Eigenvalues and eigenvectors of a symmetric tridiagonal matrix by the implicit QL
// method: d the diagonal, e[i] = T(i, i + 1); z (n x n, row major) receives the
// eigenvectors as columns.
void tridiagonal_eigen(std::vector<double>& d, std::vector<double> e, std::vector<double>& z) {
    const int n = static_cast<int>(d.size());
    e.resize(static_cast<std::size_t>(n), 0.0);
    e[static_cast<std::size_t>(n - 1)] = 0.0;
    z.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) z[static_cast<std::size_t>(i) * n + i] = 1.0;
    auto Z = [&](int r, int c) -> double& { return z[static_cast<std::size_t>(r) * n + c]; };
    auto D = [&](int i) -> double& { return d[static_cast<std::size_t>(i)]; };
    auto E = [&](int i) -> double& { return e[static_cast<std::size_t>(i)]; };
    for (int l = 0; l < n; ++l) {
        int iter = 0, m;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::abs(D(m)) + std::abs(D(m + 1));
                if (std::abs(E(m)) <= 1e-15 * dd) break;
            }
            if (m != l) {
                if (++iter > 100) throw std::runtime_error("cavity: tridiagonal eigenproblem did not converge");
                double g = (D(l + 1) - D(l)) / (2.0 * E(l));
                double r = std::hypot(g, 1.0);
                g = D(m) - D(l) + E(l) / (g + std::copysign(r, g));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; --i) {
                    double f = s * E(i);
                    const double b = c * E(i);
                    r = std::hypot(f, g);
                    E(i + 1) = r;
                    if (r == 0.0) {
                        D(i + 1) -= p;
                        E(m) = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = D(i + 1) - p;
                    r = (D(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    D(i + 1) = g + p;
                    g = c * r - b;
                    for (int k = 0; k < n; ++k) {
                        f = Z(k, i + 1);
                        Z(k, i + 1) = s * Z(k, i) + c * f;
                        Z(k, i) = c * Z(k, i) - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                D(l) -= p;
                E(l) = g;
                E(m) = 0.0;
            }
        } while (m != l);
    }
}

std::vector<double> multiply(const Csr& a, const std::vector<double>& x) {
    std::vector<double> y(a.n, 0.0);
    parallel_for(a.n, [&](std::size_t r) {
        double acc = 0;
        for (auto k = a.ptr[r]; k < a.ptr[r + 1]; ++k) acc += a.val[static_cast<std::size_t>(k)] * x[static_cast<std::size_t>(a.col[static_cast<std::size_t>(k)])];
        y[r] = acc;
    });
    return y;
}

// Cholesky factorisation of a sparse symmetric positive definite matrix in envelope form,
// after reverse Cuthill-McKee ordering: exact solves, cheap for the narrow envelopes of
// two-dimensional meshes.
class EnvelopeCholesky {
public:
    explicit EnvelopeCholesky(const Csr& a) : n_(a.n) {
        // Reverse Cuthill-McKee from a node of least degree, per connected part.
        std::vector<char> seen(n_, 0);
        std::vector<std::size_t> order;
        order.reserve(n_);
        auto degree = [&](std::size_t v) { return static_cast<std::size_t>(a.ptr[v + 1] - a.ptr[v]); };
        std::vector<std::size_t> by_degree(n_);
        for (std::size_t v = 0; v < n_; ++v) by_degree[v] = v;
        std::sort(by_degree.begin(), by_degree.end(), [&](std::size_t x, std::size_t y) { return degree(x) < degree(y); });
        for (std::size_t start : by_degree) {
            if (seen[start]) continue;
            seen[start] = 1;
            std::size_t head = order.size();
            order.push_back(start);
            while (head < order.size()) {
                const std::size_t v = order[head++];
                std::vector<std::size_t> next;
                for (auto k = a.ptr[v]; k < a.ptr[v + 1]; ++k) {
                    const auto u = static_cast<std::size_t>(a.col[static_cast<std::size_t>(k)]);
                    if (!seen[u]) {
                        seen[u] = 1;
                        next.push_back(u);
                    }
                }
                std::sort(next.begin(), next.end(), [&](std::size_t x, std::size_t y) { return degree(x) < degree(y); });
                order.insert(order.end(), next.begin(), next.end());
            }
        }
        std::reverse(order.begin(), order.end());
        perm_ = order;
        inv_.assign(n_, 0);
        for (std::size_t i = 0; i < n_; ++i) inv_[perm_[i]] = i;
        // Envelope: row i holds columns first_[i] .. i.
        first_.assign(n_, 0);
        start_.assign(n_ + 1, 0);
        for (std::size_t i = 0; i < n_; ++i) {
            std::size_t f = i;
            const std::size_t v = perm_[i];
            for (auto k = a.ptr[v]; k < a.ptr[v + 1]; ++k) f = std::min(f, inv_[static_cast<std::size_t>(a.col[static_cast<std::size_t>(k)])]);
            first_[i] = f;
            start_[i + 1] = start_[i] + (i - f + 1);
        }
        l_.assign(start_[n_], 0.0);
        for (std::size_t i = 0; i < n_; ++i) {
            const std::size_t v = perm_[i];
            for (auto k = a.ptr[v]; k < a.ptr[v + 1]; ++k) {
                const std::size_t j = inv_[static_cast<std::size_t>(a.col[static_cast<std::size_t>(k)])];
                if (j <= i) at(i, j) += a.val[static_cast<std::size_t>(k)];
            }
        }
        for (std::size_t i = 0; i < n_; ++i) {
            double* li = &l_[start_[i]] - first_[i];  // li[j] = L(i, j)
            for (std::size_t j = first_[i]; j < i; ++j) {
                const double* lj = &l_[start_[j]] - first_[j];
                const std::size_t k0 = std::max(first_[i], first_[j]);
                double sum = li[j];
                for (std::size_t k = k0; k < j; ++k) sum -= li[k] * lj[k];
                li[j] = sum / lj[j];
            }
            double d = li[i];
            for (std::size_t k = first_[i]; k < i; ++k) d -= li[k] * li[k];
            if (!(d > 0)) throw std::runtime_error("cavity: the shifted matrix is not positive definite");
            li[i] = std::sqrt(d);
        }
    }
    void solve(const std::vector<double>& b, std::vector<double>& x) const {
        std::vector<double> y(n_);
        for (std::size_t i = 0; i < n_; ++i) y[i] = b[perm_[i]];
        for (std::size_t i = 0; i < n_; ++i) {  // L y = b
            const double* li = &l_[start_[i]] - first_[i];
            double sum = y[i];
            for (std::size_t k = first_[i]; k < i; ++k) sum -= li[k] * y[k];
            y[i] = sum / li[i];
        }
        for (std::size_t i = n_; i-- > 0;) {  // L^T x = y
            const double* li = &l_[start_[i]] - first_[i];
            y[i] /= li[i];
            for (std::size_t k = first_[i]; k < i; ++k) y[k] -= li[k] * y[i];
        }
        x.assign(n_, 0.0);
        for (std::size_t i = 0; i < n_; ++i) x[perm_[i]] = y[i];
    }

private:
    std::size_t n_;
    std::vector<std::size_t> perm_, inv_, first_, start_;
    std::vector<double> l_;
    double& at(std::size_t i, std::size_t j) { return l_[start_[i] + (j - first_[i])]; }
};

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

}  // namespace

Cavity::Cavity(const PolygonSet& region, const std::vector<CavityPort>& ports, double k_max) {
    if (!(k_max > 0)) throw std::runtime_error("cavity: no frequency range");
    // The region without the port circles.
    Paths64 circles;
    double rmin = INFINITY;
    for (const auto& p : ports) {
        if (!(p.radius_nm > 0)) throw std::runtime_error("cavity: a port without a radius");
        circles.push_back(circle_path(Point64(std::llround(p.x_nm), std::llround(p.y_nm)), p.radius_nm, std::min(1000.0, 0.01 * p.radius_nm)));
        rmin = std::min(rmin, p.radius_nm);
    }
    const PolygonSet domain = circles.empty() ? region : difference(region, union_polygons(circles));
    const double k_eig = 3.0 * k_max;  // modes kept: up to three times the highest wavenumber
    MeshSettings ms;
    ms.boundary_tolerance_nm = std::min(500.0, 0.005 * rmin);
    ms.max_edge_nm = 2 * kPi / (2 * k_max) / 20 * 1e9;  // linear elements: modes to twice k_max within about 1 %
    ms.min_size_nm = std::isfinite(rmin) ? 0.5 * rmin : ms.max_edge_nm / 4;
    ms.grading = 0.3;
    ms.min_angle_deg = 25.0;
    ms.outer_size_nm = std::min(ms.max_edge_nm, 20 * ms.min_size_nm);
    // Fine at the ports only: the outline and its other holes take a quarter of the largest
    // edge, which resolves the modes.
    ms.boundary_size = [&ports, &ms](double x, double y) {
        double size = ms.max_edge_nm / 4;
        for (const auto& p : ports) {
            const double d = std::hypot(x - p.x_nm, y - p.y_nm) - p.radius_nm;
            size = std::min(size, 0.5 * p.radius_nm + 0.3 * std::max(d, 0.0));
        }
        return size;
    };
    ms.min_edge_nm = std::max(0.02 * ms.min_size_nm, 1.0);
    std::vector<double> x, y;  // m
    std::vector<std::array<std::size_t, 3>> tri;
    for (const auto& poly : domain) {
        const TriMesh m = mesh_region(poly, {}, ms);
        const std::size_t base = x.size();
        for (std::size_t i = 0; i < m.x.size(); ++i) {
            x.push_back(m.x[i] * 1e-9);
            y.push_back(m.y[i] * 1e-9);
        }
        for (const auto& t : m.triangles) tri.push_back({base + t[0], base + t[1], base + t[2]});
    }
    const std::size_t nv = x.size();
    if (nv == 0) throw std::runtime_error("cavity: empty region");
    // Port circles: their boundary edges.
    std::map<std::pair<std::size_t, std::size_t>, int> edges;
    for (const auto& t : tri)
        for (int k = 0; k < 3; ++k) ++edges[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}];
    auto port_of = [&](std::size_t v) {
        for (std::size_t p = 0; p < ports.size(); ++p) {
            const double d = std::hypot(x[v] - ports[p].x_nm * 1e-9, y[v] - ports[p].y_nm * 1e-9);
            if (std::abs(d - ports[p].radius_nm * 1e-9) < 0.05 * ports[p].radius_nm * 1e-9 + 2e-6) return static_cast<int>(p);
        }
        return -1;
    };
    std::vector<int> on_port(nv, -1);
    std::vector<std::size_t> open_index(ports.size(), SIZE_MAX);
    for (std::size_t p = 0; p < ports.size(); ++p)
        if (!ports[p].shorted) open_index[p] = open_++;
    std::vector<std::vector<std::pair<std::size_t, double>>> load(open_);  // per open port: node, weight
    std::vector<double> circumference(open_, 0.0);
    for (const auto& [e, count] : edges) {
        if (count != 1) continue;
        const int pa = port_of(e.first), pb = port_of(e.second);
        if (pa < 0 || pa != pb) continue;
        on_port[e.first] = on_port[e.second] = pa;
        const std::size_t op = open_index[static_cast<std::size_t>(pa)];
        if (op == SIZE_MAX) continue;
        const double len = std::hypot(x[e.first] - x[e.second], y[e.first] - y[e.second]);
        load[op].push_back({e.first, 0.5 * len});
        load[op].push_back({e.second, 0.5 * len});
        circumference[op] += len;
    }
    for (std::size_t p = 0; p < open_; ++p)
        if (!(circumference[p] > 0)) throw std::runtime_error("cavity: a port lies outside the region or on its edge");
    // Free nodes: all but those of shorted ports.
    std::vector<long> idx(nv, -1);
    std::size_t nf = 0;
    for (std::size_t v = 0; v < nv; ++v)
        if (on_port[v] < 0 || !ports[static_cast<std::size_t>(on_port[v])].shorted) idx[v] = static_cast<long>(nf++);
    nodes_ = nf;
    std::vector<Triplet> tk, tm;
    tk.reserve(tri.size() * 9);
    tm.reserve(tri.size() * 9);
    for (const auto& t : tri) {
        const double b[3] = {y[t[1]] - y[t[2]], y[t[2]] - y[t[0]], y[t[0]] - y[t[1]]};
        const double c[3] = {x[t[2]] - x[t[1]], x[t[0]] - x[t[2]], x[t[1]] - x[t[0]]};
        const double area = 0.5 * std::abs(b[0] * c[1] - b[1] * c[0]);
        if (!(area > 0)) continue;
        area_ += area;
        for (int i = 0; i < 3; ++i) {
            if (idx[t[i]] < 0) continue;
            for (int j = 0; j < 3; ++j) {
                if (idx[t[j]] < 0) continue;
                const auto r = static_cast<std::size_t>(idx[t[i]]), q = static_cast<std::size_t>(idx[t[j]]);
                tk.push_back({r, q, (b[i] * b[j] + c[i] * c[j]) / (4 * area)});
                tm.push_back({r, q, area / 12.0 * (i == j ? 2.0 : 1.0)});
            }
        }
    }
    const Csr K = to_csr(nf, tk), M = to_csr(nf, tm);
    std::vector<std::vector<double>> F(open_, std::vector<double>(nf, 0.0));
    for (std::size_t p = 0; p < open_; ++p)
        for (const auto& [v, w] : load[p]) F[p][static_cast<std::size_t>(idx[v])] += w / circumference[p];
    // Shift: below the first nonzero mode of the region's size.
    const double s = std::pow(kPi / std::sqrt(area_), 2) * 0.1;
    std::vector<Triplet> ts;
    ts.reserve(tk.size() + tm.size());
    for (const auto& e : tk) ts.push_back(e);
    for (const auto& e : tm) ts.push_back({e.row, e.col, s * e.value});
    const Csr S = to_csr(nf, ts);
    const EnvelopeCholesky solver(S);
    // Shift-invert Lanczos in the M inner product, fully reorthogonalised.
    const double target = k_eig * k_eig;
    std::vector<std::vector<double>> Q, MQ;
    std::vector<double> alpha, beta;
    std::vector<double> q(nf);
    {
        std::mt19937_64 rng(1);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        for (auto& v : q) v = u(rng);
        const auto mq = multiply(M, q);
        const double norm = std::sqrt(dot(q, mq));
        for (auto& v : q) v /= norm;
    }
    std::vector<double> ritz, vectors;
    const std::size_t max_steps = nf;
    std::size_t steps = 0;
    for (std::size_t j = 0; j < max_steps; ++j) {
        ++steps;
        Q.push_back(q);
        MQ.push_back(multiply(M, q));
        std::vector<double> w;
        solver.solve(MQ.back(), w);
        const double a = dot(MQ.back(), w);
        alpha.push_back(a);
        for (int pass = 0; pass < 2; ++pass)
            for (std::size_t i = 0; i <= j; ++i) {
                const double cpl = dot(MQ[i], w);
                for (std::size_t r = 0; r < nf; ++r) w[r] -= cpl * Q[i][r];
            }
        const double b = std::sqrt(std::max(dot(w, multiply(M, w)), 0.0));
        const bool exhausted = !(b > 1e-14 * std::abs(a)) || j + 1 == max_steps;
        if ((j + 1) % 10 == 0 || exhausted) {
            ritz = alpha;
            std::vector<double> e(beta.begin(), beta.end());
            tridiagonal_eigen(ritz, e, vectors);
            const std::size_t m = ritz.size();
            // Converged when every Ritz value at or below the target has a small residual
            // and three lie beyond it.
            // Residuals against the largest Ritz value (1 / s, the static mode's scale): the
            // eigenvalue's error is about the residual squared over the gap.
            std::size_t beyond = 0;
            bool converged = true;
            const double scale = *std::max_element(ritz.begin(), ritz.end());
            for (std::size_t k = 0; k < m; ++k) {
                const double theta = ritz[k];
                if (!(theta > 0)) continue;
                const double lambda = 1.0 / theta - s;
                const double residual = std::abs(b * vectors[(m - 1) * m + k]);
                if (lambda <= target) {
                    if (residual > 1e-7 * scale) converged = false;
                } else if (residual <= 1e-4 * scale) {
                    ++beyond;
                }
            }
            if (exhausted || (converged && beyond >= 3)) break;
        }
        beta.push_back(b);
        for (std::size_t r = 0; r < nf; ++r) q[r] = w[r] / b;
    }
    steps_ = steps;
    // Kept modes, ascending, with their means over the ports.
    const std::size_t m = ritz.size();
    std::vector<std::size_t> order;
    for (std::size_t k = 0; k < m; ++k)
        if (ritz[k] > 0 && 1.0 / ritz[k] - s <= target) order.push_back(k);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return ritz[a] > ritz[b]; });
    for (auto k : order) {
        lambda_.push_back(std::max(1.0 / ritz[k] - s, 0.0));
        std::vector<double> phi(nf, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            const double c = vectors[i * m + k];
            for (std::size_t r = 0; r < nf; ++r) phi[r] += c * Q[i][r];
        }
        for (std::size_t p = 0; p < open_; ++p) mean_.push_back(dot(F[p], phi));
    }
    // Static remainder of the modes not kept.
    remainder_.assign(open_ * open_, 0.0);
    for (std::size_t p = 0; p < open_; ++p) {
        std::vector<double> xs;
        solver.solve(F[p], xs);
        for (std::size_t r = 0; r < open_; ++r) {
            double v = dot(F[r], xs);
            for (std::size_t k = 0; k < lambda_.size(); ++k) v -= mean_[k * open_ + p] * mean_[k * open_ + r] / (lambda_[k] + s);
            remainder_[r * open_ + p] = v;
        }
    }
}

std::vector<std::complex<double>> Cavity::impedance(double omega, double h_m, std::complex<double> eps_r) const {
    return impedance(omega, h_m, eps_r, 0.0);
}

std::vector<std::complex<double>> Cavity::impedance(double omega, double h_m, std::complex<double> eps_r,
                                                    std::complex<double> z_surface) const {
    // Telegrapher's form of the pair: series z = j omega mu h + 2 z_s, shunt y = j omega eps / h
    // per square; k^2 = -z y, and z in place of j omega mu h.
    using cd = std::complex<double>;
    constexpr double kEps0 = 1.0 / (kMu0 * kC0 * kC0);
    const cd jwmh = cd(0.0, omega * kMu0 * h_m) + 2.0 * z_surface;
    const cd k2 = -jwmh * cd(0.0, omega * kEps0 / h_m) * eps_r;
    std::vector<cd> z(open_ * open_, 0.0);
    for (std::size_t r = 0; r < open_; ++r)
        for (std::size_t p = 0; p < open_; ++p) {
            cd v = remainder_[r * open_ + p];
            for (std::size_t k = 0; k < lambda_.size(); ++k) v += mean_[k * open_ + r] * mean_[k * open_ + p] / (lambda_[k] - k2);
            z[r * open_ + p] = jwmh * v;
        }
    return z;
}

}  // namespace kairoq
