#include "fem_util.hpp"

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace kairoq::detail {

TriangleLocator::TriangleLocator(const TriMesh& m) : m_(m) {
    if (m.triangles.empty()) return;
    double x1 = -INFINITY, y1 = -INFINITY;
    x0_ = y0_ = INFINITY;
    for (std::size_t i = 0; i < m.x.size(); ++i) {
        x0_ = std::min(x0_, m.x[i]);
        y0_ = std::min(y0_, m.y[i]);
        x1 = std::max(x1, m.x[i]);
        y1 = std::max(y1, m.y[i]);
    }
    // About two triangles per cell.
    const double area = std::max((x1 - x0_) * (y1 - y0_), 1.0);
    cell_ = std::max(std::sqrt(2.0 * area / static_cast<double>(m.triangles.size())), 1.0);
    nx_ = static_cast<std::size_t>((x1 - x0_) / cell_) + 1;
    ny_ = static_cast<std::size_t>((y1 - y0_) / cell_) + 1;
    cells_.resize(nx_ * ny_);
    for (std::size_t t = 0; t < m.triangles.size(); ++t) {
        const auto& tr = m.triangles[t];
        const double ax = std::min({m.x[tr[0]], m.x[tr[1]], m.x[tr[2]]}), bx = std::max({m.x[tr[0]], m.x[tr[1]], m.x[tr[2]]});
        const double ay = std::min({m.y[tr[0]], m.y[tr[1]], m.y[tr[2]]}), by = std::max({m.y[tr[0]], m.y[tr[1]], m.y[tr[2]]});
        const auto i0 = static_cast<std::size_t>((ax - x0_) / cell_), i1 = static_cast<std::size_t>((bx - x0_) / cell_);
        const auto j0 = static_cast<std::size_t>((ay - y0_) / cell_), j1 = static_cast<std::size_t>((by - y0_) / cell_);
        for (std::size_t j = j0; j <= std::min(j1, ny_ - 1); ++j)
            for (std::size_t i = i0; i <= std::min(i1, nx_ - 1); ++i) cells_[j * nx_ + i].push_back(static_cast<std::uint32_t>(t));
    }
}

long TriangleLocator::locate(double x, double y, std::array<double, 3>& w, double slack_nm) const {
    if (cells_.empty()) return -1;
    auto bary = [&](std::uint32_t t, std::array<double, 3>& b) {
        const auto& tr = m_.triangles[t];
        const double x1 = m_.x[tr[0]], y1 = m_.y[tr[0]], x2 = m_.x[tr[1]], y2 = m_.y[tr[1]], x3 = m_.x[tr[2]], y3 = m_.y[tr[2]];
        const double d = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        b[0] = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / d;
        b[1] = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / d;
        b[2] = 1.0 - b[0] - b[1];
    };
    const long ci = static_cast<long>(std::floor((x - x0_) / cell_)), cj = static_cast<long>(std::floor((y - y0_) / cell_));
    const long reach = slack_nm > 0 ? static_cast<long>(std::ceil(slack_nm / cell_)) : 0;
    long best = -1;
    double best_out = INFINITY;
    for (long j = cj - reach; j <= cj + reach; ++j) {
        if (j < 0 || j >= static_cast<long>(ny_)) continue;
        for (long i = ci - reach; i <= ci + reach; ++i) {
            if (i < 0 || i >= static_cast<long>(nx_)) continue;
            for (auto t : cells_[static_cast<std::size_t>(j) * nx_ + static_cast<std::size_t>(i)]) {
                std::array<double, 3> b;
                bary(t, b);
                const double out = -std::min({b[0], b[1], b[2], 0.0});
                if (out <= 1e-12) {
                    w = b;
                    return t;
                }
                if (out < best_out) {
                    best_out = out;
                    best = t;
                }
            }
        }
    }
    if (best < 0 || slack_nm <= 0) return -1;
    // Nearest triangle: clamp the coordinates and check the distance.
    std::array<double, 3> b;
    bary(static_cast<std::uint32_t>(best), b);
    for (auto& v : b) v = std::max(v, 0.0);
    const double s = b[0] + b[1] + b[2];
    for (auto& v : b) v /= s;
    const auto& tr = m_.triangles[static_cast<std::size_t>(best)];
    const double px = b[0] * m_.x[tr[0]] + b[1] * m_.x[tr[1]] + b[2] * m_.x[tr[2]];
    const double py = b[0] * m_.y[tr[0]] + b[1] * m_.y[tr[1]] + b[2] * m_.y[tr[2]];
    if (std::hypot(px - x, py - y) > slack_nm) return -1;
    w = b;
    return best;
}

std::vector<double> zz_estimate(const TriMesh& m, const std::vector<double>& u, const std::vector<double>& g,
                                const std::vector<char>& skip, double& eta2_total, double& energy, double& area_total,
                                const std::vector<std::uint8_t>* material) {
    const std::size_t nt = m.triangles.size(), nv = m.x.size();
    std::vector<double> gx(nt), gy(nt), area(nt), eta2(nt, 0.0);
    // Accumulators per node and material (two classes).
    std::vector<double> nx(2 * nv, 0.0), ny(2 * nv, 0.0), nw(2 * nv, 0.0);
    auto cls = [&](std::size_t t) { return material ? static_cast<std::size_t>((*material)[t] ? 1 : 0) : 0; };
    for (std::size_t t = 0; t < nt; ++t) {
        const auto& tr = m.triangles[t];
        const double x1 = m.x[tr[0]], y1 = m.y[tr[0]], x2 = m.x[tr[1]], y2 = m.y[tr[1]], x3 = m.x[tr[2]], y3 = m.y[tr[2]];
        const double a2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        gx[t] = ((y2 - y3) * u[tr[0]] + (y3 - y1) * u[tr[1]] + (y1 - y2) * u[tr[2]]) / a2 * 1e9;  // per m
        gy[t] = ((x3 - x2) * u[tr[0]] + (x1 - x3) * u[tr[1]] + (x2 - x1) * u[tr[2]]) / a2 * 1e9;
        area[t] = 0.5 * a2 * 1e-18;  // m^2
        if (!skip.empty() && skip[t]) continue;
        const std::size_t c0 = cls(t) * nv;
        for (int c = 0; c < 3; ++c) {
            nx[c0 + tr[c]] += area[t] * gx[t];
            ny[c0 + tr[c]] += area[t] * gy[t];
            nw[c0 + tr[c]] += area[t];
        }
    }
    for (std::size_t t = 0; t < nt; ++t) {
        if (!skip.empty() && skip[t]) continue;
        const auto& tr = m.triangles[t];
        const std::size_t c0 = cls(t) * nv;
        double e = 0;
        for (int c = 0; c < 3; ++c) {
            const double dx = nx[c0 + tr[c]] / nw[c0 + tr[c]] - gx[t], dy = ny[c0 + tr[c]] / nw[c0 + tr[c]] - gy[t];
            e += dx * dx + dy * dy;
        }
        eta2[t] = g[t] * area[t] * e / 3.0;
        eta2_total += eta2[t];
        energy += g[t] * area[t] * (gx[t] * gx[t] + gy[t] * gy[t]);
        area_total += area[t];
    }
    return eta2;
}

std::size_t mark_and_seed(const std::vector<RefinePart>& parts, double eta2_total, double energy, double target,
                          std::size_t ntri, std::size_t max_triangles, double min_edge_nm,
                          std::vector<std::vector<std::array<double, 3>>>& seeds, bool& budget_reached, double margin) {
    // Splitting a triangle n ways per side divides its squared estimate by about n^2
    // (eta_T ~ h^2 for linear elements, over n^2 children); n follows from the overall
    // reduction needed.
    const double goal = std::pow(margin * target, 2) * energy;
    const int n = std::clamp(static_cast<int>(std::ceil(std::sqrt(eta2_total / goal))), 2, 3);
    std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> order;
    for (std::size_t k = 0; k < parts.size(); ++k)
        for (std::size_t t = 0; t < parts[k].eta2->size(); ++t) order.push_back({(*parts[k].eta2)[t], {k, t}});
    std::sort(order.begin(), order.end(), [](const auto& a, const auto& c) { return a.first > c.first; });
    std::vector<std::vector<char>> mark(parts.size());
    for (std::size_t k = 0; k < parts.size(); ++k) mark[k].assign(parts[k].eta2->size(), 0);
    // Each lattice point adds about three triangles once the quality refinement has
    // graded the surroundings.
    const std::size_t per_mark = 3 * static_cast<std::size_t>(n * (n + 1) / 2);
    const double h_min = std::max(min_edge_nm * 4.0, 1.0);
    auto longest = [](const TriMesh& m, const std::array<std::uint32_t, 3>& tr) {
        double l = 0;
        for (int c = 0; c < 3; ++c) l = std::max(l, std::hypot(m.x[tr[c]] - m.x[tr[(c + 1) % 3]], m.y[tr[c]] - m.y[tr[(c + 1) % 3]]));
        return l;
    };
    double predicted = eta2_total;
    std::size_t marked = 0;
    budget_reached = false;
    for (const auto& [e2, kt] : order) {
        if (predicted <= goal) break;
        if (ntri + (marked + 1) * per_mark > max_triangles) {
            budget_reached = true;
            break;
        }
        const TriMesh& m = *parts[kt.first].mesh;
        if (longest(m, m.triangles[kt.second]) / n < h_min) continue;
        mark[kt.first][kt.second] = 1;
        predicted -= e2 * (1.0 - 1.0 / (n * n));
        ++marked;
    }
    for (const auto& p : parts) seeds[p.owner].clear();
    for (std::size_t k = 0; k < parts.size(); ++k) {
        const TriMesh& m = *parts[k].mesh;
        auto& sd = seeds[parts[k].owner];
        for (std::size_t v = 0; v < m.x.size(); ++v) sd.push_back({m.x[v], m.y[v], 0.0});
        for (std::size_t t = 0; t < m.triangles.size(); ++t) {
            if (!mark[k][t]) continue;
            const auto& tr = m.triangles[t];
            const double l = longest(m, tr);
            // Barycentric lattice shifted off the edges: n(n+1)/2 interior points.
            for (int i = 0; i < n; ++i)
                for (int j = 0; i + j < n; ++j) {
                    const double l1 = (i + 1.0 / 3.0) / n, l2 = (j + 1.0 / 3.0) / n, l3 = 1.0 - l1 - l2;
                    sd.push_back({l1 * m.x[tr[0]] + l2 * m.x[tr[1]] + l3 * m.x[tr[2]],
                                  l1 * m.y[tr[0]] + l2 * m.y[tr[1]] + l3 * m.y[tr[2]], l / n});
                }
        }
    }
    return marked;
}

namespace {
using AmgBackend = amgcl::backend::builtin<double>;
// Classical (Ruge-Stueben) coarsening suits the strongly graded adaptive meshes better
// than smoothed aggregation (docs/validation/dc.md).
using AmgSolver = amgcl::make_solver<amgcl::amg<AmgBackend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::spai0>,
                                     amgcl::solver::cg<AmgBackend>>;
}  // namespace

struct SpdSolver::Impl {
    std::unique_ptr<AmgSolver> solver;
};

SpdSolver::SpdSolver(std::size_t n, const std::vector<std::ptrdiff_t>& ptr, const std::vector<std::ptrdiff_t>& col,
                     const std::vector<double>& val, double tolerance, int max_iterations)
    : n_(n), impl_(std::make_unique<Impl>()) {
#ifdef _OPENMP
    // Inside a parallel_for worker the threads are already in use.
    if (in_parallel_worker()) omp_set_num_threads(1);
#endif
    if (n == 0) return;
    AmgSolver::params prm;
    prm.solver.tol = tolerance;
    prm.solver.maxiter = static_cast<std::size_t>(max_iterations);
    const auto nn = static_cast<std::ptrdiff_t>(n);
    impl_->solver = std::make_unique<AmgSolver>(std::tie(nn, ptr, col, val), prm);
}

SpdSolver::~SpdSolver() = default;

SolveStats SpdSolver::solve(const std::vector<double>& rhs, std::vector<double>& x) const {
    SolveStats st;
    if (n_ == 0) return st;
    x.resize(n_, 0.0);
    std::size_t iters = 0;
    double error = 0;
    std::tie(iters, error) = (*impl_->solver)(rhs, x);
    st.iterations = static_cast<int>(iters);
    st.residual = error;
    return st;
}

SolveStats solve_spd(std::size_t nfree, const std::vector<std::ptrdiff_t>& ptr, const std::vector<std::ptrdiff_t>& col,
                     const std::vector<double>& val, const std::vector<std::vector<double>>& rhs,
                     std::vector<std::vector<double>>& x, double tolerance, int max_iterations) {
    SolveStats st;
    if (nfree == 0) return st;
    // Classical (Ruge-Stueben) coarsening suits the strongly graded adaptive meshes
    // better than smoothed aggregation (docs/validation/dc.md).
    using Backend = amgcl::backend::builtin<double>;
    using Solver = amgcl::make_solver<amgcl::amg<Backend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::spai0>,
                                      amgcl::solver::cg<Backend>>;
    Solver::params prm;
    prm.solver.tol = tolerance;
    prm.solver.maxiter = static_cast<std::size_t>(max_iterations);
    const auto n = static_cast<std::ptrdiff_t>(nfree);
    Solver solver(std::tie(n, ptr, col, val), prm);
    x.resize(rhs.size());
    for (std::size_t j = 0; j < rhs.size(); ++j) {
        x[j].resize(nfree, 0.0);
        std::size_t iters = 0;
        double error = 0;
        std::tie(iters, error) = solver(rhs[j], x[j]);
        st.iterations = std::max(st.iterations, static_cast<int>(iters));
        st.residual = std::max(st.residual, error);
    }
    return st;
}

}  // namespace kairoq::detail
