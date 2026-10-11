// Transmission-line cross-sections by two-dimensional finite elements (quasi-TEM).
//
// The potential of each signal conductor at 1 V, the others and the reference at 0, is
// solved with linear triangles in the cut, bounded by a grounded box. The energies give
// the Maxwell capacitance matrix: in vacuum C0, from which L = mu0 eps0 C0^-1, and with
// the dielectrics C(f) and, by first-order perturbation in the loss, G(f). The mesh is
// refined adaptively from a gradient-recovery estimate. The surface currents of the
// high-frequency limit are the charges of the vacuum solution; with the surface
// impedance of each conductor they give R(f) and the internal inductance. Method:
// docs/solvers/si.md.
#include "kairoq/cross_section.hpp"
#include "kairoq/mesh2d.hpp"

#include "fem_util.hpp"
#include "segment_grid.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace kairoq {

namespace {

using namespace detail;
using cd = std::complex<double>;

constexpr double kEps0 = 8.8541878128e-12;
constexpr double kMu0 = 1.25663706212e-6;
constexpr double kC0 = 299792458.0;

// Inverse of a small dense symmetric positive definite matrix (Gauss-Jordan).
std::vector<double> invert(std::vector<double> a, std::size_t n) {
    std::vector<double> inv(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) inv[i * n + i] = 1.0;
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::abs(a[r * n + c]) > std::abs(a[p * n + c])) p = r;
        if (std::abs(a[p * n + c]) < 1e-300) throw std::runtime_error("cross-section: singular capacitance matrix");
        for (std::size_t k = 0; k < n; ++k) {
            std::swap(a[c * n + k], a[p * n + k]);
            std::swap(inv[c * n + k], inv[p * n + k]);
        }
        const double d = a[c * n + c];
        for (std::size_t k = 0; k < n; ++k) {
            a[c * n + k] /= d;
            inv[c * n + k] /= d;
        }
        for (std::size_t r = 0; r < n; ++r) {
            if (r == c) continue;
            const double f = a[r * n + c];
            if (f == 0) continue;
            for (std::size_t k = 0; k < n; ++k) {
                a[r * n + k] -= f * a[c * n + k];
                inv[r * n + k] -= f * inv[c * n + k];
            }
        }
    }
    return inv;
}

double point_segment_distance(double px, double py, double ax, double ay, double bx, double by) {
    const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
    const double t = l2 > 0 ? std::clamp(((px - ax) * ex + (py - ay) * ey) / l2, 0.0, 1.0) : 0.0;
    return std::hypot(px - ax - t * ex, py - ay - t * ey);
}

// The mesh of the cut, all parts together.
struct CutMesh {
    std::vector<TriMesh> parts;
    std::vector<double> x, y;                        // nm
    std::vector<std::array<std::uint32_t, 3>> tri;
    std::vector<std::size_t> part_of, local_of;      // per triangle
    std::vector<int> region;                         // dielectric per triangle, -1 vacuum
    std::vector<double> weight;                      // per triangle: 1, or 2 pi r (m) for a body of revolution
    std::vector<int> owner;                          // per node: conductor index, -2 box, -1 free
    std::vector<double> boundary_length;             // per node: half the adjacent boundary edges (m)
    std::unordered_set<std::uint64_t> boundary_edges;  // (min << 32) | max
};

class Solver {
public:
    Solver(const std::vector<CsConductor>& c, const std::vector<CsDielectric>& d, const CsSettings& s)
        : cond_(c), diel_(d), s_(s) {}
    CsResult run();

private:
    const std::vector<CsConductor>& cond_;
    const std::vector<CsDielectric>& diel_;
    const CsSettings& s_;
    std::vector<std::size_t> signal_;  // input indices of the signal conductors
    std::vector<int> group_;           // per input conductor: signal number, or N for the reference
    std::size_t n_ = 0;                // signal conductors
    PolygonSet domain_;
    Rect64 box_{};
    std::vector<PolygonSet> shapes_;   // per input conductor, united
    MeshSettings ms_;
    CutMesh m_;

    void build_domain();
    void mesh(const std::vector<std::vector<std::array<double, 3>>>& seeds);
    // Solutions for each signal conductor at 1 V with permittivity eps per triangle.
    std::vector<std::vector<double>> solve(const std::vector<double>& eps) const;
    std::vector<double> capacitance(const std::vector<std::vector<double>>& phi, const std::vector<double>& eps) const;
    std::vector<double> estimate(const std::vector<std::vector<double>>& phi, const std::vector<double>& eps,
                                 double& eta2, double& energy) const;
    std::vector<double> permittivity(double f, bool imaginary) const;
    // Integral over the surface of body k of the products of the vacuum charge densities
    // of the solutions a and b (C^2/m^3 per volt squared), N x N.
    std::vector<double> surface_products(const std::vector<std::vector<double>>& phi0, int k) const;
};

void Solver::build_domain() {
    for (std::size_t i = 0; i < cond_.size(); ++i) {
        PolygonSet u = union_polygons(to_paths(cond_[i].shape));
        if (u.empty()) throw std::runtime_error("cross-section: conductor " + cond_[i].name + " has no area");
        shapes_.push_back(std::move(u));
        if (!cond_[i].reference) signal_.push_back(i);
    }
    n_ = signal_.size();
    if (n_ == 0) throw std::runtime_error("cross-section: no signal conductor");
    group_.assign(cond_.size(), static_cast<int>(n_));
    for (std::size_t k = 0; k < n_; ++k) group_[signal_[k]] = static_cast<int>(k);
    for (std::size_t i = 0; i < shapes_.size(); ++i)
        for (std::size_t j = i + 1; j < shapes_.size(); ++j)
            if (!intersection(shapes_[i], shapes_[j]).empty())
                throw std::runtime_error("cross-section: conductors " + cond_[i].name + " and " + cond_[j].name + " overlap");
    // Element size: a quarter of each signal conductor's smaller dimension (of a thin
    // strip, a fiftieth of its width) at that conductor, growing by 0.3 per unit distance.
    // Conductors of a small weight start three times coarser.
    std::vector<std::pair<Rect64, double>> sig;
    for (std::size_t k = 0; k < signal_.size(); ++k) {
        const Rect64 b = bounds(shapes_[signal_[k]]);
        const double w = double(b.right - b.left), h = double(b.bottom - b.top);
        const double coarse = k < s_.weights.size() && s_.weights[k] < 1.0 ? 3.0 : 1.0;
        sig.push_back({b, coarse * std::max(0.25 * std::max(std::min(w, h), 0.02 * std::max(w, h)), 1.0)});
    }
    Rect64 r = bounds(shapes_.front());
    auto grow = [&](const Rect64& b) {
        r.left = std::min(r.left, b.left);
        r.top = std::min(r.top, b.top);
        r.right = std::max(r.right, b.right);
        r.bottom = std::max(r.bottom, b.bottom);
    };
    for (const auto& s : shapes_) grow(bounds(s));
    const double ext = std::max(double(r.right - r.left), double(r.bottom - r.top));
    const auto m = static_cast<std::int64_t>(std::llround(s_.margin * ext));
    box_ = Rect64(r.left - m, r.top - m, r.right + m, r.bottom + m);
    if (s_.axisymmetric) {
        if (r.left < 0) throw std::runtime_error("cross-section: a body of revolution lies at non-negative radii");
        box_.left = 0;
    }
    const double boxsize = std::max(double(box_.right - box_.left), double(box_.bottom - box_.top));
    const double hmax = boxsize / 8;
    auto size = [sig, hmax](double x, double y) {
        double h = hmax;
        for (const auto& [b, h0] : sig) {
            const double dx = std::max({double(b.left) - x, 0.0, x - double(b.right)});
            const double dy = std::max({double(b.top) - y, 0.0, y - double(b.bottom)});
            h = std::min(h, h0 + 0.3 * std::hypot(dx, dy));
        }
        return h;
    };
    // Outlines subdivided to the local size, so that each segment carries its own size.
    auto subdivide = [&](const Path64& ring) {
        Path64 out;
        for (std::size_t k = 0; k < ring.size(); ++k) {
            const Point64 a = ring[k], b = ring[(k + 1) % ring.size()];
            out.push_back(a);
            std::vector<std::pair<double, double>> stack{{0.0, 1.0}};
            std::vector<double> cuts;
            while (!stack.empty()) {
                const auto [t0, t1] = stack.back();
                stack.pop_back();
                const double tm = 0.5 * (t0 + t1);
                const double mx = a.x + tm * double(b.x - a.x), my = a.y + tm * double(b.y - a.y);
                const double len = (t1 - t0) * std::hypot(double(b.x - a.x), double(b.y - a.y));
                if (len <= 1.5 * size(mx, my) || len < 4.0) continue;
                cuts.push_back(tm);
                stack.push_back({t0, tm});
                stack.push_back({tm, t1});
            }
            std::sort(cuts.begin(), cuts.end());
            for (double t : cuts) out.emplace_back(std::llround(a.x + t * double(b.x - a.x)), std::llround(a.y + t * double(b.y - a.y)));
        }
        return out;
    };
    for (auto& sh : shapes_)
        for (auto& p : sh) {
            p.outer = subdivide(p.outer);
            for (auto& h : p.holes) h = subdivide(h);
        }
    const Path64 box = subdivide({Point64(box_.left, box_.top), Point64(box_.right, box_.top), Point64(box_.right, box_.bottom),
                                  Point64(box_.left, box_.bottom)});
    PolygonSet all;
    for (const auto& s : shapes_) all = unite(all, s);
    domain_ = difference(union_polygons({box}), all);

    double smallest = INFINITY;
    for (const auto& [b, h0] : sig) smallest = std::min(smallest, h0);
    ms_.boundary_tolerance_nm = 0.0;
    ms_.min_size_nm = smallest;
    ms_.max_edge_nm = hmax;
    ms_.grading = 0.3;
    ms_.min_angle_deg = 25.0;
    ms_.min_edge_nm = std::max(smallest * 0.002, 0.2);
    ms_.constraint_size_nm = hmax;
    ms_.boundary_size = size;
}

void Solver::mesh(const std::vector<std::vector<std::array<double, 3>>>& seeds) {
    // Material interfaces, without the conductors: an interface running along a
    // conductor's face would overlap the hole's outline.
    PolygonSet constraints;
    for (const auto& d : diel_) {
        const PolygonSet part = intersection(d.shape, domain_);
        constraints.insert(constraints.end(), part.begin(), part.end());
    }
    m_ = CutMesh{};
    m_.parts.resize(domain_.size());
    parallel_for(domain_.size(), [&](std::size_t p) {
        m_.parts[p] = mesh_region(domain_[p], constraints, ms_, seeds.empty() ? std::vector<std::array<double, 3>>{} : seeds[p]);
    });
    for (std::size_t p = 0; p < m_.parts.size(); ++p) {
        const auto& t = m_.parts[p];
        const auto base = static_cast<std::uint32_t>(m_.x.size());
        m_.x.insert(m_.x.end(), t.x.begin(), t.x.end());
        m_.y.insert(m_.y.end(), t.y.begin(), t.y.end());
        for (std::size_t k = 0; k < t.triangles.size(); ++k) {
            const auto& v = t.triangles[k];
            m_.tri.push_back({base + v[0], base + v[1], base + v[2]});
            m_.part_of.push_back(p);
            m_.local_of.push_back(k);
        }
    }
    const std::size_t nt = m_.tri.size(), nv = m_.x.size();
    m_.weight.assign(nt, 1.0);
    if (s_.axisymmetric)
        for (std::size_t t = 0; t < nt; ++t) {
            const auto& v = m_.tri[t];
            m_.weight[t] = 2 * kPi * (m_.x[v[0]] + m_.x[v[1]] + m_.x[v[2]]) / 3 * 1e-9;
        }
    // Material of each triangle: the last dielectric containing its centroid.
    m_.region.assign(nt, -1);
    std::vector<Rect64> dbox;
    for (const auto& d : diel_) dbox.push_back(bounds(d.shape));
    parallel_for(nt, [&](std::size_t t) {
        const auto& v = m_.tri[t];
        const double cx = (m_.x[v[0]] + m_.x[v[1]] + m_.x[v[2]]) / 3, cy = (m_.y[v[0]] + m_.y[v[1]] + m_.y[v[2]]) / 3;
        const Point64 c(std::llround(cx), std::llround(cy));
        for (std::size_t k = diel_.size(); k-- > 0;)
            if (c.x >= dbox[k].left && c.x <= dbox[k].right && c.y >= dbox[k].top && c.y <= dbox[k].bottom && contains(diel_[k].shape, c)) {
                m_.region[t] = static_cast<int>(k);
                break;
            }
    });
    // Boundary edges (one triangle) lie on a conductor or on the box.
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edges;
    for (const auto& v : m_.tri)
        for (int k = 0; k < 3; ++k) {
            const auto a = v[k], b = v[(k + 1) % 3];
            ++edges[{std::min(a, b), std::max(a, b)}];
        }
    m_.owner.assign(nv, -1);
    m_.boundary_length.assign(nv, 0.0);
    struct Seg {
        double ax, ay, bx, by;
        int cond;
    };
    std::vector<Seg> segs;
    for (std::size_t i = 0; i < shapes_.size(); ++i)
        for (const auto& p : shapes_[i]) {
            std::vector<const Path64*> rings{&p.outer};
            for (const auto& h : p.holes) rings.push_back(&h);
            for (const auto* r : rings)
                for (std::size_t k = 0; k < r->size(); ++k) {
                    const auto& a = (*r)[k];
                    const auto& b = (*r)[(k + 1) % r->size()];
                    segs.push_back({double(a.x), double(a.y), double(b.x), double(b.y), static_cast<int>(i)});
                }
        }
    const double tol = 2.0;
    auto on_box = [&](double x, double y) {
        return (!s_.axisymmetric && std::abs(x - box_.left) < tol) || std::abs(x - box_.right) < tol ||
               std::abs(y - box_.top) < tol || std::abs(y - box_.bottom) < tol;
    };
    std::vector<std::uint32_t> bnodes;
    for (const auto& [e, count] : edges) {
        if (count != 1) continue;
        for (auto v : {e.first, e.second})
            if (m_.owner[v] == -1) {
                m_.owner[v] = -3;  // classified below
                bnodes.push_back(v);
            }
    }
    parallel_for(bnodes.size(), [&](std::size_t k) {
        const auto v = bnodes[k];
        const double x = m_.x[v], y = m_.y[v];
        if (on_box(x, y)) {
            m_.owner[v] = -2;
            return;
        }
        double best = INFINITY;
        int who = -2;
        for (const auto& s : segs) {
            const double d = point_segment_distance(x, y, s.ax, s.ay, s.bx, s.by);
            if (d < best) {
                best = d;
                who = s.cond;
            }
        }
        // On the axis of a body of revolution and off the conductors: free.
        m_.owner[v] = s_.axisymmetric && x < tol && best > tol ? -1 : who;
    });
    {
        std::vector<std::size_t> nodes(shapes_.size(), 0);
        for (std::size_t v = 0; v < nv; ++v)
            if (m_.owner[v] >= 0) ++nodes[static_cast<std::size_t>(m_.owner[v])];
        for (std::size_t i = 0; i < shapes_.size(); ++i)
            if (nodes[i] == 0)
                throw std::runtime_error("cross-section: the mesh does not reach conductor " + cond_[i].name);
    }
    for (const auto& [e, count] : edges) {
        if (count != 1) continue;
        m_.boundary_edges.insert((std::uint64_t(e.first) << 32) | e.second);
        const int a = m_.owner[e.first], b = m_.owner[e.second];
        if (a < 0 || a != b) continue;
        const double len = std::hypot(m_.x[e.first] - m_.x[e.second], m_.y[e.first] - m_.y[e.second]) * 1e-9;
        m_.boundary_length[e.first] += 0.5 * len;
        m_.boundary_length[e.second] += 0.5 * len;
    }
}

std::vector<double> Solver::permittivity(double f, bool imaginary) const {
    std::vector<double> eps(m_.tri.size(), imaginary ? 0.0 : 1.0);
    std::vector<cd> per(diel_.size());
    for (std::size_t k = 0; k < diel_.size(); ++k) per[k] = diel_[k].material.at(f);
    for (std::size_t t = 0; t < m_.tri.size(); ++t)
        if (m_.region[t] >= 0) {
            const cd e = per[static_cast<std::size_t>(m_.region[t])];
            eps[t] = imaginary ? -e.imag() : e.real();
        }
    for (std::size_t t = 0; t < eps.size(); ++t) eps[t] *= m_.weight[t];
    return eps;
}

std::vector<std::vector<double>> Solver::solve(const std::vector<double>& eps) const {
    const std::size_t nv = m_.x.size();
    std::vector<long> idx(nv, -1);
    std::size_t nf = 0;
    for (std::size_t v = 0; v < nv; ++v)
        if (m_.owner[v] == -1) idx[v] = static_cast<long>(nf++);
    std::vector<Triplet> trip;
    trip.reserve(m_.tri.size() * 9);
    // Coupling of each free node to each fixed node, per excitation.
    std::vector<std::vector<double>> rhs(n_, std::vector<double>(nf, 0.0));
    for (std::size_t t = 0; t < m_.tri.size(); ++t) {
        const auto& v = m_.tri[t];
        const double b[3] = {m_.y[v[1]] - m_.y[v[2]], m_.y[v[2]] - m_.y[v[0]], m_.y[v[0]] - m_.y[v[1]]};
        const double c[3] = {m_.x[v[2]] - m_.x[v[1]], m_.x[v[0]] - m_.x[v[2]], m_.x[v[1]] - m_.x[v[0]]};
        const double area = 0.5 * std::abs(b[0] * c[1] - b[1] * c[0]);
        if (area <= 0) continue;
        for (int i = 0; i < 3; ++i) {
            if (idx[v[i]] < 0) continue;
            for (int j = 0; j < 3; ++j) {
                const double k = eps[t] * (b[i] * b[j] + c[i] * c[j]) / (4 * area);
                if (idx[v[j]] >= 0) trip.push_back({static_cast<std::size_t>(idx[v[i]]), static_cast<std::size_t>(idx[v[j]]), k});
                else if (m_.owner[v[j]] >= 0 && group_[static_cast<std::size_t>(m_.owner[v[j]])] < static_cast<int>(n_))
                    rhs[static_cast<std::size_t>(group_[static_cast<std::size_t>(m_.owner[v[j]])])][static_cast<std::size_t>(idx[v[i]])] -= k;
            }
        }
    }
    const Csr a = to_csr(nf, trip);
    std::vector<std::vector<double>> x(n_, std::vector<double>(nf, 0.0));
    if (nf > 0) {
        const SpdSolver solver(nf, a.ptr, a.col, a.val, 1e-11, 2000);
        for (std::size_t k = 0; k < n_; ++k) solver.solve(rhs[k], x[k]);
    }
    std::vector<std::vector<double>> phi(n_, std::vector<double>(nv, 0.0));
    for (std::size_t k = 0; k < n_; ++k)
        for (std::size_t v = 0; v < nv; ++v) {
            if (idx[v] >= 0) phi[k][v] = x[k][static_cast<std::size_t>(idx[v])];
            else if (m_.owner[v] >= 0 && group_[static_cast<std::size_t>(m_.owner[v])] == static_cast<int>(k)) phi[k][v] = 1.0;
        }
    return phi;
}

// Gradient of a nodal field on a triangle (per nm).
inline std::array<double, 2> gradient(const CutMesh& m, std::size_t t, const std::vector<double>& u, double& area) {
    const auto& v = m.tri[t];
    const double b[3] = {m.y[v[1]] - m.y[v[2]], m.y[v[2]] - m.y[v[0]], m.y[v[0]] - m.y[v[1]]};
    const double c[3] = {m.x[v[2]] - m.x[v[1]], m.x[v[0]] - m.x[v[2]], m.x[v[1]] - m.x[v[0]]};
    const double a2 = b[0] * c[1] - b[1] * c[0];
    area = 0.5 * std::abs(a2);
    return {(b[0] * u[v[0]] + b[1] * u[v[1]] + b[2] * u[v[2]]) / a2, (c[0] * u[v[0]] + c[1] * u[v[1]] + c[2] * u[v[2]]) / a2};
}

std::vector<double> Solver::capacitance(const std::vector<std::vector<double>>& phi, const std::vector<double>& eps) const {
    std::vector<double> c(n_ * n_, 0.0);
    std::vector<std::array<double, 2>> g(n_);
    for (std::size_t t = 0; t < m_.tri.size(); ++t) {
        double area = 0;
        for (std::size_t k = 0; k < n_; ++k) g[k] = gradient(m_, t, phi[k], area);
        if (eps[t] == 0) continue;
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t j = i; j < n_; ++j) c[i * n_ + j] += eps[t] * (g[i][0] * g[j][0] + g[i][1] * g[j][1]) * area;
    }
    for (std::size_t i = 0; i < n_; ++i)
        for (std::size_t j = 0; j < i; ++j) c[i * n_ + j] = c[j * n_ + i];
    for (auto& v : c) v *= kEps0;
    return c;
}

// Recovered gradient, per node and material; the squared difference to each triangle's
// gradient, weighted with its permittivity, summed over the solutions relative to their
// energies.
std::vector<double> Solver::estimate(const std::vector<std::vector<double>>& phi, const std::vector<double>& eps, double& eta2,
                                     double& energy) const {
    const std::size_t nt = m_.tri.size(), nv = m_.x.size();
    std::vector<double> out(nt, 0.0);
    eta2 = 0;
    energy = 0;
    for (std::size_t k = 0; k < n_; ++k) {
        std::unordered_map<std::uint64_t, std::array<double, 3>> acc;  // node and material -> area, gx A, gy A
        acc.reserve(2 * nv);
        auto key = [&](std::uint32_t v, int region) { return (std::uint64_t(v) << 20) | std::uint64_t(region + 1); };
        std::vector<std::array<double, 2>> g(nt);
        std::vector<double> area(nt);
        double e = 0;
        for (std::size_t t = 0; t < nt; ++t) {
            g[t] = gradient(m_, t, phi[k], area[t]);
            e += eps[t] * (g[t][0] * g[t][0] + g[t][1] * g[t][1]) * area[t];
            for (auto v : m_.tri[t]) {
                auto& a = acc[key(v, m_.region[t])];
                a[0] += area[t];
                a[1] += area[t] * g[t][0];
                a[2] += area[t] * g[t][1];
            }
        }
        if (!(e > 0)) continue;
        const double wk = k < s_.weights.size() ? s_.weights[k] : 1.0;
        if (!(wk > 0)) continue;
        double sum = 0;
        for (std::size_t t = 0; t < nt; ++t) {
            double d = 0;
            for (auto v : m_.tri[t]) {
                const auto& a = acc.at(key(v, m_.region[t]));
                const double dx = a[1] / a[0] - g[t][0], dy = a[2] / a[0] - g[t][1];
                d += (dx * dx + dy * dy) / 3;
            }
            const double w = wk * eps[t] * d * area[t] / e;
            out[t] += w;
            sum += w;
        }
        eta2 += sum;
        energy += wk;
    }
    return out;
}

// Wheeler's incremental inductance rule in the domain form of the shape derivative: when
// the surface of body k moves outward with unit normal speed (velocity field V, extended
// into the field and vanishing before the other conductors), the vacuum capacitance
// changes by
//   dC_ab = eps0 * integral of [grad phi_a . grad phi_b div V - grad phi_a . (DV + DV^T) grad phi_b] dA,
// which equals the surface integral of rho_a rho_b / eps0. The area integral converges as
// the energy does, also at corners, where the charge density is singular.
std::vector<double> Solver::surface_products(const std::vector<std::vector<double>>& phi0, int k) const {
    const std::size_t nv = m_.x.size();
    auto body = [&](int owner) { return owner >= 0 ? group_[static_cast<std::size_t>(owner)] : -1; };
    // Segments of the body, in micrometres.
    std::vector<SegmentGrid::Seg> segs;
    double bx0 = INFINITY, by0 = INFINITY, bx1 = -INFINITY, by1 = -INFINITY;
    for (std::size_t i = 0; i < shapes_.size(); ++i) {
        if (group_[i] != k) continue;
        for (const auto& p : shapes_[i]) {
            std::vector<const Path64*> rings{&p.outer};
            for (const auto& h : p.holes) rings.push_back(&h);
            for (const auto* r : rings)
                for (std::size_t q = 0; q < r->size(); ++q) {
                    const auto& a = (*r)[q];
                    const auto& b = (*r)[(q + 1) % r->size()];
                    segs.push_back({a.x * 1e-3, a.y * 1e-3, b.x * 1e-3, b.y * 1e-3, true, 0.0});
                    bx0 = std::min({bx0, a.x * 1e-3, b.x * 1e-3});
                    by0 = std::min({by0, a.y * 1e-3, b.y * 1e-3});
                    bx1 = std::max({bx1, a.x * 1e-3, b.x * 1e-3});
                    by1 = std::max({by1, a.y * 1e-3, b.y * 1e-3});
                }
        }
    }
    std::vector<double> out(n_ * n_, 0.0);
    if (segs.empty()) return out;
    const auto segs_copy = segs;
    // Distance to the nearest other conductor or the box: V vanishes before it.
    double dmin = INFINITY;
    {
        const SegmentGrid grid(segs, bx0, by0, bx1, by1, std::max({(bx1 - bx0) / 64, (by1 - by0) / 64, 1e-3}));
        for (std::size_t v = 0; v < nv; ++v) {
            if (m_.owner[v] == -1 || m_.owner[v] == -3 || body(m_.owner[v]) == k) continue;
            const double x = m_.x[v] * 1e-3, y = m_.y[v] * 1e-3;
            const double dx = std::max({bx0 - x, 0.0, x - bx1}), dy = std::max({by0 - y, 0.0, y - by1});
            if (std::hypot(dx, dy) >= dmin) continue;
            for (const auto& sg : segs_copy) dmin = std::min(dmin, point_segment_distance(x, y, sg.x0, sg.y0, sg.x1, sg.y1));
        }
    }
    const double a = 0.45 * dmin;  // V = 1 within a, falling to 0 at 2a
    std::vector<std::array<double, 2>> V(nv, {0.0, 0.0});
    // On the surface: the motion that moves each adjacent face by one along its normal.
    std::vector<std::vector<std::array<double, 2>>> normals(nv);
    for (std::size_t t = 0; t < m_.tri.size(); ++t) {
        const auto& v = m_.tri[t];
        for (int e = 0; e < 3; ++e) {
            const auto p = v[e], q = v[(e + 1) % 3], r = v[(e + 2) % 3];
            if (m_.owner[p] < 0 || m_.owner[q] < 0 || body(m_.owner[p]) != k || body(m_.owner[q]) != k) continue;
            if (!m_.boundary_edges.count((std::uint64_t(std::min(p, q)) << 32) | std::max(p, q))) continue;
            // An edge of the body's surface: both ends on it and the triangle outside it.
            const double ex = m_.x[q] - m_.x[p], ey = m_.y[q] - m_.y[p], len = std::hypot(ex, ey);
            double nx = -ey / len, ny = ex / len;
            if ((m_.x[r] - m_.x[p]) * nx + (m_.y[r] - m_.y[p]) * ny < 0) {
                nx = -nx;
                ny = -ny;
            }
            normals[p].push_back({nx, ny});
            normals[q].push_back({nx, ny});
        }
    }
    for (std::size_t v = 0; v < nv; ++v) {
        if (m_.owner[v] < 0 || body(m_.owner[v]) != k || normals[v].empty()) continue;
        // Distinct normals of the faces meeting at the node.
        std::vector<std::array<double, 2>> ns;
        for (const auto& n : normals[v]) {
            bool seen = false;
            for (const auto& m : ns) seen |= n[0] * m[0] + n[1] * m[1] > 0.999999;
            if (!seen) ns.push_back(n);
        }
        if (ns.size() == 1) {
            V[v] = ns[0];
        } else {
            // Two faces: V . n1 = V . n2 = 1; nearly parallel faces: their mean normal.
            const auto& n1 = ns[0];
            const auto& n2 = ns[1];
            const double det = n1[0] * n2[1] - n1[1] * n2[0];
            if (std::abs(det) > 1e-3) V[v] = {(n2[1] - n1[1]) / det, (n1[0] - n2[0]) / det};
            else {
                const double mx = n1[0] + n2[0], my = n1[1] + n2[1], ml = std::hypot(mx, my);
                V[v] = {mx / ml, my / ml};
            }
        }
    }
    // In the field: the direction away from the nearest surface point, tapered.
    const SegmentGrid grid(segs, bx0 - 2 * a, by0 - 2 * a, bx1 + 2 * a, by1 + 2 * a, std::max(a / 2, 1e-3));
    parallel_for(nv, [&](std::size_t v) {
        if (m_.owner[v] != -1) return;
        const double x = m_.x[v] * 1e-3, y = m_.y[v] * 1e-3;
        if (x < bx0 - 2 * a || x > bx1 + 2 * a || y < by0 - 2 * a || y > by1 + 2 * a) return;
        const auto near = grid.nearest(x, y, 2 * a);
        if (!near) return;
        const auto& sg = segs_copy[near->first];
        const double px = sg.x0 + near->second * (sg.x1 - sg.x0), py = sg.y0 + near->second * (sg.y1 - sg.y0);
        const double d = std::hypot(x - px, y - py);
        if (!(d > 0) || d >= 2 * a) return;
        const double chi = d <= a ? 1.0 : (2 * a - d) / a;
        V[v] = {chi * (x - px) / d, chi * (y - py) / d};
    });
    // The area integral; coordinates in nm, so per nm of motion.
    std::vector<std::array<double, 2>> g(n_);
    for (std::size_t t = 0; t < m_.tri.size(); ++t) {
        const auto& v = m_.tri[t];
        if (V[v[0]][0] == 0 && V[v[0]][1] == 0 && V[v[1]][0] == 0 && V[v[1]][1] == 0 && V[v[2]][0] == 0 && V[v[2]][1] == 0) continue;
        std::vector<double> vx{V[v[0]][0], V[v[1]][0], V[v[2]][0]}, vy{V[v[0]][1], V[v[1]][1], V[v[2]][1]};
        // Gradients of the velocity components on the triangle.
        const double b[3] = {m_.y[v[1]] - m_.y[v[2]], m_.y[v[2]] - m_.y[v[0]], m_.y[v[0]] - m_.y[v[1]]};
        const double c[3] = {m_.x[v[2]] - m_.x[v[1]], m_.x[v[0]] - m_.x[v[2]], m_.x[v[1]] - m_.x[v[0]]};
        const double a2 = b[0] * c[1] - b[1] * c[0];
        const double area = 0.5 * std::abs(a2);
        const double dxx = (b[0] * vx[0] + b[1] * vx[1] + b[2] * vx[2]) / a2, dxy = (c[0] * vx[0] + c[1] * vx[1] + c[2] * vx[2]) / a2;
        const double dyx = (b[0] * vy[0] + b[1] * vy[1] + b[2] * vy[2]) / a2, dyy = (c[0] * vy[0] + c[1] * vy[1] + c[2] * vy[2]) / a2;
        const double div = dxx + dyy;
        double unused = 0;
        for (std::size_t q = 0; q < n_; ++q) g[q] = gradient(m_, t, phi0[q], unused);
        for (std::size_t p = 0; p < n_; ++p)
            for (std::size_t q = 0; q < n_; ++q) {
                const auto& ga = g[p];
                const auto& gb = g[q];
                // (DV + DV^T) gb
                const double sx = 2 * dxx * gb[0] + (dxy + dyx) * gb[1], sy = (dxy + dyx) * gb[0] + 2 * dyy * gb[1];
                out[p * n_ + q] += ((ga[0] * gb[0] + ga[1] * gb[1]) * div - (ga[0] * sx + ga[1] * sy)) * area;
            }
    }
    // dC per metre of motion is eps0 * out * 1e9; the surface integral is eps0 times that.
    for (auto& x : out) x *= kEps0 * kEps0 * 1e9;
    return out;
}

CsResult Solver::run() {
    const auto t0 = Clock::now();
    CsResult res;
    build_domain();
    for (auto i : signal_) res.conductors.push_back(cond_[i].name);
    double fref = 1e9;
    for (double f : s_.frequencies_hz)
        if (f > 0) {
            fref = f;
            break;
        }
    // Adaptive passes on the dielectric problem at the first frequency.
    std::vector<std::vector<std::array<double, 3>>> seeds;
    std::vector<std::vector<double>> phi;
    std::vector<double> eps;
    for (int pass = 0;; ++pass) {
        mesh(seeds);
        eps = permittivity(fref, false);
        phi = solve(eps);
        double eta2 = 0, energy = 0;
        const auto eta = estimate(phi, eps, eta2, energy);
        res.estimated_error = std::sqrt(eta2 / std::max(energy, 1.0));
        res.passes = pass + 1;
        res.triangles = m_.tri.size();
        if (res.estimated_error <= s_.target_error || pass + 1 >= s_.max_passes || m_.tri.size() >= s_.max_triangles) break;
        std::vector<std::vector<double>> eta_part(m_.parts.size());
        for (std::size_t p = 0; p < m_.parts.size(); ++p) eta_part[p].assign(m_.parts[p].triangles.size(), 0.0);
        for (std::size_t t = 0; t < m_.tri.size(); ++t) eta_part[m_.part_of[t]][m_.local_of[t]] = eta[t];
        std::vector<RefinePart> parts;
        for (std::size_t p = 0; p < m_.parts.size(); ++p) parts.push_back({&m_.parts[p], &eta_part[p], p});
        seeds.assign(m_.parts.size(), {});
        bool budget = false;
        const std::size_t marked = mark_and_seed(parts, eta2, energy, s_.target_error, m_.tri.size(), s_.max_triangles,
                                                 ms_.min_edge_nm, seeds, budget);
        if (marked == 0) break;
    }
    if (res.estimated_error > s_.target_error)
        res.warnings.push_back("estimated error " + std::to_string(res.estimated_error) + " above the target " +
                               std::to_string(s_.target_error));

    // Vacuum: capacitance, external inductance and the high-frequency currents.
    const std::vector<double>& vac = m_.weight;
    const auto phi0 = solve(vac);
    res.c0 = capacitance(phi0, vac);
    const auto c0inv = invert(res.c0, n_);
    res.l_ext.resize(n_ * n_);
    for (std::size_t k = 0; k < n_ * n_; ++k) res.l_ext[k] = kMu0 * kEps0 * c0inv[k];
    // Per body (each signal conductor, and the reference): the products of the
    // high-frequency surface currents of unit currents, integrated over its surface,
    // its area and its material.
    struct Body {
        std::vector<double> m;
        double area = 0;
        ConductorModel metal;
    };
    std::vector<Body> bodies(n_ + 1);
    for (auto& b : bodies) b.m.assign(n_ * n_, 0.0);
    bool have_reference = false;
    for (std::size_t i = 0; i < cond_.size(); ++i) {
        auto& b = bodies[static_cast<std::size_t>(group_[i])];
        b.area += area_nm2(shapes_[i]) * 1e-18;
        b.metal = cond_[i].metal;
        if (cond_[i].reference) have_reference = true;
    }
    for (std::size_t k = 0; k <= n_; ++k) {
        if (k == n_ && !have_reference) continue;
        const auto dc = surface_products(phi0, static_cast<int>(k));
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t jj = 0; jj < n_; ++jj) {
                double v = 0;
                for (std::size_t a = 0; a < n_; ++a)
                    for (std::size_t b = 0; b < n_; ++b) v += c0inv[a * n_ + i] * dc[a * n_ + b] * c0inv[b * n_ + jj];
                bodies[k].m[i * n_ + jj] = v;
            }
    }
    // Effective thickness of each body: its surface impedance tends to 1 / (sigma t_eff)
    // at low frequency, which gives its DC resistance for its own current.
    std::vector<double> teff(n_ + 1, 0.0);
    for (std::size_t k = 0; k <= n_; ++k) {
        if (k == n_ && !have_reference) continue;
        double self = 0;
        if (k < n_) self = bodies[k].m[k * n_ + k];
        else {
            for (std::size_t i = 0; i < n_; ++i) self += bodies[k].m[i * n_ + i];
            self /= double(n_);
        }
        teff[k] = self > 0 ? bodies[k].area * self : 0.0;
    }

    // The model: bodies with their effective thickness, and the field energy per region.
    CsModel& model = res.model;
    model.n = n_;
    model.l_ext = res.l_ext;
    model.c_vacuum = res.c0;
    for (std::size_t k = 0; k <= n_; ++k) {
        if (k == n_ && !have_reference) continue;
        model.bodies.push_back({bodies[k].m, teff[k], bodies[k].metal});
    }
    {
        std::vector<std::vector<double>> energy(diel_.size() + 1, std::vector<double>(n_ * n_, 0.0));  // last: vacuum
        std::vector<std::array<double, 2>> g(n_);
        for (std::size_t t = 0; t < m_.tri.size(); ++t) {
            double area = 0;
            for (std::size_t k = 0; k < n_; ++k) g[k] = gradient(m_, t, phi[k], area);
            auto& e = energy[m_.region[t] >= 0 ? static_cast<std::size_t>(m_.region[t]) : diel_.size()];
            for (std::size_t i = 0; i < n_; ++i)
                for (std::size_t j = 0; j < n_; ++j) e[i * n_ + j] += (g[i][0] * g[j][0] + g[i][1] * g[j][1]) * area * kEps0 * m_.weight[t];
        }
        for (std::size_t k = 0; k <= diel_.size(); ++k) {
            CsModel::Region r;
            r.energy = std::move(energy[k]);
            r.vacuum = k == diel_.size();
            if (!r.vacuum) r.material = diel_[k].material;
            model.regions.push_back(std::move(r));
        }
    }
    std::vector<double> fs(s_.frequencies_hz.begin(), s_.frequencies_hz.end());
    if (!s_.exact_dispersion) {
        for (double f : fs) res.frequencies.push_back(model.at(f));
    } else {
        // Dielectric solutions at anchor frequencies, interpolated in log f.
        double fmin = INFINITY, fmax = 0;
        for (double f : fs)
            if (f > 0) {
                fmin = std::min(fmin, f);
                fmax = std::max(fmax, f);
            }
        std::vector<double> anchors;
        if (fmax > fmin) {
            const int na = std::max(2, static_cast<int>(std::ceil(2 * std::log10(fmax / fmin))) + 1);
            for (int k = 0; k < na; ++k) anchors.push_back(fmin * std::pow(fmax / fmin, double(k) / (na - 1)));
        } else {
            anchors.push_back(fmax > 0 ? fmax : fref);
        }
        std::vector<std::vector<double>> ca, ga;  // per anchor: C, and G / omega
        for (double fa : anchors) {
            const auto er = permittivity(fa, false), ei = permittivity(fa, true);
            const auto p = (fa == fref) ? phi : solve(er);
            ca.push_back(capacitance(p, er));
            ga.push_back(capacitance(p, ei));
        }
        auto interp = [&](const std::vector<std::vector<double>>& a, double f) {
            if (anchors.size() == 1 || !(f > anchors.front())) return a.front();
            if (f >= anchors.back()) return a.back();
            std::size_t k = 0;
            while (k + 2 < anchors.size() && f > anchors[k + 1]) ++k;
            const double w = std::log(f / anchors[k]) / std::log(anchors[k + 1] / anchors[k]);
            std::vector<double> out(a[k].size());
            for (std::size_t i = 0; i < out.size(); ++i) out[i] = (1 - w) * a[k][i] + w * a[k + 1][i];
            return out;
        };
        for (double f : fs) {
            CsFrequency r = model.at(f);
            r.c = interp(ca, f);
            const auto gw = interp(ga, f);
            for (std::size_t k = 0; k < n_ * n_; ++k) r.g[k] = 2 * kPi * f * gw[k];
            res.frequencies.push_back(std::move(r));
        }
    }
    if (!have_reference)
        res.warnings.push_back("no reference conductor: the return is the far boundary, and its loss is not included");
    res.seconds = since(t0);
    return res;
}

}  // namespace

CsFrequency CsModel::at(double f) const {
    CsFrequency r;
    r.frequency_hz = f;
    const std::size_t nn = n * n;
    const double w = 2 * kPi * f;
    r.c.assign(nn, 0.0);
    r.g.assign(nn, 0.0);
    for (const auto& reg : regions) {
        const cd e = reg.vacuum ? cd(1.0, 0.0) : reg.material.at(f);
        for (std::size_t i = 0; i < nn; ++i) {
            r.c[i] += e.real() * reg.energy[i];
            r.g[i] += w * -e.imag() * reg.energy[i];
        }
    }
    std::vector<cd> z(nn, 0.0);
    for (const auto& b : bodies) {
        const double sigma = b.metal.conductivity_s_per_m;
        if (!(b.teff_m > 0) || !(sigma > 0)) continue;
        cd zs;
        if (!(f > 0)) {
            zs = 1.0 / (sigma * b.teff_m);
        } else {
            const double delta = 1.0 / std::sqrt(kPi * f * kMu0 * sigma);
            const cd x = cd(1, 1) * b.teff_m / delta;
            const cd ex = std::exp(-2.0 * x);
            zs = b.metal.roughness_factor(f) * cd(1, 1) / (sigma * delta) * (1.0 + ex) / (1.0 - ex);
        }
        for (std::size_t i = 0; i < nn; ++i) z[i] += b.m[i] * zs;
    }
    r.r.resize(nn);
    r.l.resize(nn);
    for (std::size_t i = 0; i < nn; ++i) {
        r.r[i] = z[i].real();
        r.l[i] = l_ext[i] + (w > 0 ? z[i].imag() / w : 0.0);
    }
    // Microstrip dispersion: C -> er C0 - (er C0 - C) / (1 + P), which raises the self
    // conductor's effective permittivity from e0 to the formula's e_f and leaves a filling
    // by one dielectric unchanged; then C / s and L_ext s with s = z_ratio sqrt(e_f / e0),
    // which gives the formula's impedance.
    if (microstrip.enabled && f > 0 && c_vacuum.size() == nn && microstrip.height_m > 0) {
        const std::size_t si = microstrip.self * n + microstrip.self;
        const double er = microstrip.substrate.at(f).real();
        const double e0 = r.c[si] / c_vacuum[si];
        if (er > e0 && e0 > 0) {
            double ef = e0, zr = 1.0;
            // Held at the formulas' range, f h <= 20 GHz mm, beyond it.
            kirschning_jansen(microstrip.width_m / microstrip.height_m, std::min(f * microstrip.height_m * 1e-6, kMicrostripFnMax), er, e0,
                              ef, zr);
            const double fill = (er - e0) / (er - ef);  // 1 + P
            const double s = zr * std::sqrt(ef / e0);
            const double c_self = r.c[si];
            for (std::size_t i = 0; i < nn; ++i) r.c[i] = (er * c_vacuum[i] - (er * c_vacuum[i] - r.c[i]) / fill) / s;
            const double gs = r.c[si] / c_self;
            for (std::size_t i = 0; i < nn; ++i) {
                r.g[i] *= gs;
                r.l[i] += l_ext[i] * (s - 1.0);
            }
        }
    }
    return r;
}

void kirschning_jansen(double u, double fn, double er, double e0, double& e_f, double& z_ratio) {
    using std::exp;
    using std::pow;
    const double p1 = 0.27488 + (0.6315 + 0.525 / pow(1 + 0.0157 * fn, 20)) * u - 0.065683 * exp(-8.7513 * u);
    const double p2 = 0.33622 * (1 - exp(-0.03442 * er));
    const double p3 = 0.0363 * exp(-4.6 * u) * (1 - exp(-pow(fn / 38.7, 4.97)));
    const double p4 = 1 + 2.751 * (1 - exp(-pow(er / 15.916, 8)));
    const double p = p1 * p2 * pow((0.1844 + p3 * p4) * fn, 1.5763);
    e_f = er - (er - e0) / (1 + p);
    const double r1 = std::min(0.03891 * pow(er, 1.4), 20.0);
    const double r2 = std::min(0.2671 * pow(u, 7), 20.0);
    const double r3 = 4.766 * exp(-3.228 * pow(u, 0.641));
    const double r4 = 0.016 + pow(0.0514 * er, 4.524);
    const double r5 = pow(fn / 28.843, 12);
    const double r6 = std::min(22.2 * pow(u, 1.92), 20.0);
    const double r7 = 1.206 - 0.3144 * exp(-r1) * (1 - exp(-r2));
    const double r8 = 1 + 1.275 * (1 - exp(-0.004625 * r3 * pow(er, 1.674) * pow(fn / 18.365, 2.745)));
    const double r9 = 5.086 * r4 * r5 / (0.3838 + 0.386 * r4) * exp(-r6) / (1 + 1.2992 * r5) * pow(er - 1, 6) /
                      (1 + 10 * pow(er - 1, 6));
    const double r10 = 0.00044 * pow(er, 2.136) + 0.0184;
    const double r11 = pow(fn / 19.47, 6) / (1 + 0.0962 * pow(fn / 19.47, 6));
    const double r12 = 1 / (1 + 0.00245 * u * u);
    const double r13 = 0.9408 * pow(e_f, r8) - 0.9603;
    const double r14 = (0.9408 - r9) * pow(e0, r8) - 0.9603;
    const double r15 = 0.707 * r10 * pow(fn / 12.3, 1.097);
    const double r16 = 1 + 0.0503 * er * er * r11 * (1 - exp(-pow(u / 15, 6)));
    const double r17 = r7 * (1 - 1.1241 * r12 / r16 * exp(-0.026 * pow(fn, 1.15656) - r15));
    z_ratio = pow(r13 / r14, r17);
}

std::complex<double> DielectricModel::at(double f) const {
    if (!wideband || loss_tangent <= 0) return cd(epsilon_r, -epsilon_r * loss_tangent);
    // Djordjevic-Sarkar: eps(w) = eps_inf + d / ln(w2 / w1) * ln((w2 + j w) / (w1 + j w)),
    // fitted to epsilon_r and the loss tangent at the given frequency.
    const double w0 = 2 * kPi * frequency_hz, w = 2 * kPi * f, w1 = omega_low, w2 = omega_high;
    const double span = std::log(w2 / w1);
    const double d = epsilon_r * loss_tangent * span / (std::atan(w0 / w1) - std::atan(w0 / w2));
    const double einf = epsilon_r - d / span * 0.5 * std::log((w2 * w2 + w0 * w0) / (w1 * w1 + w0 * w0));
    return einf + d / span * std::log(cd(w2, w) / cd(w1, w));
}

double ConductorModel::roughness_factor(double f) const {
    if (!(f > 0) || !(conductivity_s_per_m > 0)) return 1.0;
    const double delta = 1.0 / std::sqrt(kPi * f * kMu0 * conductivity_s_per_m);
    switch (roughness) {
    case Roughness::Hammerstad:
        return 1.0 + 2.0 / kPi * std::atan(1.4 * std::pow(rms_m / delta, 2));
    case Roughness::Huray:
        if (!(huray_radius_m > 0)) return 1.0;
        return 1.0 + huray_ratio / (1.0 + delta / huray_radius_m + delta * delta / (2 * huray_radius_m * huray_radius_m));
    default:
        return 1.0;
    }
}

CsResult solve_cross_section(const std::vector<CsConductor>& conductors, const std::vector<CsDielectric>& dielectrics,
                             const CsSettings& settings) {
    return Solver(conductors, dielectrics, settings).run();
}

LineQuantities line_quantities(const CsFrequency& f, std::size_t n, std::size_t i) {
    LineQuantities q;
    const double w = 2 * kPi * f.frequency_hz;
    const std::size_t k = i * n + i;
    const double r = f.r[k], l = f.l[k], g = f.g[k], c = f.c[k];
    if (!(w > 0)) {
        q.z0 = std::sqrt(l / c);
        q.epsilon_eff = kC0 * kC0 * l * c;
        q.delay_s_per_m = std::sqrt(l * c);
        return q;
    }
    const cd zs(r, w * l), ys(g, w * c);
    q.z0 = std::sqrt(zs / ys);
    const cd gamma = std::sqrt(zs * ys);
    const double beta = std::abs(gamma.imag());
    q.epsilon_eff = std::pow(beta * kC0 / w, 2);
    q.delay_s_per_m = beta / w;
    const double db = 20.0 / std::log(10.0);
    q.loss_db_per_m = db * std::abs(gamma.real());
    const double z = std::abs(q.z0);
    q.conductor_loss_db_per_m = db * r / (2 * z);
    q.dielectric_loss_db_per_m = db * g * z / 2;
    return q;
}

PairQuantities pair_quantities(const CsFrequency& f, std::size_t n, std::size_t i, std::size_t j) {
    PairQuantities p;
    // The symmetric pair: self terms averaged; the mutual capacitance of the Maxwell
    // matrix is negative.
    const double ls = 0.5 * (f.l[i * n + i] + f.l[j * n + j]), lm = f.l[i * n + j];
    const double cs = 0.5 * (f.c[i * n + i] + f.c[j * n + j]), cm = f.c[i * n + j];
    p.z_odd = std::sqrt((ls - lm) / (cs - cm));
    p.z_even = std::sqrt((ls + lm) / (cs + cm));
    p.z_differential = 2 * p.z_odd;
    p.z_common = p.z_even / 2;
    p.epsilon_eff_odd = kC0 * kC0 * (ls - lm) * (cs - cm);
    p.epsilon_eff_even = kC0 * kC0 * (ls + lm) * (cs + cm);
    return p;
}

// ------------------------------------------------------------------------ JSON

namespace {

using nlohmann::json;

PolygonSet shape_from_json(const json& j, double tol_nm) {
    auto mm = [](double v) { return v * kNmPerMm; };
    if (j.contains("rect_mm")) {
        const auto r = j["rect_mm"].get<std::vector<double>>();
        if (r.size() != 4) throw std::runtime_error("rect_mm needs [x0, y0, x1, y1]");
        const auto x0 = std::llround(mm(std::min(r[0], r[2]))), x1 = std::llround(mm(std::max(r[0], r[2])));
        const auto y0 = std::llround(mm(std::min(r[1], r[3]))), y1 = std::llround(mm(std::max(r[1], r[3])));
        return union_polygons({{Point64(x0, y0), Point64(x1, y0), Point64(x1, y1), Point64(x0, y1)}});
    }
    if (j.contains("polygon_mm")) {
        Path64 p;
        for (const auto& v : j["polygon_mm"]) p.emplace_back(std::llround(mm(v[0].get<double>())), std::llround(mm(v[1].get<double>())));
        return union_polygons({p});
    }
    if (j.contains("circle_mm")) {
        const auto c = j["circle_mm"].get<std::vector<double>>();
        return union_polygons({circle_path(Point64(std::llround(mm(c[0])), std::llround(mm(c[1]))), mm(c[2]), tol_nm)});
    }
    if (j.contains("ring_mm")) {
        const auto c = j["ring_mm"].get<std::vector<double>>();
        const Point64 o(std::llround(mm(c[0])), std::llround(mm(c[1])));
        return difference(union_polygons({circle_path(o, mm(c[3]), tol_nm)}), union_polygons({circle_path(o, mm(c[2]), tol_nm)}));
    }
    throw std::runtime_error("a shape needs rect_mm, polygon_mm, circle_mm or ring_mm");
}

json matrix(const std::vector<double>& m) { return m; }

}  // namespace

std::string solve_cross_section_json(const std::string& text) {
    const json spec = json::parse(text);
    const double tol = spec.value("arc_tolerance_um", 0.05) * 1000.0;
    std::vector<CsConductor> conductors;
    for (const auto& c : spec.at("conductors")) {
        CsConductor k;
        k.name = c.value("name", "C" + std::to_string(conductors.size() + 1));
        k.shape = shape_from_json(c, tol);
        k.reference = c.value("reference", false);
        k.metal.conductivity_s_per_m = c.value("conductivity_s_per_m", k.metal.conductivity_s_per_m);
        if (c.contains("roughness")) {
            const auto& r = c["roughness"];
            const std::string model = r.value("model", "none");
            if (model == "hammerstad") {
                k.metal.roughness = ConductorModel::Roughness::Hammerstad;
                k.metal.rms_m = r.at("rms_um").get<double>() * 1e-6;
            } else if (model == "huray") {
                k.metal.roughness = ConductorModel::Roughness::Huray;
                k.metal.huray_radius_m = r.at("radius_um").get<double>() * 1e-6;
                k.metal.huray_ratio = r.at("ratio").get<double>();
            } else if (model != "none") {
                throw std::runtime_error("roughness model is none, hammerstad or huray");
            }
        }
        conductors.push_back(std::move(k));
    }
    std::vector<CsDielectric> dielectrics;
    for (const auto& d : spec.value("dielectrics", json::array())) {
        CsDielectric k;
        k.shape = shape_from_json(d, tol);
        k.material.epsilon_r = d.value("epsilon_r", 1.0);
        k.material.loss_tangent = d.value("loss_tangent", 0.0);
        k.material.frequency_hz = d.value("frequency_hz", 1.0e9);
        k.material.wideband = d.value("wideband", true);
        dielectrics.push_back(std::move(k));
    }
    CsSettings s;
    s.target_error = spec.value("target_error", s.target_error);
    s.max_triangles = spec.value("max_triangles", s.max_triangles);
    s.max_passes = spec.value("max_passes", s.max_passes);
    s.margin = spec.value("margin", s.margin);
    s.exact_dispersion = spec.value("exact_dispersion", false);
    s.axisymmetric = spec.value("axisymmetric", false);
    if (spec.contains("frequencies_hz")) s.frequencies_hz = spec["frequencies_hz"].get<std::vector<double>>();
    CsResult r = solve_cross_section(conductors, dielectrics, s);
    const std::size_t n = r.conductors.size();
    if (spec.contains("microstrip")) {
        // Dispersion of a microstrip: {width_mm, height_mm, epsilon_r, loss_tangent?,
        // frequency_hz?, wideband?}, for the first signal conductor; the frequencies'
        // parameters are recomputed with it.
        const auto& m = spec["microstrip"];
        auto& ms = r.model.microstrip;
        ms.enabled = true;
        ms.width_m = m.at("width_mm").get<double>() * 1e-3;
        ms.height_m = m.at("height_mm").get<double>() * 1e-3;
        ms.substrate.epsilon_r = m.at("epsilon_r").get<double>();
        ms.substrate.loss_tangent = m.value("loss_tangent", 0.0);
        ms.substrate.frequency_hz = m.value("frequency_hz", 1.0e9);
        ms.substrate.wideband = m.value("wideband", true);
        for (auto& fq : r.frequencies) fq = r.model.at(fq.frequency_hz);
    }
    json out;
    out["conductors"] = r.conductors;
    out["c0"] = matrix(r.c0);
    out["l_ext"] = matrix(r.l_ext);
    out["estimated_error"] = r.estimated_error;
    out["triangles"] = r.triangles;
    out["passes"] = r.passes;
    out["seconds"] = r.seconds;
    out["warnings"] = r.warnings;
    json fr = json::array();
    for (const auto& f : r.frequencies) {
        json e{{"frequency_hz", f.frequency_hz}, {"r", f.r}, {"l", f.l}, {"g", f.g}, {"c", f.c}};
        json lines = json::array();
        for (std::size_t i = 0; i < n; ++i) {
            const auto q = line_quantities(f, n, i);
            lines.push_back({{"z0", {q.z0.real(), q.z0.imag()}},
                             {"epsilon_eff", q.epsilon_eff},
                             {"delay_ps_per_mm", q.delay_s_per_m * 1e9},
                             {"loss_db_per_m", q.loss_db_per_m},
                             {"conductor_loss_db_per_m", q.conductor_loss_db_per_m},
                             {"dielectric_loss_db_per_m", q.dielectric_loss_db_per_m}});
        }
        e["lines"] = lines;
        if (n == 2) {
            const auto p = pair_quantities(f, n, 0, 1);
            e["pair"] = {{"z_odd", p.z_odd},
                         {"z_even", p.z_even},
                         {"z_differential", p.z_differential},
                         {"z_common", p.z_common},
                         {"epsilon_eff_odd", p.epsilon_eff_odd},
                         {"epsilon_eff_even", p.epsilon_eff_even}};
        }
        fr.push_back(e);
    }
    out["frequencies"] = fr;
    return out.dump();
}

}  // namespace kairoq
