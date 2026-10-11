// Finite element utilities shared by the solvers (internal).
#pragma once

#include "kairoq/dc.hpp"
#include "kairoq/mesh2d.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <thread>
#include <vector>

namespace kairoq::detail {

using Clock = std::chrono::steady_clock;
inline double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

struct UnionFind {
    std::vector<std::size_t> parent;
    explicit UnionFind(std::size_t n = 0) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    std::size_t add() {
        parent.push_back(parent.size());
        return parent.size() - 1;
    }
    std::size_t find(std::size_t a) {
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        return a;
    }
    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[std::max(a, b)] = std::min(a, b);
    }
};


inline bool boxes_overlap(const Rect64& a, const Rect64& b) {
    return a.left <= b.right && b.left <= a.right && a.top <= b.bottom && b.top <= a.bottom;
}

// Hole wall as a polygon of at least 16 sides (equal area): an equipotential ring needs no
// finer outline, and every vertex refines the mesh around it.
inline Path64 hole_outline(const Drill& d, double tol) {
    const double r = d.diameter_nm / 2;
    tol = std::max(tol, r * (1.0 - std::cos(kPi / 16.0)) * 0.999);
    if (d.slot_end) return capsule(d.at, *d.slot_end, r, tol);
    return circle_path(d.at, r, tol);
}

// Port areas per copper layer.
struct PortShape {
    std::map<int, PolygonSet> by_layer;
};

inline PortShape resolve_port(const Board& b, const PortSpec& p, double tol) {
    PortShape s;
    auto add_pad = [&](const Pad& pad) {
        auto& dst = s.by_layer[pad.layer];
        dst = unite(dst, pad.shape);
    };
    for (const auto& ref : p.pads) {
        const auto dot = ref.find('.');
        if (dot == std::string::npos) throw std::runtime_error("port " + p.name + ": pad '" + ref + "' is not REF.PIN");
        const std::string comp = ref.substr(0, dot), pin = ref.substr(dot + 1);
        bool found = false;
        for (const auto& pad : b.pads)
            if (pad.component == comp && pad.pin == pin) {
                add_pad(pad);
                found = true;
            }
        if (!found) throw std::runtime_error("port " + p.name + ": no pad " + ref);
    }
    if (!p.component.empty()) {
        // With a net: the component's pads on that net; without: all its pads.
        const int net = p.net.empty() ? -1 : b.find_net(p.net);
        if (!p.net.empty() && net < 0) throw std::runtime_error(p.name + ": no net '" + p.net + "'");
        bool found = false;
        for (const auto& pad : b.pads)
            if (pad.component == p.component && (net < 0 || pad.net == net)) {
                add_pad(pad);
                found = true;
            }
        if (!found)
            throw std::runtime_error(p.name + ": component " + p.component + " has no pad" + (p.net.empty() ? "" : " on " + p.net));
    }
    if (p.at) {
        if (p.layer < 1 || p.layer > static_cast<int>(b.copper.size()))
            throw std::runtime_error("port " + p.name + ": layer " + std::to_string(p.layer) + " does not exist");
        if (p.diameter_nm <= 0) throw std::runtime_error("port " + p.name + ": contact diameter must be positive");
        auto& dst = s.by_layer[p.layer];
        dst = unite(dst, union_polygons({circle_path(*p.at, p.diameter_nm / 2, tol)}));
    }
    if (s.by_layer.empty()) throw std::runtime_error("port " + p.name + " selects nothing");
    return s;
}

// Nodes of a mesh near polygons, by a uniform grid.
class NodeGrid {
public:
    NodeGrid(const TriMesh& m, double cell) : m_(m), cell_(cell) {
        if (m.x.empty()) return;
        x0_ = *std::min_element(m.x.begin(), m.x.end());
        y0_ = *std::min_element(m.y.begin(), m.y.end());
        const double x1 = *std::max_element(m.x.begin(), m.x.end());
        const double y1 = *std::max_element(m.y.begin(), m.y.end());
        cell_ = std::max(cell_, std::max(x1 - x0_, y1 - y0_) / 2000.0);
        nx_ = static_cast<std::size_t>((x1 - x0_) / cell_) + 1;
        ny_ = static_cast<std::size_t>((y1 - y0_) / cell_) + 1;
        cells_.resize(nx_ * ny_);
        for (std::size_t i = 0; i < m.x.size(); ++i)
            cells_[cell_of(m.x[i], m.y[i])].push_back(static_cast<std::uint32_t>(i));
    }
    // Nodes inside the polygons (boundary included, with a 2 nm margin).
    std::vector<std::uint32_t> inside(const PolygonSet& polys) const {
        std::vector<std::uint32_t> out;
        if (cells_.empty() || polys.empty()) return out;
        const PolygonSet grown = union_polygons(Clipper2Lib::InflatePaths(to_paths(polys), 2.0, Clipper2Lib::JoinType::Miter,
                                                                          Clipper2Lib::EndType::Polygon));
        const Rect64 b = bounds(grown);
        const auto ix0 = index(b.left, x0_, nx_), ix1 = index(b.right, x0_, nx_);
        const auto iy0 = index(b.top, y0_, ny_), iy1 = index(b.bottom, y0_, ny_);
        for (auto iy = iy0; iy <= iy1; ++iy)
            for (auto ix = ix0; ix <= ix1; ++ix)
                for (auto n : cells_[iy * nx_ + ix])
                    if (contains(grown, Point64(std::llround(m_.x[n]), std::llround(m_.y[n])))) out.push_back(n);
        return out;
    }

private:
    std::size_t index(double v, double v0, std::size_t n) const {
        const double i = std::floor((v - v0) / cell_);
        return static_cast<std::size_t>(std::clamp(i, 0.0, static_cast<double>(n - 1)));
    }
    std::size_t cell_of(double x, double y) const { return index(y, y0_, ny_) * nx_ + index(x, x0_, nx_); }
    const TriMesh& m_;
    double cell_;
    double x0_ = 0, y0_ = 0;
    std::size_t nx_ = 0, ny_ = 0;
    std::vector<std::vector<std::uint32_t>> cells_;
};

// True on the worker threads of parallel_for: a nested parallel_for runs serially.
inline bool& in_parallel_worker() {
    thread_local bool inside = false;
    return inside;
}

// Run f(i) for i in [0, n) on all hardware threads (serially when nested).
template <typename F>
inline void parallel_for(std::size_t n, F f) {
    if (in_parallel_worker()) {
        for (std::size_t i = 0; i < n; ++i) f(i);
        return;
    }
    const unsigned threads = std::max(1u, std::min<unsigned>(std::thread::hardware_concurrency(), static_cast<unsigned>(n)));
    std::atomic<std::size_t> next{0};
    std::vector<std::exception_ptr> errors(threads);
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            in_parallel_worker() = true;
            try {
                for (std::size_t i = next++; i < n; i = next++) f(i);
            } catch (...) {
                errors[t] = std::current_exception();
            }
        });
    }
    for (auto& th : pool) th.join();
    for (auto& e : errors)
        if (e) std::rethrow_exception(e);
}

struct Triplet {
    std::size_t row, col;
    double value;
};

// Compressed sparse rows from triplets (duplicates summed).
struct Csr {
    std::size_t n = 0;
    std::vector<std::ptrdiff_t> ptr;
    std::vector<std::ptrdiff_t> col;
    std::vector<double> val;
};

inline Csr to_csr(std::size_t n, const std::vector<Triplet>& t) {
    // Counting sort by row, then each row sorted by column and duplicates summed.
    std::vector<std::size_t> start(n + 1, 0);
    for (const auto& e : t) ++start[e.row + 1];
    for (std::size_t r = 0; r < n; ++r) start[r + 1] += start[r];
    std::vector<std::pair<std::size_t, double>> entries(t.size());
    {
        std::vector<std::size_t> fill(start.begin(), start.end() - 1);
        for (const auto& e : t) entries[fill[e.row]++] = {e.col, e.value};
    }
    std::vector<std::size_t> unique(n, 0);
    parallel_for((n + 4095) / 4096, [&](std::size_t chunk) {
        for (std::size_t r = chunk * 4096; r < std::min(n, (chunk + 1) * 4096); ++r) {
            auto first = entries.begin() + static_cast<std::ptrdiff_t>(start[r]);
            auto last = entries.begin() + static_cast<std::ptrdiff_t>(start[r + 1]);
            std::sort(first, last, [](const auto& a, const auto& b) { return a.first < b.first; });
            std::size_t u = 0;
            for (auto it = first; it != last; ++it) {
                if (u > 0 && (first + static_cast<std::ptrdiff_t>(u - 1))->first == it->first)
                    (first + static_cast<std::ptrdiff_t>(u - 1))->second += it->second;
                else
                    *(first + static_cast<std::ptrdiff_t>(u++)) = *it;
            }
            unique[r] = u;
        }
    });
    Csr m;
    m.n = n;
    m.ptr.assign(n + 1, 0);
    for (std::size_t r = 0; r < n; ++r) m.ptr[r + 1] = m.ptr[r] + static_cast<std::ptrdiff_t>(unique[r]);
    m.col.resize(static_cast<std::size_t>(m.ptr[n]));
    m.val.resize(static_cast<std::size_t>(m.ptr[n]));
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t k = 0; k < unique[r]; ++k) {
            m.col[static_cast<std::size_t>(m.ptr[r]) + k] = static_cast<std::ptrdiff_t>(entries[start[r] + k].first);
            m.val[static_cast<std::size_t>(m.ptr[r]) + k] = entries[start[r] + k].second;
        }
    return m;
}

// ---------------------------------------------------------------------- point location

// Triangles of a mesh on a uniform grid: the triangle containing a point and its
// barycentric coordinates.
class TriangleLocator {
public:
    explicit TriangleLocator(const TriMesh& m);
    // Triangle containing (x, y) (nm), or the nearest one when the point lies outside
    // the mesh by less than `slack_nm`; -1 otherwise. w receives barycentric weights.
    long locate(double x, double y, std::array<double, 3>& w, double slack_nm = 0.0) const;

private:
    const TriMesh& m_;
    double x0_ = 0, y0_ = 0, cell_ = 1;
    std::size_t nx_ = 0, ny_ = 0;
    std::vector<std::vector<std::uint32_t>> cells_;
};

// ---------------------------------------------------------------------- estimation

// Gradient-recovery (Zienkiewicz-Zhu) estimate on one mesh. `value` gives the nodal
// solution, `conductance` the sheet conductance of each triangle (W/K or S); triangles
// with `skip` set are left out. Returns the squared estimate per triangle and adds to
// the totals the estimate, the energy (sum of conductance |grad u|^2 A) and the area.
// With `material` (0 or 1 per triangle), the gradient is recovered within each material
// separately: across a change of conductance the gradient is discontinuous, and
// averaging over it would report an error that refinement cannot remove.
std::vector<double> zz_estimate(const TriMesh& m, const std::vector<double>& value, const std::vector<double>& conductance,
                                const std::vector<char>& skip, double& eta2_total, double& energy, double& area_total,
                                const std::vector<std::uint8_t>* material = nullptr);

struct RefinePart {
    const TriMesh* mesh;
    const std::vector<double>* eta2;
    std::size_t owner;  // index of the seed list that receives the new vertices
};

// Marks triangles in descending order of their estimate until the predicted total meets
// (margin target)^2 energy, and writes for each owner the previous vertices plus lattice
// points in its marked triangles. Returns the number marked. A margin near 1 refines less
// per pass; whether that saves a pass or costs one depends on the problem.
std::size_t mark_and_seed(const std::vector<RefinePart>& parts, double eta2_total, double energy, double target,
                          std::size_t triangles, std::size_t max_triangles, double min_edge_nm,
                          std::vector<std::vector<std::array<double, 3>>>& seeds, bool& budget_reached, double margin = 0.8);

// ---------------------------------------------------------------------- linear algebra

struct SolveStats {
    int iterations = 0;
    double residual = 0.0;
};

// Conjugate gradients with classical algebraic multigrid, set up once for a matrix and
// used for any number of right-hand sides (x is the initial guess and the solution).
class SpdSolver {
public:
    SpdSolver(std::size_t n, const std::vector<std::ptrdiff_t>& ptr, const std::vector<std::ptrdiff_t>& col,
              const std::vector<double>& val, double tolerance, int max_iterations);
    ~SpdSolver();
    SpdSolver(const SpdSolver&) = delete;
    SpdSolver& operator=(const SpdSolver&) = delete;
    SolveStats solve(const std::vector<double>& rhs, std::vector<double>& x) const;
    std::size_t size() const { return n_; }

private:
    struct Impl;
    std::size_t n_ = 0;
    std::unique_ptr<Impl> impl_;
};

// Conjugate gradients with classical algebraic multigrid; one hierarchy for all
// right-hand sides.
SolveStats solve_spd(std::size_t n, const std::vector<std::ptrdiff_t>& ptr, const std::vector<std::ptrdiff_t>& col,
                     const std::vector<double>& val, const std::vector<std::vector<double>>& rhs,
                     std::vector<std::vector<double>>& x, double tolerance, int max_iterations);

}  // namespace kairoq::detail
