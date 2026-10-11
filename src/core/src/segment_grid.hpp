// Line segments on a uniform grid (internal): distances and point-in-region tests.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace kairoq::detail {

// Line segments on a uniform grid: distance queries and point-in-region tests.
class SegmentGrid {
public:
    struct Seg {
        double x0, y0, x1, y1;
        bool boundary;  // part of the region outline (used for inside tests)
        double h;       // element size wanted at the segment
    };

    SegmentGrid(std::vector<Seg> segs, double x0, double y0, double x1, double y1, double cell)
        : segs_(std::move(segs)), x0_(x0), y0_(y0) {
        cell_ = std::max(cell, std::max(x1 - x0, y1 - y0) / 1024.0);
        nx_ = static_cast<long>((x1 - x0) / cell_) + 1;
        ny_ = static_cast<long>((y1 - y0) / cell_) + 1;
        cells_.resize(static_cast<std::size_t>(nx_ * ny_));
        for (std::size_t i = 0; i < segs_.size(); ++i) {
            const auto& g = segs_[i];
            const long ax = ix(std::min(g.x0, g.x1)), bx = ix(std::max(g.x0, g.x1));
            const long ay = iy(std::min(g.y0, g.y1)), by = iy(std::max(g.y0, g.y1));
            for (long y = ay; y <= by; ++y)
                for (long x = ax; x <= bx; ++x) cells_[static_cast<std::size_t>(y * nx_ + x)].push_back(static_cast<std::uint32_t>(i));
        }
        // A cell that no boundary segment touches lies wholly inside or outside the region:
        // its state, from the crossings right of it along its middle, ends an inside test.
        state_.assign(cells_.size(), kMixed);
        for (long cy = 0; cy < ny_; ++cy) {
            const double y = y0_ + (static_cast<double>(cy) + 0.5) * cell_;
            bool in = false;
            for (long cx = nx_ - 1; cx >= 0; --cx) {
                const auto c = static_cast<std::size_t>(cy * nx_ + cx);
                bool mixed = false;
                for (auto i : cells_[c]) {
                    if (!segs_[i].boundary) continue;
                    mixed = true;
                    if (crosses(segs_[i], y) && ix(crossing(segs_[i], y)) == cx) in = !in;
                }
                if (!mixed) state_[c] = in ? kInside : kOutside;
            }
        }
    }

    // True if a segment lies within distance r of (x, y).
    bool any_within(double x, double y, double r) const {
        const double r2 = r * r;
        for (long cy = iy(y - r); cy <= iy(y + r); ++cy)
            for (long cx = ix(x - r); cx <= ix(x + r); ++cx)
                for (auto i : cells_[static_cast<std::size_t>(cy * nx_ + cx)])
                    if (dist2(segs_[i], x, y) <= r2) return true;
        return false;
    }

    // Nearest segment within distance r of (x, y): its index and the parameter of the
    // closest point along it.
    std::optional<std::pair<std::size_t, double>> nearest(double x, double y, double r) const {
        std::optional<std::pair<std::size_t, double>> best;
        double best_d2 = r * r;
        for (long cy = iy(y - r); cy <= iy(y + r); ++cy)
            for (long cx = ix(x - r); cx <= ix(x + r); ++cx)
                for (auto i : cells_[static_cast<std::size_t>(cy * nx_ + cx)]) {
                    const auto& g = segs_[i];
                    const double dx = g.x1 - g.x0, dy = g.y1 - g.y0;
                    const double l2 = dx * dx + dy * dy;
                    const double t = l2 > 0 ? std::clamp(((x - g.x0) * dx + (y - g.y0) * dy) / l2, 0.0, 1.0) : 0.0;
                    const double ex = g.x0 + t * dx - x, ey = g.y0 + t * dy - y;
                    if (ex * ex + ey * ey <= best_d2) {
                        best_d2 = ex * ex + ey * ey;
                        best = std::make_pair(static_cast<std::size_t>(i), t);
                    }
                }
        return best;
    }

    // True if the size field at (x, y), min over segments of h_s + grading * d_s, is
    // below `size`. `h_floor` is the smallest h_s.
    bool finer_than(double x, double y, double size, double grading, double h_floor) const {
        const double r = (size - h_floor) / grading;
        if (r <= 0) return false;
        for (long cy = iy(y - r); cy <= iy(y + r); ++cy)
            for (long cx = ix(x - r); cx <= ix(x + r); ++cx)
                for (auto i : cells_[static_cast<std::size_t>(cy * nx_ + cx)]) {
                    const double lim = (size - segs_[i].h) / grading;
                    if (lim > 0 && dist2(segs_[i], x, y) <= lim * lim) return true;
                }
        return false;
    }

    // Point in the region bounded by the boundary segments (even-odd), by counting the
    // crossings of a ray towards +x; each crossing is counted in the cell containing it.
    bool inside(double x, double y) const {
        const long cy = iy(y);
        if (y < y0_ || cy < 0) return false;
        // The ray ends at the first cell clear of the boundary, for points on the grid.
        const bool on_grid = x >= x0_ && x < x0_ + static_cast<double>(nx_) * cell_ && y < y0_ + static_cast<double>(ny_) * cell_;
        bool in = false;
        for (long cx = std::max(0L, ix(x)); cx < nx_; ++cx) {
            const auto c = static_cast<std::size_t>(cy * nx_ + cx);
            if (on_grid && state_[c] != kMixed) return in != (state_[c] == kInside);
            for (auto i : cells_[c]) {
                const auto& g = segs_[i];
                if (!g.boundary) continue;
                if (!crosses(g, y)) continue;
                const double xc = crossing(g, y);
                if (xc <= x) continue;
                // Each crossing is counted in one cell, found as the segments were
                // registered: bounds computed otherwise can lose a crossing on a cell edge.
                if (ix(xc) == cx) in = !in;
            }
        }
        return in;
    }

private:
    enum : std::uint8_t { kMixed, kOutside, kInside };
    static bool crosses(const Seg& g, double y) { return (g.y0 > y) != (g.y1 > y); }
    static double crossing(const Seg& g, double y) { return g.x0 + (y - g.y0) * (g.x1 - g.x0) / (g.y1 - g.y0); }
    static double dist2(const Seg& g, double x, double y) {
        const double dx = g.x1 - g.x0, dy = g.y1 - g.y0;
        const double l2 = dx * dx + dy * dy;
        double t = l2 > 0 ? ((x - g.x0) * dx + (y - g.y0) * dy) / l2 : 0.0;
        t = std::clamp(t, 0.0, 1.0);
        const double ex = g.x0 + t * dx - x, ey = g.y0 + t * dy - y;
        return ex * ex + ey * ey;
    }
    long ix(double x) const { return std::clamp(static_cast<long>(std::floor((x - x0_) / cell_)), 0L, nx_ - 1); }
    long iy(double y) const { return std::clamp(static_cast<long>(std::floor((y - y0_) / cell_)), 0L, ny_ - 1); }

    std::vector<Seg> segs_;
    double x0_, y0_, cell_;
    long nx_ = 0, ny_ = 0;
    std::vector<std::vector<std::uint32_t>> cells_;
    std::vector<std::uint8_t> state_;
};

}  // namespace kairoq::detail
