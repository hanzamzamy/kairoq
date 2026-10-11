// Time-domain response of a channel with drivers and receivers (docs/solvers/timedomain.md).
//
// The channel is its scattering matrix on the linear sweep, turned into band-limited
// impulse responses: the waves leaving the ports are the convolution of those entering,
// b = s * a. At each time step the waves entering are solved, by Newton's method, together
// with the terminations at the ports: sources, loads and IBIS buffers with their package
// and die capacitance (backward Euler).
#include "kairoq/timedomain.hpp"
#include "kairoq/ibis.hpp"

#include "fem_util.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

using namespace detail;
namespace fs = std::filesystem;

// Linear feedback shift register of x^n + x^m + 1, from all ones.
std::vector<int> prbs(int n, int m, std::size_t count) {
    std::uint64_t reg = (std::uint64_t(1) << n) - 1;
    std::vector<int> out;
    for (std::size_t i = 0; i < count; ++i) {
        const int bit = static_cast<int>(((reg >> (n - 1)) ^ (reg >> (m - 1))) & 1);
        reg = ((reg << 1) | std::uint64_t(bit)) & ((std::uint64_t(1) << n) - 1);
        out.push_back(bit);
    }
    return out;
}

// A port's termination: its current into the network at a port voltage, given its own
// unknowns, and their equations.
struct Termination {
    virtual ~Termination() = default;
    virtual std::size_t unknowns() const { return 0; }
    // Current into the network at voltage v; internal residuals into r.
    virtual double current(double v, const double* x, double* r, double t) const = 0;
    virtual void accept(double v, const double* x) { (void)v, (void)x; }
    virtual void initial(double* x) const { (void)x; }
};

struct Open : Termination {
    double current(double, const double*, double*, double) const override { return 0.0; }
};

// A source behind a resistance, or a load (resistor to a voltage, capacitor), backward Euler.
struct Linear : Termination {
    double r = 0, c = 0, dt = 1;
    std::function<double(double)> source;  // V at time t; empty: a load
    double v_term = 0, v_prev = 0;
    double current(double v, const double*, double*, double t) const override {
        double i = 0;
        if (source) i += (source(t) - v) / std::max(r, 1e-6);
        else if (r > 0) i -= (v - v_term) / r;
        i -= c * (v - v_prev) / dt;
        return i;
    }
    void accept(double v, const double*) override { v_prev = v; }
};

// An IBIS buffer behind its package: pin (the port) - R, L - die (C_comp, the buffer's
// current); C_pkg at the pin. Unknowns: the die voltage and the inductor's current.
struct Ibis : Termination {
    const IbisModel* model = nullptr;
    int corner = 0;
    double rp = 0, lp = 0, cp = 0, dt = 1;
    std::vector<std::pair<double, int>> edges;  // time, rising (1) or falling (0)
    IbisSwitching up, down;
    bool driver = false;
    double v_pin_prev = 0, v_die_prev = 0, il_prev = 0;
    std::size_t unknowns() const override { return 2; }
    void weights(double t, double& ku, double& kd) const {
        ku = 0;
        kd = driver ? 1 : 0;
        if (!driver) return;
        for (const auto& [te, rising] : edges) {
            if (t < te) break;
            const auto& s = rising ? up : down;
            const double tl = t - te;
            if (s.time_s.empty()) continue;
            if (tl >= s.time_s.back()) {
                ku = s.ku.back();
                kd = s.kd.back();
                continue;
            }
            const auto it = std::upper_bound(s.time_s.begin(), s.time_s.end(), tl);
            const std::size_t k = std::max<std::size_t>(1, static_cast<std::size_t>(it - s.time_s.begin()));
            const double w = (tl - s.time_s[k - 1]) / (s.time_s[k] - s.time_s[k - 1]);
            ku = s.ku[k - 1] + w * (s.ku[k] - s.ku[k - 1]);
            kd = s.kd[k - 1] + w * (s.kd[k] - s.kd[k - 1]);
        }
    }
    double current(double v, const double* x, double* r, double t) const override {
        const double vd = x[0], il = x[1];
        double ku, kd;
        weights(t, ku, kd);
        const double cc = model->c_comp.at(corner);
        // At the die: the inductor's current feeds the buffer and C_comp.
        r[0] = il - model->current(vd, ku, kd, corner) - cc * (vd - v_die_prev) / dt;
        // Along the package.
        r[1] = v - vd - rp * il - lp * (il - il_prev) / dt;
        // Into the network: minus what the pin's capacitance and the package take.
        return -il - cp * (v - v_pin_prev) / dt;
    }
    void accept(double v, const double* x) override {
        v_pin_prev = v;
        v_die_prev = x[0];
        il_prev = x[1];
    }
};

IbisSwitching from_corner(const IbisModel& m, bool rising, int corner, std::vector<std::string>& w) {
    return switching_coefficients(m, rising, corner, &w);
}


}  // namespace

std::vector<int> Stimulus::sequence() const {
    if (pattern == "edge") return {1};
    if (pattern == "pulse") return {1, 0};
    static const std::map<std::string, std::pair<int, int>> taps = {
        {"prbs7", {7, 6}}, {"prbs9", {9, 5}}, {"prbs15", {15, 14}}, {"prbs23", {23, 18}}, {"prbs31", {31, 28}}};
    if (auto it = taps.find(pattern); it != taps.end()) return prbs(it->second.first, it->second.second, bits);
    std::vector<int> out;
    for (char c : pattern) {
        if (c == '0' || c == '1') out.push_back(c - '0');
        else throw std::runtime_error("stimulus: pattern is edge, pulse, prbsN or a string of 0 and 1");
    }
    if (out.empty()) throw std::runtime_error("stimulus: empty pattern");
    return out;
}

SimulationResult simulate(const ChannelResult& ch, const SimulationSpec& spec, const std::string& base_dir) {
    const auto t0 = Clock::now();
    SimulationResult res;
    res.name = spec.name;
    res.channel = ch.name;
    res.ports = ch.ports;
    const std::size_t np = ch.ports.size(), nf = ch.frequency_hz.size();
    if (nf < 2) throw std::runtime_error("simulation " + spec.name + ": the channel has no sweep");
    const double fmax = ch.frequency_hz.back(), df = ch.frequency_hz[1] - ch.frequency_hz[0];
    const double sz = std::sqrt(ch.impedance_ohm);  // the channel's reference
    const double dt = 1.0 / (2.0 * double(spec.oversample) * fmax);
    res.dt = dt;
    // Duration: the stimulus and 2 ns after it, at most half the sweep's period.
    double stim = 0;
    for (const auto& d : spec.drivers) stim = std::max(stim, d.stimulus.delay_s + double(d.stimulus.sequence().size()) / d.stimulus.bit_rate);
    double duration = spec.duration_s > 0 ? spec.duration_s : stim + 2e-9;
    const double period = 1.0 / df;
    if (duration > 0.5 * period) {
        if (stim > 0.5 * period || spec.duration_s > 0)
            res.warnings.push_back("duration limited to half the sweep's period (" + std::to_string(0.5 * period * 1e9) +
                                   " ns); more sweep points allow longer");
        duration = 0.5 * period;
    }
    const std::size_t nt = static_cast<std::size_t>(std::ceil(duration / dt)) + 1;

    // Impulse responses: the sweep tapered to zero over its last quarter, the lowest
    // frequency's value extended to DC. Band-limited, an instantaneous part (the
    // reflection at a port's own junction) spreads to either side of t = 0; the part
    // before it (within 3 / f_max, wrapped to the end of the period) is added at t = 0, so
    // that the response stays causal and keeps its area.
    const std::size_t nperiod = static_cast<std::size_t>(std::llround(period / dt));
    const std::size_t nw = std::min<std::size_t>(nperiod / 4, static_cast<std::size_t>(std::ceil(3.0 / (fmax * dt))));
    std::vector<std::vector<double>> s(np * np, std::vector<double>(nt, 0.0));
    parallel_for(np * np, [&](std::size_t ij) {
        auto& out = s[ij];
        const double s0 = ch.s[0][ij].real();
        std::vector<double> before(nw, s0);  // times -nw dt .. -dt
        for (std::size_t n = 0; n < nt; ++n) out[n] = s0;
        for (std::size_t k = 0; k < nf; ++k) {
            const double f = ch.frequency_hz[k];
            const double w = f <= 0.75 * fmax ? 1.0 : 0.5 * (1 + std::cos(kPi * (f - 0.75 * fmax) / (0.25 * fmax)));
            const cplx v = 2.0 * w * ch.s[k][ij];
            const cplx rot = std::exp(cplx(0, 2 * kPi * f * dt));
            cplx ph = 1.0;
            for (std::size_t n = 0; n < nt; ++n) {
                out[n] += (v * ph).real();
                ph *= rot;
            }
            ph = std::exp(cplx(0, -2 * kPi * f * dt * double(nw)));
            for (std::size_t n = 0; n < nw; ++n) {
                before[n] += (v * ph).real();
                ph *= rot;
            }
        }
        for (double x : before) out[0] += x;
        for (auto& x : out) x *= df * dt;
    });

    // Terminations.
    std::vector<std::unique_ptr<Termination>> term(np);
    for (std::size_t p = 0; p < np; ++p) term[p] = std::make_unique<Open>();
    std::map<std::string, IbisFile> files;
    auto ibis = [&](const std::string& file, const std::string& model) -> const IbisModel* {
        fs::path path(file);
        if (path.is_relative()) path = fs::path(base_dir) / path;
        auto it = files.find(path.string());
        if (it == files.end()) {
            it = files.emplace(path.string(), read_ibis(path.string())).first;
            for (const auto& w : it->second.warnings) res.warnings.push_back(file + ": " + w);
        }
        const IbisModel* m = it->second.model(model);
        if (!m) throw std::runtime_error(file + " has no model " + model);
        return m;
    };
    auto package = [&](Ibis& t, const std::string& file, bool use, int corner) {
        if (!use) return;
        const auto& pk = files.at((fs::path(file).is_relative() ? fs::path(base_dir) / file : fs::path(file)).string()).package;
        t.rp = pk.r.at(corner);
        t.lp = pk.l.at(corner);
        t.cp = pk.c.at(corner);
    };
    for (const auto& d : spec.drivers) {
        if (d.port >= np) throw std::runtime_error("simulation " + spec.name + ": driver port out of range");
        const auto bits = d.stimulus.sequence();
        const double ui = 1.0 / d.stimulus.bit_rate;
        // Transitions: (time, rising).
        std::vector<std::pair<double, int>> edges;
        int level = 0;
        for (std::size_t k = 0; k < bits.size(); ++k)
            if (bits[k] != level) {
                edges.push_back({d.stimulus.delay_s + double(k) * ui, bits[k]});
                level = bits[k];
            }
        if (!d.ibis_model.empty()) {
            auto t = std::make_unique<Ibis>();
            t->model = ibis(d.ibis_file, d.ibis_model);
            t->corner = d.corner;
            t->dt = dt;
            t->driver = true;
            t->edges = edges;
            t->up = from_corner(*t->model, true, d.corner, res.warnings);
            t->down = from_corner(*t->model, false, d.corner, res.warnings);
            package(*t, d.ibis_file, d.package, d.corner);
            term[d.port] = std::move(t);
        } else {
            auto t = std::make_unique<Linear>();
            t->r = d.r_ohm;
            t->dt = dt;
            const double rise = 1.25 * d.rise_s, fall = 1.25 * (d.fall_s >= 0 ? d.fall_s : d.rise_s);
            const double lo = d.v_low, hi = d.v_high;
            t->source = [edges, rise, fall, lo, hi](double tt) {
                double v = lo;
                for (const auto& [te, up] : edges) {
                    if (tt < te) break;
                    const double span = up ? rise : fall;
                    const double frac = span > 0 ? std::min(1.0, (tt - te) / span) : 1.0;
                    const double from = v, to = up ? hi : lo;
                    v = from + frac * (to - from);
                }
                return v;
            };
            if (d.rise_s < 0.9 / fmax)
                res.warnings.push_back("driver at " + ch.ports[d.port] + ": its edges are faster than the sweep to " +
                                       std::to_string(fmax * 1e-9) + " GHz resolves");
            term[d.port] = std::move(t);
        }
    }
    for (const auto& r : spec.receivers) {
        if (r.port >= np) throw std::runtime_error("simulation " + spec.name + ": receiver port out of range");
        if (!r.ibis_model.empty()) {
            auto t = std::make_unique<Ibis>();
            t->model = ibis(r.ibis_file, r.ibis_model);
            t->corner = r.corner;
            t->dt = dt;
            package(*t, r.ibis_file, r.package, r.corner);
            term[r.port] = std::move(t);
        } else {
            auto t = std::make_unique<Linear>();
            t->r = r.r_ohm;
            t->c = r.c_f;
            t->v_term = r.v_term;
            t->dt = dt;
            term[r.port] = std::move(t);
        }
    }
    // Unknowns: the incident waves, then each termination's own.
    std::vector<std::size_t> offset(np);
    std::size_t nx = np;
    for (std::size_t p = 0; p < np; ++p) {
        offset[p] = nx;
        nx += term[p]->unknowns();
    }
    std::vector<double> x(nx, 0.0);
    std::vector<std::vector<double>> a(np, std::vector<double>(nt, 0.0));
    res.voltage.assign(np, std::vector<double>(nt, 0.0));
    std::vector<double> h(np), r(nx), r2(nx), J(nx * nx);
    std::size_t failed = 0;
    for (std::size_t n = 0; n < nt; ++n) {
        const double t = double(n) * dt;
        res.time_s.push_back(t);
        // History of the convolution.
        for (std::size_t i = 0; i < np; ++i) {
            double acc = 0;
            for (std::size_t j = 0; j < np; ++j) {
                const auto& sij = s[i * np + j];
                const auto& aj = a[j];
                for (std::size_t m = 1; m <= n; ++m) acc += sij[m] * aj[n - m];
            }
            h[i] = acc;
        }
        auto residual = [&](const std::vector<double>& xv, std::vector<double>& out) {
            for (std::size_t i = 0; i < np; ++i) {
                double b = h[i];
                for (std::size_t j = 0; j < np; ++j) b += s[i * np + j][0] * xv[j];
                const double v = sz * (xv[i] + b), cur = (xv[i] - b) / sz;
                const double* own = term[i]->unknowns() ? &xv[offset[i]] : nullptr;
                double* rown = term[i]->unknowns() ? &out[offset[i]] : nullptr;
                out[i] = cur - term[i]->current(v, own, rown, t);
            }
        };
        bool ok = false;
        for (int it = 0; it < 30; ++it) {
            residual(x, r);
            double norm = 0;
            for (double v : r) norm = std::max(norm, std::abs(v));
            if (norm < 1e-12) {
                ok = true;
                break;
            }
            for (std::size_t c = 0; c < nx; ++c) {
                auto xp = x;
                const double e = 1e-7 * std::max(1.0, std::abs(x[c]));
                xp[c] += e;
                residual(xp, r2);
                for (std::size_t q = 0; q < nx; ++q) J[q * nx + c] = (r2[q] - r[q]) / e;
            }
            // Solve J dx = -r.
            std::vector<double> A = J, bvec(nx);
            for (std::size_t q = 0; q < nx; ++q) bvec[q] = -r[q];
            for (std::size_t c = 0; c < nx; ++c) {
                std::size_t piv = c;
                for (std::size_t q = c + 1; q < nx; ++q)
                    if (std::abs(A[q * nx + c]) > std::abs(A[piv * nx + c])) piv = q;
                if (A[piv * nx + c] == 0) break;
                if (piv != c) {
                    for (std::size_t k = 0; k < nx; ++k) std::swap(A[c * nx + k], A[piv * nx + k]);
                    std::swap(bvec[c], bvec[piv]);
                }
                for (std::size_t q = c + 1; q < nx; ++q) {
                    const double f = A[q * nx + c] / A[c * nx + c];
                    for (std::size_t k = c; k < nx; ++k) A[q * nx + k] -= f * A[c * nx + k];
                    bvec[q] -= f * bvec[c];
                }
            }
            double step = 0;
            for (std::size_t c = nx; c-- > 0;) {
                double v = bvec[c];
                for (std::size_t k = c + 1; k < nx; ++k) v -= A[c * nx + k] * bvec[k];
                bvec[c] = v / A[c * nx + c];
                step = std::max(step, std::abs(bvec[c]));
            }
            for (std::size_t c = 0; c < nx; ++c) x[c] += bvec[c];
            if (step < 1e-12) {
                ok = true;
                break;
            }
        }
        if (!ok) ++failed;
        for (std::size_t i = 0; i < np; ++i) {
            a[i][n] = x[i];
            double b = h[i];
            for (std::size_t j = 0; j < np; ++j) b += s[i * np + j][0] * x[j];
            const double v = sz * (x[i] + b);
            res.voltage[i][n] = v;
            term[i]->accept(v, term[i]->unknowns() ? &x[offset[i]] : nullptr);
        }
    }
    if (failed) res.warnings.push_back(std::to_string(failed) + " time steps did not converge");
    // Metrics: extremes; for single edges, the overshoot beyond the final value and the
    // undershoot below the initial (or final) one.
    bool edges_only = true;
    for (const auto& d : spec.drivers) edges_only &= d.stimulus.pattern == "edge";
    for (std::size_t p = 0; p < np; ++p) {
        PortMetrics m;
        m.port = ch.ports[p];
        const auto& v = res.voltage[p];
        m.v_min = *std::min_element(v.begin(), v.end());
        m.v_max = *std::max_element(v.begin(), v.end());
        m.v_final = v.back();
        if (edges_only) {
            m.overshoot = std::max(0.0, m.v_max - std::max(m.v_final, v.front()));
            m.undershoot = std::max(0.0, std::min(m.v_final, v.front()) - m.v_min);
        }
        res.metrics.push_back(m);
    }
    res.seconds = since(t0);
    return res;
}

std::string simulations_summary(const std::vector<SimulationResult>& rs) {
    std::ostringstream o;
    char line[256];
    for (const auto& r : rs) {
        std::snprintf(line, sizeof line, "simulation %s (channel %s): %zu steps of %.2f ps, %.1f s\n", r.name.c_str(), r.channel.c_str(),
                      r.time_s.size(), r.dt * 1e12, r.seconds);
        o << line;
        for (const auto& m : r.metrics) {
            std::snprintf(line, sizeof line, "  %-16s min %8.4f V  max %8.4f V  final %8.4f V  overshoot %7.4f V  undershoot %7.4f V\n",
                          m.port.c_str(), m.v_min, m.v_max, m.v_final, m.overshoot, m.undershoot);
            o << line;
        }
        for (const auto& w : r.warnings) o << "  warning: " << w << "\n";
    }
    return o.str();
}

std::string simulations_json(const std::vector<SimulationResult>& rs) {
    using nlohmann::json;
    json out = json::array();
    for (const auto& r : rs) {
        json m = json::array();
        for (const auto& x : r.metrics)
            m.push_back({{"port", x.port}, {"v_min", x.v_min}, {"v_max", x.v_max}, {"v_final", x.v_final}, {"overshoot", x.overshoot},
                         {"undershoot", x.undershoot}});
        out.push_back({{"name", r.name}, {"channel", r.channel}, {"ports", r.ports}, {"time_s", r.time_s}, {"voltage", r.voltage},
                       {"metrics", m}, {"dt", r.dt}, {"seconds", r.seconds}, {"warnings", r.warnings}});
    }
    return out.dump();
}

}  // namespace kairoq
