// De-embedding of S-parameters (docs/solvers/rf.md).
#include "kairoq/deembed.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace fs = std::filesystem;

namespace {

constexpr double kLight = 299792458.0;

// Inverse of a small complex matrix (Gauss-Jordan, partial pivoting).
std::vector<cplx> inverse(std::vector<cplx> a, std::size_t n) {
    std::vector<cplx> inv(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) inv[i * n + i] = 1.0;
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::abs(a[r * n + c]) > std::abs(a[p * n + c])) p = r;
        if (std::abs(a[p * n + c]) == 0.0) throw std::runtime_error("de-embedding: singular matrix");
        for (std::size_t k = 0; k < n; ++k) {
            std::swap(a[c * n + k], a[p * n + k]);
            std::swap(inv[c * n + k], inv[p * n + k]);
        }
        const cplx d = a[c * n + c];
        for (std::size_t k = 0; k < n; ++k) {
            a[c * n + k] /= d;
            inv[c * n + k] /= d;
        }
        for (std::size_t r = 0; r < n; ++r) {
            if (r == c) continue;
            const cplx f = a[r * n + c];
            if (f == 0.0) continue;
            for (std::size_t k = 0; k < n; ++k) {
                a[r * n + k] -= f * a[c * n + k];
                inv[r * n + k] -= f * inv[c * n + k];
            }
        }
    }
    return inv;
}

std::array<cplx, 4> abcd_from_s(const std::array<cplx, 4>& s, double z0) {
    const cplx s11 = s[0], s12 = s[1], s21 = s[2], s22 = s[3];
    const cplx d = 2.0 * s21;
    return {((1.0 + s11) * (1.0 - s22) + s12 * s21) / d, z0 * ((1.0 + s11) * (1.0 + s22) - s12 * s21) / d,
            ((1.0 - s11) * (1.0 - s22) - s12 * s21) / (d * z0), ((1.0 - s11) * (1.0 + s22) + s12 * s21) / d};
}

std::array<cplx, 4> inverse2(const std::array<cplx, 4>& m) {
    const cplx det = m[0] * m[3] - m[1] * m[2];
    return {m[3] / det, -m[1] / det, -m[2] / det, m[0] / det};
}

void check_port(const Touchstone& t, std::size_t port) {
    if (port >= t.ports) throw std::runtime_error("de-embedding: port " + std::to_string(port + 1) + " of a " + std::to_string(t.ports) + "-port");
}

}  // namespace

Touchstone renormalise(const Touchstone& t, double z) {
    if (!(z > 0)) throw std::runtime_error("de-embedding: the reference must be positive");
    const std::size_t n = t.ports;
    const double rho = (z - t.reference_ohm) / (z + t.reference_ohm);
    Touchstone out = t;
    out.reference_ohm = z;
    for (std::size_t k = 0; k < t.s.size(); ++k) {
        // S' = (I - rho S)^-1 (S - rho I).
        std::vector<cplx> a(n * n), b(n * n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                a[i * n + j] = (i == j ? 1.0 : 0.0) - rho * t.s[k][i * n + j];
                b[i * n + j] = t.s[k][i * n + j] - (i == j ? rho : 0.0);
            }
        const auto ai = inverse(a, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                cplx v = 0;
                for (std::size_t m = 0; m < n; ++m) v += ai[i * n + m] * b[m * n + j];
                out.s[k][i * n + j] = v;
            }
    }
    return out;
}

Touchstone connect_fixture(const Touchstone& t, std::size_t p, const std::vector<std::array<cplx, 4>>& x) {
    check_port(t, p);
    if (x.size() != t.s.size()) throw std::runtime_error("de-embedding: the fixture has another number of frequencies");
    const std::size_t n = t.ports;
    Touchstone out = t;
    for (std::size_t k = 0; k < t.s.size(); ++k) {
        const auto& s = t.s[k];
        auto& o = out.s[k];
        const cplx x11 = x[k][0], x12 = x[k][1], x21 = x[k][2], x22 = x[k][3];
        const cplx den = 1.0 - x22 * s[p * n + p];
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                if (i == p && j == p) o[i * n + j] = x11 + x12 * x21 * s[p * n + p] / den;
                else if (i == p) o[i * n + j] = x12 * s[p * n + j] / den;
                else if (j == p) o[i * n + j] = s[i * n + p] * x21 / den;
                else o[i * n + j] = s[i * n + j] + s[i * n + p] * x22 * s[p * n + j] / den;
            }
    }
    return out;
}

Touchstone add_line(const Touchstone& t, std::size_t port, const LineFixture& l) {
    if (!(l.z_ohm > 0) || !(l.epsilon_eff >= 1.0)) throw std::runtime_error("de-embedding: a line needs z_ohm > 0 and epsilon_eff >= 1");
    std::vector<std::array<cplx, 4>> x;
    for (double f : t.frequency_hz) {
        const double alpha = l.loss_db_per_m_1ghz / (20 / std::log(10.0)) * std::sqrt(f / 1e9);  // Np/m
        const cplx gamma(alpha, 2 * kPi * f * std::sqrt(l.epsilon_eff) / kLight);
        const auto s = s_from_abcd(line_abcd(l.z_ohm, gamma, l.length_m), t.reference_ohm);
        x.push_back({s[0], s[1], s[2], s[3]});
    }
    return connect_fixture(t, port, x);
}

std::vector<std::array<cplx, 4>> half_thru_abcd(const Touchstone& thru, const std::vector<double>& freqs) {
    if (thru.ports != 2) throw std::runtime_error("de-embedding: a 2x-thru is a two-port");
    const double z0 = thru.reference_ohm;
    std::vector<std::array<cplx, 4>> out;
    std::array<cplx, 4> prev = {1.0, 0.0, 0.0, 1.0};  // the half at zero frequency
    for (double f : freqs) {
        const auto v = thru.at(f);
        const auto m = abcd_from_s({v[0], v[1], v[2], v[3]}, z0);
        // Square roots of a 2 x 2 matrix: (M + s I) / t with s^2 = det M, t^2 = tr M + 2 s.
        const cplx det = m[0] * m[3] - m[1] * m[2];
        std::array<cplx, 4> best{};
        double dist = INFINITY;
        for (const double ss : {1.0, -1.0}) {
            const cplx s = ss * std::sqrt(det);
            const cplx t2 = m[0] + m[3] + 2.0 * s;
            if (std::abs(t2) < 1e-12) continue;
            for (const double ts : {1.0, -1.0}) {
                const cplx t = ts * std::sqrt(t2);
                const std::array<cplx, 4> r = {(m[0] + s) / t, m[1] / t, m[2] / t, (m[3] + s) / t};
                // Distance to the previous root, B and C made dimensionless.
                const double d = std::abs(r[0] - prev[0]) + std::abs((r[1] - prev[1]) / z0) + std::abs((r[2] - prev[2]) * z0) + std::abs(r[3] - prev[3]);
                if (d < dist) {
                    dist = d;
                    best = r;
                }
            }
        }
        if (!(dist < INFINITY)) throw std::runtime_error("de-embedding: the 2x-thru has no square root at " + std::to_string(f) + " Hz");
        out.push_back(best);
        prev = best;
    }
    return out;
}

Touchstone deembed_thru(const Touchstone& t, std::size_t port, const Touchstone& thru) {
    check_port(t, port);
    if (t.frequency_hz.back() > thru.frequency_hz.back() * (1 + 1e-9) || t.frequency_hz.front() < thru.frequency_hz.front() * (1 - 1e-9))
        throw std::runtime_error("de-embedding: the 2x-thru does not cover the frequencies");
    const Touchstone th = std::abs(thru.reference_ohm - t.reference_ohm) > 1e-9 ? renormalise(thru, t.reference_ohm) : thru;
    const auto half = half_thru_abcd(th, t.frequency_hz);
    std::vector<std::array<cplx, 4>> x;
    for (const auto& h : half) {
        const auto s = s_from_abcd(inverse2(h), t.reference_ohm);
        x.push_back({s[0], s[1], s[2], s[3]});
    }
    return connect_fixture(t, port, x);
}

std::string touchstone_text(const Touchstone& t, const std::string& comment) {
    std::ostringstream o;
    if (!comment.empty()) o << "! " << comment << "\n";
    char buf[64];
    std::snprintf(buf, sizeof buf, "# Hz S RI R %.9g\n", t.reference_ohm);
    o << buf;
    const std::size_t n = t.ports;
    for (std::size_t k = 0; k < t.frequency_hz.size(); ++k) {
        std::snprintf(buf, sizeof buf, "%.9g", t.frequency_hz[k]);
        o << buf;
        std::size_t count = 0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                const cplx v = n == 2 ? t.s[k][j * n + i] : t.s[k][i * n + j];  // two-ports: S11 S21 S12 S22
                std::snprintf(buf, sizeof buf, " %.12g %.12g", v.real(), v.imag());
                o << buf;
                if (n > 2 && ++count % 4 == 0 && !(i == n - 1 && j == n - 1)) o << "\n";
            }
        o << "\n";
    }
    return o.str();
}

DeembedResult deembed(std::string_view yaml, const std::string& base_dir) {
    const YAML::Node root = YAML::Load(std::string(yaml));
    auto path = [&](const std::string& p) { return fs::path(p).is_absolute() ? p : (fs::path(base_dir) / p).string(); };
    for (const auto& kv : root) {
        const std::string k = kv.first.as<std::string>();
        if (k != "input" && k != "operations" && k != "output") throw std::runtime_error("deembed: unknown key '" + k + "'");
    }
    if (!root["input"]) throw std::runtime_error("deembed: needs an input (a Touchstone file)");
    DeembedResult r;
    r.data = read_touchstone(path(root["input"].as<std::string>()));
    char line[200];
    std::size_t i = 0;
    for (const auto& op : root["operations"]) {
        const std::string where = "deembed operation " + std::to_string(++i);
        if (!op.IsMap() || op.size() != 1) throw std::runtime_error(where + ": one of renormalise_ohm, line, thru");
        const std::string kind = op.begin()->first.as<std::string>();
        const YAML::Node v = op.begin()->second;
        if (kind == "renormalise_ohm") {
            r.data = renormalise(r.data, v.as<double>());
            std::snprintf(line, sizeof line, "reference changed to %.6g ohm", v.as<double>());
        } else if (kind == "line") {
            for (const auto& kv : v) {
                const std::string k = kv.first.as<std::string>();
                if (k != "port" && k != "length_mm" && k != "z_ohm" && k != "epsilon_eff" && k != "loss_db_per_m_1ghz")
                    throw std::runtime_error(where + ": unknown key '" + k + "'");
            }
            if (!v["port"] || !v["length_mm"] || !v["epsilon_eff"]) throw std::runtime_error(where + ": a line needs port, length_mm and epsilon_eff");
            LineFixture l;
            l.length_m = v["length_mm"].as<double>() * 1e-3;
            l.z_ohm = v["z_ohm"] ? v["z_ohm"].as<double>() : r.data.reference_ohm;
            l.epsilon_eff = v["epsilon_eff"].as<double>();
            l.loss_db_per_m_1ghz = v["loss_db_per_m_1ghz"] ? v["loss_db_per_m_1ghz"].as<double>() : 0.0;
            const std::size_t port = v["port"].as<std::size_t>();
            if (port == 0) throw std::runtime_error(where + ": ports count from 1");
            r.data = add_line(r.data, port - 1, l);
            std::snprintf(line, sizeof line, "port %zu: line of %.6g mm %s (%.6g ohm, effective permittivity %.6g)", port,
                          std::abs(l.length_m) * 1e3, l.length_m < 0 ? "removed" : "added", l.z_ohm, l.epsilon_eff);
        } else if (kind == "thru") {
            for (const auto& kv : v) {
                const std::string k = kv.first.as<std::string>();
                if (k != "port" && k != "file") throw std::runtime_error(where + ": unknown key '" + k + "'");
            }
            if (!v["port"] || !v["file"]) throw std::runtime_error(where + ": a thru needs port and file");
            const std::size_t port = v["port"].as<std::size_t>();
            if (port == 0) throw std::runtime_error(where + ": ports count from 1");
            r.data = deembed_thru(r.data, port - 1, read_touchstone(path(v["file"].as<std::string>())));
            std::snprintf(line, sizeof line, "port %zu: half of the 2x-thru %s removed", port, v["file"].as<std::string>().c_str());
        } else {
            throw std::runtime_error(where + ": unknown operation '" + kind + "'");
        }
        r.steps.push_back(line);
    }
    if (root["output"]) {
        r.output = path(root["output"].as<std::string>());
        std::ofstream f(r.output);
        f << touchstone_text(r.data, "kairoq deembed of " + root["input"].as<std::string>());
        if (!f) throw std::runtime_error("deembed: cannot write " + r.output);
    }
    return r;
}

}  // namespace kairoq
