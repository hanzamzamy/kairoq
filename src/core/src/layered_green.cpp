// Electrostatic Green's function of a stratified medium (layered_green.hpp).
//
// In each slab the transform of the potential obeys phi'' = k^2 phi. Looking down from a
// height, the ratio Y = D / phi (D = eps dphi/dz) follows from the slabs below; looking up,
// Y = -D / phi from those above. A unit charge at zs gives phi(zs) = 1 / (Y_down + Y_up), and
// across a slab of thickness t the potential falls by 1 / (cosh kt + y sinh kt), with y the
// normalised ratio at the far side. The potential in space is
// G(rho) = 1 / (2 pi) int_0^inf phi(k) J0(k rho) k dk.

#include "layered_green.hpp"

#include "fem_util.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace kairoq::detail {

namespace {

constexpr double kEps0 = 8.8541878128e-12;
constexpr double kPi = 3.14159265358979323846;

// Interfaces of the slabs: z[i] is the bottom of slab i, z[n] the top of the last.
std::vector<double> interfaces(const Stratified& s) {
    std::vector<double> z(s.thickness.size() + 1, 0.0);
    for (std::size_t i = 0; i < s.thickness.size(); ++i) z[i + 1] = z[i] + s.thickness[i];
    return z;
}

// Normalised ratio across a slab segment of length t: from y at its near end to the far end.
double carry(double y, double kt) {
    const double T = std::tanh(kt);
    if (std::isinf(y)) return 1.0 / T;
    return (T + y) / (1.0 + y * T);
}

// Y_down / k (F/m) at height z: the slabs below, the ground or the half-space under them.
double y_down(const Stratified& s, const std::vector<double>& z, double k, double h) {
    double Y = s.ground ? INFINITY : kEps0 * s.epsilon_below;  // Y / k, continuous across interfaces
    for (std::size_t i = 0; i < s.thickness.size(); ++i) {
        if (h <= z[i]) break;
        const double e = kEps0 * s.epsilon_r[i];
        const double top = std::min(h, z[i + 1]);
        Y = e * carry(std::isinf(Y) ? INFINITY : Y / e, k * (top - z[i]));
    }
    if (h > z.back()) {
        const double e = kEps0 * s.epsilon_above;
        Y = e * carry(std::isinf(Y) ? INFINITY : Y / e, k * (h - z.back()));
    }
    return Y;
}

// Y_up / k (F/m) at height z: the slabs above and the half-space over them.
double y_up(const Stratified& s, const std::vector<double>& z, double k, double h) {
    double Y = kEps0 * s.epsilon_above;
    for (std::size_t i = s.thickness.size(); i-- > 0;) {
        if (h >= z[i + 1]) break;
        const double e = kEps0 * s.epsilon_r[i];
        const double bottom = std::max(h, z[i]);
        Y = e * carry(Y / e, k * (z[i + 1] - bottom));
    }
    if (h < 0 && !s.ground) {
        const double e = kEps0 * s.epsilon_below;
        Y = e * carry(Y / e, k * (0.0 - h));
    }
    return Y;
}

// Permittivity (F/m) at a height between two levels a < b (its middle decides).
double eps_at(const Stratified& s, const std::vector<double>& z, double a, double b) {
    const double m = 0.5 * (a + b);
    if (m >= z.back()) return kEps0 * s.epsilon_above;
    if (m < 0) return kEps0 * s.epsilon_below;
    for (std::size_t i = 0; i < s.thickness.size(); ++i)
        if (m < z[i + 1]) return kEps0 * s.epsilon_r[i];
    return kEps0 * s.epsilon_above;
}

// A height within rounding of an interface is taken as on it.
double snap(const std::vector<double>& z, double h) {
    for (double x : z)
        if (std::abs(h - x) <= 1e-12) return x;
    return h;
}

}  // namespace

double spectral_scaled(const Stratified& s, double k, double z, double zs) {
    const auto zi = interfaces(s);
    z = snap(zi, z);
    zs = snap(zi, zs);
    if (s.ground && (z <= 0 || zs <= 0)) return 0.0;
    const double yd = y_down(s, zi, k, zs), yu = y_up(s, zi, k, zs);
    double phi = std::isinf(yd) ? 0.0 : 1.0 / (k * (yd + yu));
    if (phi == 0.0 || z == zs) return phi;
    // Levels between zs and z: the interfaces strictly between them.
    std::vector<double> lv{zs};
    for (double x : zi)
        if ((x > std::min(z, zs)) && (x < std::max(z, zs))) lv.push_back(x);
    if (z > zs) {
        std::sort(lv.begin(), lv.end());
    } else {
        std::sort(lv.begin(), lv.end(), std::greater<double>());
    }
    lv.push_back(z);
    for (std::size_t i = 0; i + 1 < lv.size(); ++i) {
        const double a = lv[i], b = lv[i + 1], t = std::abs(b - a);
        if (t == 0) continue;
        const double e = eps_at(s, zi, std::min(a, b), std::max(a, b));
        // y at the far end, looking away from the source.
        const double Y = z > zs ? y_up(s, zi, k, b) : y_down(s, zi, k, b);
        if (std::isinf(Y)) return 0.0;
        const double y = Y / e, q = std::exp(-2.0 * k * t);
        phi *= 2.0 / ((1.0 + q) + y * (1.0 - q));  // (cosh kt + y sinh kt)^-1 times exp(kt)
    }
    return phi;
}

LayeredGreen::LayeredGreen(const Stratified& s, double z, double zs, double rho_max) {
    if (s.thickness.size() != s.epsilon_r.size()) throw std::runtime_error("layered medium: thickness and permittivity differ in length");
    const auto zi = interfaces(s);
    z = snap(zi, z);
    zs = snap(zi, zs);
    dz_ = std::abs(z - zs);
    // Smallest distance that shapes the transform: slabs, and the source and field points to
    // the interfaces.
    double t_min = INFINITY, height = zi.back();
    for (double t : s.thickness)
        if (t > 0) t_min = std::min(t_min, t);
    for (double x : zi)
        for (double p : {z, zs})
            if (std::abs(p - x) > 0) t_min = std::min(t_min, std::abs(p - x));
    if (!std::isfinite(t_min)) t_min = std::max(dz_, 1e-6);
    height = std::max({height, z, zs, t_min});
    const double k_big = 60.0 / t_min;
    c_ = k_big * spectral_scaled(s, k_big, z, zs) / (2 * kPi);
    const double C = 2 * kPi * c_;
    if (z == zs) {
        // Images of the nearest interface, from the sheet between the medium on its other
        // side (eps_o) and the slab towards it (eps_m, thickness d) with eps_f beyond:
        // k phi = C (1 + G u) / (1 - G K u), u = exp(-2 k d), G = (eps_m - eps_f) / (eps_m + eps_f)
        // (-1 for the ground), K = (eps_m - eps_o) / (eps_m + eps_o); the n-th image is
        // C (1 + K) G^n K^(n-1) at depth 2 n d. The other interfaces are left to the remainder.
        const double tiny = 1e-12;
        double up = INFINITY, down = INFINITY;
        for (double x : zi) {
            if (x > z + tiny) up = std::min(up, x - z);
            if (x < z - tiny) down = std::min(down, z - x);
        }
        if (std::isfinite(up) || std::isfinite(down)) {
            const bool above = up <= down;
            const double d = above ? up : down;
            const double em = above ? eps_at(s, zi, z, z + d) : eps_at(s, zi, z - d, z);
            const double eo = above ? eps_at(s, zi, z - tiny, z) : eps_at(s, zi, z, z + tiny);
            double G;
            if (!above && s.ground && std::abs(z - d) < tiny) {
                G = -1.0;
            } else {
                const double ef = above ? eps_at(s, zi, z + d, z + d + tiny) : eps_at(s, zi, z - d - tiny, z - d);
                G = (em - ef) / (em + ef);
            }
            const double K = (em - eo) / (em + eo);
            double term = c_ * (1 + K) * G;  // n = 1
            for (int n = 1; n <= 60 && std::abs(term) > 1e-7 * std::abs(c_); ++n) {
                images_.push_back({term, 2 * n * d});
                term *= G * K;
            }
        }
    }
    auto f = [&](double k) {
        double v = k * spectral_scaled(s, k, z, zs) - C;
        for (const auto& im : images_) v -= 2 * kPi * im.coefficient * std::exp(-k * im.depth);
        return v * std::exp(-k * dz_);
    };

    // 8-point Gauss-Legendre on [0, 1].
    static const double gx[8] = {0.0198550717512319, 0.1016667612931866, 0.2372337950418355, 0.4082826787521751,
                                 0.5917173212478249, 0.7627662049581645, 0.8983332387068134, 0.9801449282487681};
    static const double gw[8] = {0.0506142681451881, 0.1111905172266872, 0.1568533229389436, 0.1813418916891810,
                                 0.1813418916891810, 0.1568533229389436, 0.1111905172266872, 0.0506142681451881};
    auto integrate = [&](double rho) {
        const double w_max = std::min(0.5 / t_min, 0.25 / height);
        double sum = 0, k = 0, w = rho > 0 ? std::min(kPi / rho, w_max) : w_max;
        int quiet = 0;
        for (int p = 0; p < 400000 && quiet < 4; ++p) {
            // Panels grow where the transform no longer varies on their scale.
            const double wp = rho > 0 ? std::min(kPi / rho, std::max(w_max, 0.5 * k)) : std::max(w_max, 0.25 * k);
            w = wp;
            double part = 0, peak = 0;
            for (int g = 0; g < 8; ++g) {
                const double kk = k + w * gx[g];
                const double v = f(kk);
                peak = std::max(peak, std::abs(v));
                part += gw[g] * w * v * (rho > 0 ? std::cyl_bessel_j(0.0, kk * rho) : 1.0);
            }
            sum += part;
            k += w;
            quiet = peak < 1e-12 * std::abs(C) ? quiet + 1 : 0;
        }
        return sum / (2 * kPi);
    };
    rem0_ = integrate(0.0);
    rho_min_ = 0.02 * t_min;
    const double hi = std::max(rho_max, 10 * rho_min_);
    const int n = static_cast<int>(std::ceil(32 * std::log10(hi / rho_min_))) + 1;
    log_rho_.resize(n);
    rem_.resize(n);
    for (int i = 0; i < n; ++i) log_rho_[i] = std::log(rho_min_) + (std::log(hi) - std::log(rho_min_)) * i / (n - 1);
    a_ = height;
    parallel_for(static_cast<std::size_t>(n), [&](std::size_t i) {
        const double rho = std::exp(log_rho_[i]);
        rem_[i] = integrate(rho) * std::sqrt(rho * rho + a_ * a_);
    });
}

double LayeredGreen::remainder(double rho) const {
    const double scale = 1.0 / std::sqrt(rho * rho + a_ * a_);
    if (rho <= rho_min_) {
        const double u = rho / rho_min_;
        return rem0_ + u * (rem_.front() * scale - rem0_);
    }
    const double x = std::log(rho);
    const std::size_t n = log_rho_.size();
    if (x >= log_rho_.back()) return rem_.back() * scale;
    const std::size_t i = std::upper_bound(log_rho_.begin(), log_rho_.end(), x) - log_rho_.begin();
    const double u = (x - log_rho_[i - 1]) / (log_rho_[i] - log_rho_[i - 1]);
    // Catmull-Rom between i - 1 and i, ends extended linearly.
    const double p1 = rem_[i - 1], p2 = rem_[i];
    const double p0 = i >= 2 ? rem_[i - 2] : 2 * p1 - p2, p3 = i + 1 < n ? rem_[i + 1] : 2 * p2 - p1;
    const double v = p1 + 0.5 * u * (p2 - p0 + u * (2 * p0 - 5 * p1 + 4 * p2 - p3 + u * (3 * (p1 - p2) + p3 - p0)));
    return v * scale;
}

double LayeredGreen::operator()(double rho) const {
    const double R = std::sqrt(rho * rho + dz_ * dz_);
    double v = c_ / R + remainder(rho);
    for (const auto& im : images_) v += im.coefficient / std::sqrt(rho * rho + (dz_ + im.depth) * (dz_ + im.depth));
    return v;
}

}  // namespace kairoq::detail
