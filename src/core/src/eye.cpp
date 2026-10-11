// Eye diagrams by superposition of the single-bit response (docs/solvers/eye.md).
//
// The transfer from the driver's source to the receiver follows from the channel's
// scattering matrix with the terminations at each frequency; with the CTLE and the
// spectrum of one bit it gives the single-bit response. The eye of a bit pattern is its
// superposition; the statistical eye takes, at each sampling phase, the distribution of the
// intersymbol interference from the cursors, Gaussian noise, and jitter as a mixture over
// sampling offsets.
#include "kairoq/eye.hpp"
#include "kairoq/timedomain.hpp"

#include "fem_util.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

using namespace detail;

// Solves A x = b for small complex systems (partial pivoting).
std::vector<cplx> csolve(std::vector<cplx> A, std::vector<cplx> b, std::size_t n) {
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::abs(A[r * n + c]) > std::abs(A[p * n + c])) p = r;
        if (std::abs(A[p * n + c]) == 0.0) throw std::runtime_error("eye: singular network");
        for (std::size_t k = 0; k < n; ++k) std::swap(A[c * n + k], A[p * n + k]);
        std::swap(b[c], b[p]);
        for (std::size_t r = c + 1; r < n; ++r) {
            const cplx f = A[r * n + c] / A[c * n + c];
            for (std::size_t k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
            b[r] -= f * b[c];
        }
    }
    for (std::size_t c = n; c-- > 0;) {
        for (std::size_t k = c + 1; k < n; ++k) b[c] -= A[c * n + k] * b[k];
        b[c] /= A[c * n + c];
    }
    return b;
}

// Upper tail of the standard normal distribution.
double tail(double x) { return 0.5 * std::erfc(x / std::sqrt(2.0)); }

// Probability of a standard normal variable within [a, b], accurate in the tails.
double normal_mass(double a, double b) { return a >= 0 ? tail(a) - tail(b) : b <= 0 ? tail(-b) - tail(-a) : 1.0 - tail(-a) - tail(b); }

long floor_div(long a, long b) { return static_cast<long>(std::floor(double(a) / double(b))); }

// Adds weight w at value x to a histogram of bins lo + i step, split linearly between the
// two nearest bins (which keeps the mean).
void deposit(std::vector<double>& h, double lo, double step, double x, double w) {
    const double u = (x - lo) / step;
    const double fl = std::floor(u);
    const long i = static_cast<long>(fl);
    const double fr = u - fl;
    const long n = long(h.size());
    if (i + 1 <= 0) {
        h.front() += w;
    } else if (i >= n - 1) {
        h.back() += w;
    } else {
        h[static_cast<std::size_t>(i)] += w * (1 - fr);
        h[static_cast<std::size_t>(i + 1)] += w * fr;
    }
}

}  // namespace

EyeResult eye(const ChannelResult& ch, const EyeSpec& sp) {
    const auto t0 = Clock::now();
    EyeResult res;
    res.name = sp.name;
    res.channel = ch.name;
    const std::size_t np = ch.ports.size(), nf = ch.frequency_hz.size();
    if (nf < 2) throw std::runtime_error("eye " + sp.name + ": the channel has no sweep");
    const bool diff_tx = sp.tx_n != SIZE_MAX, diff_rx = sp.rx_n != SIZE_MAX;
    if (sp.tx_p >= np || sp.rx_p >= np || (diff_tx && sp.tx_n >= np) || (diff_rx && sp.rx_n >= np))
        throw std::runtime_error("eye " + sp.name + ": port beyond the channel's " + std::to_string(np));
    const double ui = 1.0 / sp.bit_rate;
    res.ui = ui;
    const double fmax = ch.frequency_hz.back(), df = ch.frequency_hz[1] - ch.frequency_hz[0];
    const double period = 1.0 / df;
    if (0.75 * fmax < sp.bit_rate) res.warnings.push_back("the sweep's untapered band (0.75 f_max) ends below the bit rate");
    const double rise = 1.25 * std::max(sp.rise_s, 0.0);
    if (sp.rise_s < 0.9 / fmax) res.warnings.push_back("the edges are faster than the sweep resolves (0.9 / f_max)");
    const double z0 = ch.impedance_ohm;

    // Transfer from a unit source (differential: +1/2 and -1/2 on the legs) to the
    // receiver's voltage (differential: the legs' difference), with the CTLE.
    std::vector<cplx> H(nf);
    double h0 = 0;  // the channel's DC transfer: its lowest frequency's, real (as in simulate())
    parallel_for(nf, [&](std::size_t k) {
        const double f = ch.frequency_hz[k], w = 2 * kPi * f;
        const auto& S = ch.s[k];
        // Z = Z0 (I + S)(I - S)^-1, by rows: z_i^T (I - S) = Z0 (e_i + S_i)^T.
        std::vector<cplx> At(np * np), Z(np * np);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j) At[i * np + j] = (i == j ? 1.0 : 0.0) - S[j * np + i];
        for (std::size_t i = 0; i < np; ++i) {
            std::vector<cplx> rhs(np);
            for (std::size_t j = 0; j < np; ++j) rhs[j] = z0 * ((i == j ? 1.0 : 0.0) + S[i * np + j]);
            const auto row = csolve(At, rhs, np);
            for (std::size_t j = 0; j < np; ++j) Z[i * np + j] = row[j];
        }
        // Terminations: admittance yt to the reference behind source vs; (I + Z Yt) V = Z Yt Vs.
        std::vector<cplx> yt(np, 0.0), vs(np, 0.0);
        auto load = [&](std::size_t p) { yt[p] = 1.0 / sp.r_term_ohm + cplx(0, w * sp.c_term_f); };
        load(sp.rx_p);
        if (diff_rx) load(sp.rx_n);
        yt[sp.tx_p] = 1.0 / sp.r_source_ohm;
        vs[sp.tx_p] = diff_tx ? 0.5 : 1.0;
        if (diff_tx) {
            yt[sp.tx_n] = 1.0 / sp.r_source_ohm;
            vs[sp.tx_n] = -0.5;
        }
        std::vector<cplx> M(np * np), b(np, 0.0);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j) {
                M[i * np + j] = (i == j ? 1.0 : 0.0) + Z[i * np + j] * yt[j];
                b[i] += Z[i * np + j] * yt[j] * vs[j];
            }
        const auto v = csolve(M, b, np);
        cplx h = v[sp.rx_p] - (diff_rx ? v[sp.rx_n] : 0.0);
        if (k == 0) h0 = h.real();
        if (sp.ctle.f_zero_hz > 0) {
            const cplx jf(0, f);
            cplx c = std::pow(10.0, sp.ctle.dc_gain_db / 20) * (1.0 + jf / sp.ctle.f_zero_hz);
            if (sp.ctle.f_pole1_hz > 0) c /= 1.0 + jf / sp.ctle.f_pole1_hz;
            if (sp.ctle.f_pole2_hz > 0) c /= 1.0 + jf / sp.ctle.f_pole2_hz;
            h *= c;
        }
        H[k] = h;
    });
    // At DC the CTLE is its DC gain.
    if (sp.ctle.f_zero_hz > 0) h0 *= std::pow(10.0, sp.ctle.dc_gain_db / 20);

    // The single-bit response over the sweep's period at UI / samples_per_ui; its first
    // tenth is the time before the bit (the band-limited edges spread to either side).
    const long spu = static_cast<long>(sp.samples_per_ui);
    const double dt = ui / double(spu);
    if (period < 8 * ui) throw std::runtime_error("eye " + sp.name + ": the sweep's period is under 8 UI; use more points");
    const std::size_t nper = static_cast<std::size_t>(std::floor(period / dt));
    const std::size_t npre = nper / 10;
    std::vector<double> pulse(nper, 0.0);  // index m: time (m - npre) dt; for a unit source swing
    {
        std::vector<cplx> hp(nf);
        for (std::size_t k = 0; k < nf; ++k) {
            const double f = ch.frequency_hz[k];
            const cplx x(0, 2 * kPi * f);
            // One bit: a ramp step up at 0 and down at UI.
            auto edge = [&](double ts) {
                if (rise <= 0) return cplx(std::exp(-x * ts) / x);
                return cplx((std::exp(-x * ts) - std::exp(-x * (ts + rise))) / (x * x * rise));
            };
            const double w = f <= 0.75 * fmax ? 1.0 : 0.5 * (1 + std::cos(kPi * (f - 0.75 * fmax) / (0.25 * fmax)));
            hp[k] = w * H[k] * (edge(0) - edge(ui));
        }
        parallel_for(nper, [&](std::size_t m) {
            const double t = (double(m) - double(npre)) * dt;
            double v = h0 * ui;  // DC: the bit's area
            for (std::size_t k = 0; k < nf; ++k) v += 2 * (hp[k] * std::exp(cplx(0, 2 * kPi * ch.frequency_hz[k] * t))).real();
            pulse[m] = v * df;
        });
    }
    // Transmitter FFE: tap i weights the response shifted by (i - main) UI.
    double ffe_sum = 1.0;
    if (!sp.ffe.empty()) {
        std::vector<double> eq(nper, 0.0);
        ffe_sum = 0;
        for (std::size_t i = 0; i < sp.ffe.size(); ++i) {
            ffe_sum += sp.ffe[i];
            const long shift = (long(i) - long(sp.ffe_main)) * spu;
            for (std::size_t m = 0; m < nper; ++m) {
                const long src = long(m) - shift;
                if (src >= 0 && src < long(nper)) eq[m] += sp.ffe[i] * pulse[static_cast<std::size_t>(src)];
            }
        }
        pulse.swap(eq);
    }
    // Symbols are +-amp about the middle level: half the swing, about the DC response to
    // half the swing for a single-ended driver (levels 0 and swing).
    const double amp = 0.5 * sp.swing_v;
    const double mid = diff_tx ? 0.0 : amp * h0 * ffe_sum;
    res.threshold_v = mid;
    std::size_t peak = 0;
    for (std::size_t m = 0; m < nper; ++m)
        if (std::abs(pulse[m]) > std::abs(pulse[peak])) peak = m;
    res.main_index = peak;
    res.pulse_t.resize(nper);
    res.pulse_v.resize(nper);
    for (std::size_t m = 0; m < nper; ++m) {
        res.pulse_t[m] = (double(m) - double(npre)) * dt;
        res.pulse_v[m] = sp.swing_v * pulse[m];
    }
    const double pmax = std::abs(pulse[peak]);
    if (!(pmax > 0)) throw std::runtime_error("eye " + sp.name + ": no response at the receiver");
    // The span of the response: from the first to the last sample above 1e-4 of its peak.
    std::size_t first = 0, last = nper - 1;
    while (first < peak && std::abs(pulse[first]) < 1e-4 * pmax) ++first;
    while (last > peak && std::abs(pulse[last]) < 1e-4 * pmax) --last;
    if (first == 0 || last + 1 == nper)
        res.warnings.push_back("the single-bit response lasts nearly the sweep's period; a finer sweep is needed");

    // Cursors of a sampling time (pulse sample n): c_k = p(n + k spu) for the bits k UI
    // before (k > 0) and after (k < 0) the sampled one, within the response's span; k0 is
    // the index of the main cursor.
    auto cursors_at = [&](long n, std::vector<double>& c, long& k0) {
        const long kmin = std::min(0L, -floor_div(n - long(first), spu)), kmax = std::max(0L, floor_div(long(last) - n, spu));
        c.clear();
        for (long k = kmin; k <= kmax; ++k) {
            const long m = n + k * spu;
            c.push_back(m >= 0 && m < long(nper) ? pulse[static_cast<std::size_t>(m)] : 0.0);
        }
        k0 = -kmin;
    };
    // DFE: taps cancelling the post-cursors at the peak, each limited.
    std::vector<double> dfe(sp.dfe_taps, 0.0);
    {
        std::vector<double> c;
        long k0;
        cursors_at(long(peak), c, k0);
        const double lim = sp.dfe_limit * std::abs(c[static_cast<std::size_t>(k0)]);
        for (std::size_t k = 0; k < sp.dfe_taps; ++k) {
            const std::size_t idx = static_cast<std::size_t>(k0) + 1 + k;
            if (idx < c.size()) dfe[k] = std::clamp(c[idx], -lim, lim);
        }
        for (double d : dfe) res.dfe.push_back(sp.swing_v * d);
    }
    auto equalised = [&](long n, std::vector<double>& c, long& k0) {
        cursors_at(n, c, k0);
        for (std::size_t k = 0; k < dfe.size(); ++k) {
            const std::size_t idx = static_cast<std::size_t>(k0) + 1 + k;
            if (idx < c.size()) c[idx] -= dfe[k];
        }
    };

    // Statistical eye over two UI of phases centred on the peak, which hold the eye's
    // opening wherever it lies relative to the peak.
    const long nph = 2 * spu;
    const long n_start = long(peak) - spu;
    // Jitter: sampling offsets (samples) and their probabilities; Gaussian RJ to 9 sigma,
    // each offset taking the mass of its sample interval, and dual-Dirac DJ.
    std::vector<std::pair<long, double>> jit;
    {
        std::map<long, double> acc;
        const std::vector<double> djs = sp.dj_pp_s > 0 ? std::vector<double>{-0.5 * sp.dj_pp_s, 0.5 * sp.dj_pp_s} : std::vector<double>{0.0};
        for (const double dj : djs) {
            const double share = 1.0 / double(djs.size());
            if (sp.rj_rms_s > 0) {
                const long r = static_cast<long>(std::ceil(9 * sp.rj_rms_s / dt));
                const long c = std::lround(dj / dt);
                for (long o = -r; o <= r; ++o)
                    acc[c + o] += share * normal_mass((double(o) - 0.5) * dt / sp.rj_rms_s, (double(o) + 0.5) * dt / sp.rj_rms_s);
            } else {
                acc[std::lround(dj / dt)] += share;
            }
        }
        for (const auto& [o, w] : acc)
            if (w > 1e-3 * sp.ber) jit.push_back({o, w});
    }
    // The ISI of a sampling time: the main cursor and a histogram of the sum of the other
    // cursors with equiprobable signs, on bins lo + i step.
    struct Isi {
        double main = 0, worst = 0, lo = 0, step = 1;
        std::vector<double> w;
    };
    const std::size_t nb = 2048;
    auto isi_at = [&](long n) {
        std::vector<double> c;
        long k0;
        equalised(n, c, k0);
        Isi d;
        d.main = amp * c[static_cast<std::size_t>(k0)];
        for (std::size_t k = 0; k < c.size(); ++k)
            if (long(k) != k0) d.worst += amp * std::abs(c[k]);
        if (!(d.worst > 0)) {
            d.w = {1.0};
            return d;
        }
        d.lo = -d.worst;
        d.step = 2 * d.worst / double(nb - 1);
        d.w.assign(nb, 0.0);
        deposit(d.w, d.lo, d.step, 0.0, 1.0);
        std::vector<double> nw(nb);
        for (std::size_t k = 0; k < c.size(); ++k) {
            const double x = amp * c[k];
            if (long(k) == k0 || x == 0) continue;
            std::fill(nw.begin(), nw.end(), 0.0);
            for (std::size_t i = 0; i < nb; ++i) {
                if (d.w[i] == 0) continue;
                const double v = d.lo + double(i) * d.step;
                deposit(nw, d.lo, d.step, v + x, 0.5 * d.w[i]);
                deposit(nw, d.lo, d.step, v - x, 0.5 * d.w[i]);
            }
            d.w.swap(nw);
        }
        return d;
    };
    std::vector<long> needed;
    for (long j = 0; j < nph; ++j)
        for (const auto& [off, wgt] : jit) needed.push_back(n_start + j + off);
    std::sort(needed.begin(), needed.end());
    needed.erase(std::unique(needed.begin(), needed.end()), needed.end());
    std::vector<Isi> isis(needed.size());
    parallel_for(needed.size(), [&](std::size_t i) { isis[i] = isi_at(needed[i]); });
    auto isi = [&](long n) -> const Isi& { return isis[static_cast<std::size_t>(std::lower_bound(needed.begin(), needed.end(), n) - needed.begin())]; };

    const double sn = sp.noise_rms_v;
    res.phase_t.resize(static_cast<std::size_t>(nph));
    res.upper_v.assign(static_cast<std::size_t>(nph), 0.0);
    res.lower_v.assign(static_cast<std::size_t>(nph), 0.0);
    std::vector<double> pda(static_cast<std::size_t>(nph), 0.0);
    parallel_for(static_cast<std::size_t>(nph), [&](std::size_t ju) {
        const long j = long(ju);
        res.phase_t[ju] = double(j - spu) * dt;
        const Isi& d0 = isi(n_start + j);
        pda[ju] = 2 * (d0.main - d0.worst);
        // The distribution of a one's level, a mixture over the jitter's offsets, on a grid
        // with noise convolved; the inner edge of the ones is where its lower tail reaches
        // the BER. Zeros mirror the ones about the middle.
        double vlo = INFINITY, vhi = -INFINITY;
        for (const auto& [off, wgt] : jit) {
            const Isi& d = isi(n_start + j + off);
            vlo = std::min(vlo, d.main - d.worst);
            vhi = std::max(vhi, d.main + d.worst);
        }
        vlo -= 10 * sn + 1e-9;
        vhi += 10 * sn + 1e-9;
        const std::size_t ng = 8192;
        const double g = (vhi - vlo) / double(ng - 1);
        std::vector<double> dist(ng, 0.0);
        for (const auto& [off, wgt] : jit) {
            const Isi& d = isi(n_start + j + off);
            for (std::size_t i = 0; i < d.w.size(); ++i)
                if (d.w[i] > 0) deposit(dist, vlo, g, d.main + d.lo + double(i) * d.step, wgt * d.w[i]);
        }
        double one = vlo;
        if (sn > 0) {
            // P(one below x) = sum_i dist_i Phi((x - v_i) / sigma), by bisection.
            std::vector<std::pair<double, double>> pts;
            for (std::size_t i = 0; i < ng; ++i)
                if (dist[i] > 0) pts.push_back({vlo + double(i) * g, dist[i]});
            auto below = [&](double x) {
                double p = 0;
                for (const auto& [v, w] : pts) p += w * tail((v - x) / sn);
                return p;
            };
            double a = vlo, b = vhi;
            for (int it = 0; it < 60; ++it) {
                const double m = 0.5 * (a + b);
                (below(m) <= sp.ber ? a : b) = m;
            }
            one = a;
        } else {
            // The cumulative probability from below, each bin's mass spread over [v - g/2, v + g/2].
            double cum = 0;
            for (std::size_t i = 0; i < ng; ++i) {
                if (cum + dist[i] > sp.ber) {
                    one = vlo + (double(i) - 0.5 + (sp.ber - cum) / dist[i]) * g;
                    break;
                }
                cum += dist[i];
                one = vlo + (double(i) + 0.5) * g;
            }
        }
        res.lower_v[ju] = mid + one;
        res.upper_v[ju] = mid - one;
    });
    // Height and width of the opening.
    std::size_t best = 0;
    for (std::size_t j = 0; j < res.lower_v.size(); ++j)
        if (res.lower_v[j] - res.upper_v[j] > res.lower_v[best] - res.upper_v[best]) best = j;
    res.eye_height_v = std::max(0.0, res.lower_v[best] - res.upper_v[best]);
    res.best_phase_s = res.phase_t[best];
    if (res.eye_height_v > 0) {
        std::size_t a = best, b = best;
        while (a > 0 && res.lower_v[a - 1] > res.upper_v[a - 1]) --a;
        while (b + 1 < res.lower_v.size() && res.lower_v[b + 1] > res.upper_v[b + 1]) ++b;
        res.eye_width_ui = std::min(1.0, double(b - a + 1) / double(spu));
    }
    res.peak_distortion_v = std::max(0.0, *std::max_element(pda.begin(), pda.end()));
    {
        std::vector<double> c;
        long k0;
        equalised(n_start + long(best), c, k0);
        res.pre_cursors = static_cast<std::size_t>(k0);
        for (double v : c) res.cursors.push_back(sp.swing_v * v);
    }

    // Eye of the pattern: superposition of the single-bit response, folded over two UI
    // centred on the best phase, after the response's span has filled.
    {
        Stimulus st;
        st.pattern = sp.pattern;
        st.bits = sp.bits;
        std::vector<int> seq = st.sequence();
        while (seq.size() < sp.bits) seq.insert(seq.end(), seq.begin(), seq.end());
        seq.resize(sp.bits);
        const long span_bits = (long(last) - long(first)) / spu + 2;
        const long nsamp = long(seq.size()) * spu;
        std::vector<double> wave(static_cast<std::size_t>(nsamp), mid);
        for (std::size_t k = 0; k < seq.size(); ++k) {
            const double s = seq[k] ? amp : -amp;
            const long base = long(k) * spu - long(npre);
            for (long m = long(first); m <= long(last); ++m) {
                const long n = base + m;
                if (n >= 0 && n < nsamp) wave[static_cast<std::size_t>(n)] += s * pulse[static_cast<std::size_t>(m)];
            }
        }
        // The sample of a bit's best phase, from the bit's start.
        const long offset = long(peak) - long(npre) - spu + long(best);
        const std::size_t nv = 128, nt2 = static_cast<std::size_t>(2 * spu);
        double vmin = INFINITY, vmax = -INFINITY;
        std::vector<long> centres;
        for (long k = span_bits; k + span_bits < long(seq.size()); ++k) {
            const long centre = k * spu + offset;
            if (centre - spu < 0 || centre + spu >= nsamp) continue;
            centres.push_back(centre);
            for (long i = -spu; i < spu; ++i) {
                vmin = std::min(vmin, wave[static_cast<std::size_t>(centre + i)]);
                vmax = std::max(vmax, wave[static_cast<std::size_t>(centre + i)]);
            }
        }
        if (centres.empty()) {
            res.warnings.push_back("the pattern is too short for its eye: it needs more bits than twice the response's span");
        } else {
            const double margin = 0.05 * (vmax - vmin);
            vmin -= margin;
            vmax += margin;
            res.eye_v.resize(nv);
            for (std::size_t i = 0; i < nv; ++i) res.eye_v[i] = vmin + (double(i) + 0.5) * (vmax - vmin) / double(nv);
            res.eye_t.resize(nt2);
            for (std::size_t i = 0; i < nt2; ++i) res.eye_t[i] = (double(i) - double(spu)) * dt;
            res.density.assign(nt2, std::vector<double>(nv, 0.0));
            // Each trace marks the bins it passes in a column: from its sample to the
            // midpoint towards each neighbour.
            auto bin = [&](double v) { return std::clamp(static_cast<long>(std::floor((v - vmin) / (vmax - vmin) * double(nv))), 0L, long(nv) - 1); };
            for (const long centre : centres)
                for (long i = 0; i < long(nt2); ++i) {
                    const long n = centre - spu + i;
                    const double v = wave[static_cast<std::size_t>(n)];
                    const double prev = 0.5 * (v + wave[static_cast<std::size_t>(std::max(n - 1, 0L))]);
                    const double next = 0.5 * (v + wave[static_cast<std::size_t>(std::min(n + 1, nsamp - 1))]);
                    const long lo = bin(std::min({v, prev, next})), hi = bin(std::max({v, prev, next}));
                    for (long b = lo; b <= hi; ++b) res.density[static_cast<std::size_t>(i)][static_cast<std::size_t>(b)] += 1.0 / double(centres.size());
                }
        }
    }
    res.seconds = since(t0);
    return res;
}

std::string eyes_summary(const std::vector<EyeResult>& rs) {
    std::ostringstream o;
    char line[256];
    for (const auto& r : rs) {
        std::snprintf(line, sizeof line, "eye %s (channel %s), %.4g Gbit/s: height %.4f V and width %.3f UI at the BER; worst case %.4f V; %.1f s\n",
                      r.name.c_str(), r.channel.c_str(), 1e-9 / r.ui, r.eye_height_v, r.eye_width_ui, r.peak_distortion_v, r.seconds);
        o << line;
        std::snprintf(line, sizeof line, "  main cursor %.4f V, %zu pre-cursors, %zu post-cursors; best phase %+.1f ps from the peak",
                      r.pre_cursors < r.cursors.size() ? r.cursors[r.pre_cursors] : 0.0, r.pre_cursors,
                      r.cursors.empty() ? std::size_t(0) : r.cursors.size() - r.pre_cursors - 1, r.best_phase_s * 1e12);
        o << line;
        if (!r.dfe.empty()) {
            o << "; DFE taps";
            for (double d : r.dfe) {
                std::snprintf(line, sizeof line, " %.4f", d);
                o << line;
            }
            o << " V";
        }
        o << "\n";
        for (const auto& w : r.warnings) o << "  warning: " << w << "\n";
    }
    return o.str();
}

std::string eyes_json(const std::vector<EyeResult>& rs) {
    using nlohmann::json;
    json out = json::array();
    for (const auto& r : rs)
        out.push_back({{"name", r.name}, {"channel", r.channel}, {"ui_s", r.ui}, {"pulse_t", r.pulse_t}, {"pulse_v", r.pulse_v},
                       {"main_index", r.main_index}, {"cursors", r.cursors}, {"pre_cursors", r.pre_cursors}, {"dfe", r.dfe},
                       {"eye_t", r.eye_t}, {"eye_v", r.eye_v}, {"density", r.density}, {"phase_t", r.phase_t}, {"upper_v", r.upper_v},
                       {"lower_v", r.lower_v}, {"eye_height_v", r.eye_height_v}, {"eye_width_ui", r.eye_width_ui},
                       {"best_phase_s", r.best_phase_s}, {"peak_distortion_v", r.peak_distortion_v}, {"threshold_v", r.threshold_v},
                       {"seconds", r.seconds}, {"warnings", r.warnings}});
    return out.dump();
}

std::string eye_svg(const EyeResult& r) {
    // The pattern's density (darker where more traces pass) and the statistical eye's
    // inner edges at the BER, over the central UI.
    const std::size_t nt = r.eye_t.size(), nv = r.eye_v.size();
    if (nt < 2 || nv < 2) return "";
    const double W = 640, Hh = 400, L = 64, R = 12, T = 24, B = 36;
    const double cw = (W - L - R) / double(nt), chh = (Hh - T - B) / double(nv);
    const double dv = r.eye_v[1] - r.eye_v[0];
    const double t0 = r.eye_t.front(), t1 = r.eye_t.back() + (r.eye_t[1] - r.eye_t[0]);
    const double v0 = r.eye_v.front() - 0.5 * dv, v1 = r.eye_v.back() + 0.5 * dv;
    auto X = [&](double t) { return L + (W - L - R) * (t - t0) / (t1 - t0); };
    auto Y = [&](double v) { return Hh - B - (Hh - T - B) * (v - v0) / (v1 - v0); };
    double dmax = 0;
    for (const auto& c : r.density)
        for (double x : c) dmax = std::max(dmax, x);
    std::ostringstream o;
    char buf[256];
    std::snprintf(buf, sizeof buf, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" font-family=\"sans-serif\" font-size=\"11\">\n", W, Hh);
    o << buf << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    for (std::size_t i = 0; i < nt; ++i)
        for (std::size_t k = 0; k < nv; ++k) {
            const double d = r.density[i][k];
            if (d <= 0) continue;
            const double a = 0.15 + 0.85 * std::log10(1 + 9 * d / dmax);
            std::snprintf(buf, sizeof buf, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.2f\" height=\"%.2f\" fill=\"#1f4e9a\" fill-opacity=\"%.3f\"/>\n",
                          L + double(i) * cw, Hh - B - double(k + 1) * chh, cw + 0.3, chh + 0.3, a);
            o << buf;
        }
    for (const auto* edge : {&r.lower_v, &r.upper_v}) {
        o << "<polyline fill=\"none\" stroke=\"#d62728\" stroke-width=\"1.5\" points=\"";
        for (std::size_t j = 0; j < r.phase_t.size(); ++j) {
            const double t = r.phase_t[j] - r.best_phase_s;
            if (t < t0 || t > t1) continue;
            const double v = std::clamp((*edge)[j], v0, v1);
            std::snprintf(buf, sizeof buf, "%.1f,%.1f ", X(t), Y(v));
            o << buf;
        }
        o << "\"/>\n";
    }
    std::snprintf(buf, sizeof buf, "<text x=\"%.0f\" y=\"16\">%s: height %.3f V, width %.2f UI at the BER</text>\n", L, r.name.c_str(),
                  r.eye_height_v, r.eye_width_ui);
    o << buf;
    std::snprintf(buf, sizeof buf, "<text x=\"4\" y=\"%.0f\">%.3f V</text>\n<text x=\"4\" y=\"%.0f\">%.3f V</text>\n", Y(v1) + 10, v1, Y(v0), v0);
    o << buf;
    std::snprintf(buf, sizeof buf, "<text x=\"%.0f\" y=\"%.0f\">%.1f ps</text>\n<text x=\"%.0f\" y=\"%.0f\" text-anchor=\"end\">%.1f ps</text>\n", L,
                  Hh - 14, t0 * 1e12, W - R, Hh - 14, t1 * 1e12);
    o << buf << "</svg>\n";
    return o.str();
}

}  // namespace kairoq
