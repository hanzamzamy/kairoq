// IBIS buffer models: the reader of the keywords that describe a model's behaviour at its
// pin, and its switching coefficients (docs/solvers/timedomain.md).
#include "kairoq/ibis.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r") - a + 1);
}

// Keyword names compare without regard to case, spaces or underscores.
std::string norm(std::string s) {
    std::string out;
    for (char c : lower(s))
        if (c != ' ' && c != '_') out += c;
    return out;
}

// A number with an optional scale (T G M k m u n p f) and units after it; "NA" is missing.
std::optional<double> number(const std::string& word) {
    if (word.empty() || lower(word) == "na") return std::nullopt;
    std::size_t used = 0;
    double v;
    try {
        v = std::stod(word, &used);
    } catch (...) {
        return std::nullopt;
    }
    if (used < word.size()) {
        switch (word[used]) {
        case 'T': v *= 1e12; break;
        case 'G': v *= 1e9; break;
        case 'M': v *= 1e6; break;
        case 'k': case 'K': v *= 1e3; break;
        case 'm': v *= 1e-3; break;
        case 'u': v *= 1e-6; break;
        case 'n': v *= 1e-9; break;
        case 'p': v *= 1e-12; break;
        case 'f': v *= 1e-15; break;
        default: break;
        }
    }
    return v;
}

std::vector<std::string> words(const std::string& s) {
    std::istringstream in(s);
    std::vector<std::string> out;
    std::string w;
    while (in >> w) out.push_back(w);
    return out;
}

// Typical, minimum and maximum from up to three words; missing ones take the typical.
IbisValue triple(const std::vector<std::string>& w, std::size_t from) {
    IbisValue v;
    const auto t = from < w.size() ? number(w[from]) : std::nullopt;
    if (!t) throw std::runtime_error("IBIS: no typical value");
    v.typ = *t;
    const auto mn = from + 1 < w.size() ? number(w[from + 1]) : std::nullopt;
    const auto mx = from + 2 < w.size() ? number(w[from + 2]) : std::nullopt;
    v.min = mn.value_or(v.typ);
    v.max = mx.value_or(v.typ);
    return v;
}

// dV/dt values are given as a fraction "dv/dt".
double ratio(const std::string& w) {
    const auto slash = w.find('/');
    if (slash == std::string::npos) return number(w).value_or(0.0);
    const auto a = number(w.substr(0, slash)), b = number(w.substr(slash + 1));
    return a && b && *b != 0 ? *a / *b : 0.0;
}

}  // namespace

double IbisTable::at(double v, int corner) const {
    const auto& ys = y[static_cast<std::size_t>(std::clamp(corner, 0, 2))];
    if (x.empty()) return 0.0;
    if (x.size() == 1) return ys[0];
    std::size_t k = 1;
    if (v > x.front()) {
        k = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), v) - x.begin());
        k = std::clamp<std::size_t>(k, 1, x.size() - 1);
    }
    const double w = (v - x[k - 1]) / (x[k] - x[k - 1]);
    return ys[k - 1] + w * (ys[k] - ys[k - 1]);
}

double IbisModel::vpu(int c) const {
    if (pullup_ref) return pullup_ref->at(c);
    return voltage_range.at(c);
}

double IbisModel::vpd(int c) const { return pulldown_ref ? pulldown_ref->at(c) : 0.0; }

double IbisModel::clamp_current(double v, int c) const {
    const double vpc = power_clamp_ref ? power_clamp_ref->at(c) : vpu(c);
    const double vgc = gnd_clamp_ref ? gnd_clamp_ref->at(c) : 0.0;
    return gnd_clamp.at(v - vgc, c) + power_clamp.at(vpc - v, c);
}

double IbisModel::current(double v, double ku, double kd, int c) const {
    return ku * pullup.at(vpu(c) - v, c) + kd * pulldown.at(v - vpd(c), c) + clamp_current(v, c);
}

const IbisModel* IbisFile::model(const std::string& name) const {
    for (const auto& m : models)
        if (lower(m.name) == lower(name)) return &m;
    return nullptr;
}

IbisFile parse_ibis(const std::string& text) {
    IbisFile f;
    char comment = '|';
    std::istringstream in(text);
    std::string line, keyword;
    IbisModel* model = nullptr;
    IbisTable* table = nullptr;
    IbisWaveform* wave = nullptr;
    std::set<std::string> ignored;
    auto add_row = [&](const std::vector<std::string>& w) {
        if (!table || w.size() < 2) return;
        const auto x = number(w[0]), t = number(w[1]);
        if (!x || !t) return;
        table->x.push_back(*x);
        table->y[0].push_back(*t);
        table->y[1].push_back(w.size() > 2 ? number(w[2]).value_or(*t) : *t);
        table->y[2].push_back(w.size() > 3 ? number(w[3]).value_or(*t) : *t);
    };
    while (std::getline(in, line)) {
        if (const auto c = line.find(comment); c != std::string::npos) line = line.substr(0, c);
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '[') {
            const auto close = line.find(']');
            if (close == std::string::npos) continue;
            keyword = norm(line.substr(1, close - 1));
            const std::string rest = trim(line.substr(close + 1));
            table = nullptr;
            if (keyword == "commentchar") {
                if (!rest.empty()) comment = rest[0];
            } else if (keyword == "ibisver") {
                f.version = rest;
            } else if (keyword == "component") {
                f.component = rest;
            } else if (keyword == "model") {
                f.models.push_back({});
                model = &f.models.back();
                model->name = rest;
                wave = nullptr;
            } else if (keyword == "end") {
                break;
            } else if (model) {
                const auto w = words(rest);
                if (keyword == "voltagerange") model->voltage_range = triple(w, 0);
                else if (keyword == "pullupreference") model->pullup_ref = triple(w, 0);
                else if (keyword == "pulldownreference") model->pulldown_ref = triple(w, 0);
                else if (keyword == "powerclampreference") model->power_clamp_ref = triple(w, 0);
                else if (keyword == "gndclampreference") model->gnd_clamp_ref = triple(w, 0);
                else if (keyword == "pullup") table = &model->pullup;
                else if (keyword == "pulldown") table = &model->pulldown;
                else if (keyword == "powerclamp") table = &model->power_clamp;
                else if (keyword == "gndclamp") table = &model->gnd_clamp;
                else if (keyword == "ramp") wave = nullptr;
                else if (keyword == "risingwaveform" || keyword == "fallingwaveform") {
                    auto& list = keyword == "risingwaveform" ? model->rising : model->falling;
                    list.push_back({});
                    wave = &list.back();
                    table = &wave->v;
                } else {
                    ignored.insert(keyword);
                }
            } else {
                ignored.insert(keyword);
            }
            continue;
        }
        // Lines within a keyword.
        const auto w = words(line);
        if (keyword == "package" && w.size() >= 2) {
            const std::string k = lower(w[0]);
            if (k == "r_pkg") f.package.r = triple(w, 1);
            else if (k == "l_pkg") f.package.l = triple(w, 1);
            else if (k == "c_pkg") f.package.c = triple(w, 1);
        } else if (keyword == "pin" && w.size() >= 3) {
            if (lower(w[0]) == "signal_name") continue;
            IbisPin p{w[0], w[1], w[2], std::nullopt, std::nullopt, std::nullopt};
            if (w.size() >= 6) {
                p.r = number(w[3]);
                p.l = number(w[4]);
                p.c = number(w[5]);
            }
            f.pins.push_back(p);
        } else if (keyword == "model" && model) {
            // Subparameters: "name value...", or "name = value".
            std::string k = lower(w[0]);
            std::vector<std::string> v(w.begin() + 1, w.end());
            if (!v.empty() && v[0] == "=") v.erase(v.begin());
            if (k == "model_type" && !v.empty()) model->type = v[0];
            else if (k == "c_comp" && !v.empty()) model->c_comp = triple(v, 0);
            else if (k == "vinl" && !v.empty()) model->vinl = triple(v, 0);
            else if (k == "vinh" && !v.empty()) model->vinh = triple(v, 0);
        } else if (keyword == "ramp" && model && !w.empty()) {
            const std::string k = lower(w[0]);
            if (k == "dv/dt_r" && w.size() >= 2) {
                model->ramp_rise.typ = ratio(w[1]);
                model->ramp_rise.min = w.size() > 2 ? ratio(w[2]) : model->ramp_rise.typ;
                model->ramp_rise.max = w.size() > 3 ? ratio(w[3]) : model->ramp_rise.typ;
            } else if (k == "dv/dt_f" && w.size() >= 2) {
                model->ramp_fall.typ = ratio(w[1]);
                model->ramp_fall.min = w.size() > 2 ? ratio(w[2]) : model->ramp_fall.typ;
                model->ramp_fall.max = w.size() > 3 ? ratio(w[3]) : model->ramp_fall.typ;
            } else if (k == "r_load") {
                const std::string v = w.size() >= 3 && w[1] == "=" ? w[2] : (w.size() >= 2 ? w[1] : "");
                model->ramp_r_load = number(v).value_or(50.0);
            }
        } else if (wave && (keyword == "risingwaveform" || keyword == "fallingwaveform")) {
            if (w.size() >= 2 && (w[1] == "=" || !number(w[0]))) {
                const std::string k = lower(w[0]);
                const auto v = number(w[1] == "=" && w.size() >= 3 ? w[2] : w[1]);
                if (!v) continue;
                if (k == "r_fixture") wave->r_fixture = *v;
                else if (k == "v_fixture") wave->v_fixture = *v;
                else if (k == "c_fixture") wave->c_fixture = *v;
                else if (k == "l_fixture") wave->l_fixture = *v;
                else if (k == "r_dut") wave->r_dut = *v;
                else if (k == "l_dut") wave->l_dut = *v;
                else if (k == "c_dut") wave->c_dut = *v;
                continue;
            }
            add_row(w);
        } else if (table) {
            add_row(w);
        }
    }
    for (const auto& k : ignored) f.warnings.push_back("IBIS keyword [" + k + "] is not read");
    for (auto& m : f.models)
        for (IbisTable* t : {&m.pullup, &m.pulldown, &m.power_clamp, &m.gnd_clamp}) {
            // Tables in ascending voltage.
            std::vector<std::size_t> order(t->x.size());
            for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
            std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return t->x[a] < t->x[b]; });
            IbisTable s;
            for (auto i : order) {
                s.x.push_back(t->x[i]);
                for (int c = 0; c < 3; ++c) s.y[static_cast<std::size_t>(c)].push_back(t->y[static_cast<std::size_t>(c)][i]);
            }
            *t = std::move(s);
        }
    return f;
}

IbisFile read_ibis(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return parse_ibis(ss.str());
}

IbisSwitching switching_coefficients(const IbisModel& m, bool rising, int corner, std::vector<std::string>* warnings) {
    IbisSwitching out;
    const auto& waves = rising ? m.rising : m.falling;
    const double cc = m.c_comp.at(corner);
    // The current the pin delivers into a fixture at time t, with the die capacitance.
    auto fixture_balance = [&](const IbisWaveform& w, double t, double& v, double& dvdt) {
        v = w.v.at(t, corner);
        const double h = 1e-13;
        dvdt = (w.v.at(t + h, corner) - w.v.at(t - h, corner)) / (2 * h);
        // Into the pin: minus the fixture's current, minus the die and fixture capacitance.
        return -(v - w.v_fixture) / w.r_fixture - (w.c_fixture + cc) * dvdt;
    };
    for (const auto& w : waves)
        if ((w.l_fixture > 0 || w.r_dut > 0 || w.l_dut > 0 || w.c_dut > 0) && warnings)
            warnings->push_back("model " + m.name + ": fixture inductance and DUT parasitics of the waveforms are not used");
    if (waves.size() >= 2) {
        // Two waveforms (into fixtures to different voltages): ku, kd from two equations.
        const auto& a = waves[0];
        const auto& b = waves[1];
        std::set<double> times(a.v.x.begin(), a.v.x.end());
        times.insert(b.v.x.begin(), b.v.x.end());
        for (double t : times) {
            double va, da, vb, db;
            const double ia = fixture_balance(a, t, va, da), ib = fixture_balance(b, t, vb, db);
            const double pa = m.pullup.at(m.vpu(corner) - va, corner), qa = m.pulldown.at(va - m.vpd(corner), corner);
            const double pb = m.pullup.at(m.vpu(corner) - vb, corner), qb = m.pulldown.at(vb - m.vpd(corner), corner);
            const double ra = ia - m.clamp_current(va, corner), rb = ib - m.clamp_current(vb, corner);
            const double det = pa * qb - pb * qa;
            if (std::abs(det) < 1e-30) continue;
            out.time_s.push_back(t);
            out.ku.push_back((ra * qb - rb * qa) / det);
            out.kd.push_back((pa * rb - pb * ra) / det);
        }
        return out;
    }
    if (waves.size() == 1) {
        // One waveform: the stages switch complementarily, kd = 1 - ku.
        const auto& a = waves[0];
        for (double t : a.v.x) {
            double v, d;
            const double i = fixture_balance(a, t, v, d) - m.clamp_current(v, corner);
            const double p = m.pullup.at(m.vpu(corner) - v, corner), q = m.pulldown.at(v - m.vpd(corner), corner);
            if (std::abs(p - q) < 1e-30) continue;
            const double ku = (i - q) / (p - q);
            out.time_s.push_back(t);
            out.ku.push_back(ku);
            out.kd.push_back(1 - ku);
        }
        return out;
    }
    // The ramp: its slope is that of 20 to 80 % of the swing; a linear transition over the
    // full swing at that slope.
    const double slope = (rising ? m.ramp_rise : m.ramp_fall).at(corner);
    if (!(slope > 0)) throw std::runtime_error("IBIS model " + m.name + " has neither waveforms nor a ramp");
    const double t_full = (m.vpu(corner) - m.vpd(corner)) / slope;
    if (warnings) warnings->push_back("model " + m.name + ": switching from the ramp (no V-t waveforms)");
    for (int k = 0; k <= 20; ++k) {
        const double t = t_full * k / 20.0, s = k / 20.0;
        out.time_s.push_back(t);
        out.ku.push_back(rising ? s : 1 - s);
        out.kd.push_back(rising ? 1 - s : s);
    }
    return out;
}

}  // namespace kairoq
