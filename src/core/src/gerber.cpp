#include "kairoq/gerber.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace kairoq {

using namespace Clipper2Lib;

namespace {

std::string trim(std::string_view s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return std::string(s.substr(a, b - a));
}

std::string strip_newlines(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s)
        if (c != '\r' && c != '\n') out.push_back(c);
    return out;
}

std::vector<std::string> split(std::string_view s, char sep) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == sep) {
            out.emplace_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    return out;
}

double to_double(std::string_view s) {
    const std::string t = trim(s);
    if (t.empty()) return 0.0;
    std::size_t used = 0;
    double v = std::stod(t, &used);
    return v;
}

// ---------------------------------------------------------------- macro expressions

class ExpressionParser {
public:
    ExpressionParser(std::string_view text, const std::vector<double>& vars) : s_(text), vars_(vars) {}

    double parse() {
        const double v = sum();
        skip();
        if (i_ != s_.size()) throw std::runtime_error("unexpected text in macro expression");
        return v;
    }

private:
    void skip() {
        while (i_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[i_]))) ++i_;
    }
    bool eat(char c) {
        skip();
        if (i_ < s_.size() && s_[i_] == c) {
            ++i_;
            return true;
        }
        return false;
    }
    double sum() {
        double v = product();
        while (true) {
            if (eat('+')) v += product();
            else if (eat('-')) v -= product();
            else return v;
        }
    }
    double product() {
        double v = unary();
        while (true) {
            if (eat('x') || eat('X')) v *= unary();
            else if (eat('/')) v /= unary();
            else return v;
        }
    }
    double unary() {
        if (eat('-')) return -unary();
        if (eat('+')) return unary();
        return primary();
    }
    double primary() {
        skip();
        if (eat('(')) {
            const double v = sum();
            if (!eat(')')) throw std::runtime_error("missing ')' in macro expression");
            return v;
        }
        if (i_ < s_.size() && s_[i_] == '$') {
            ++i_;
            std::size_t j = i_;
            while (j < s_.size() && std::isdigit(static_cast<unsigned char>(s_[j]))) ++j;
            const int n = std::stoi(std::string(s_.substr(i_, j - i_)));
            i_ = j;
            return n < static_cast<int>(vars_.size()) ? vars_[n] : 0.0;
        }
        std::size_t j = i_;
        while (j < s_.size() && (std::isdigit(static_cast<unsigned char>(s_[j])) || s_[j] == '.')) ++j;
        if (j == i_) throw std::runtime_error("expected a number in macro expression");
        const double v = std::stod(std::string(s_.substr(i_, j - i_)));
        i_ = j;
        return v;
    }

    std::string_view s_;
    const std::vector<double>& vars_;
    std::size_t i_ = 0;
};

// ---------------------------------------------------------------- apertures

struct MacroStatement {
    bool assignment = false;
    int variable = 0;                // assignment target
    std::vector<std::string> exprs;  // primitive code and modifiers, or the assigned expression
};

struct Macro {
    std::vector<MacroStatement> statements;
};

enum class ApertureType { Circle, Rect, Obround, Polygon, Macro, Block };

struct Aperture {
    ApertureType type = ApertureType::Circle;
    Paths64 shape;                    // at the origin, before load transformations
    double circle_diameter_nm = -1;  // circle apertures only
    std::string function;             // .AperFunction, first field
    std::vector<Feature> block;       // aperture blocks
};

Path64 rect_path(double cx, double cy, double w, double h) {
    const double x0 = cx - w / 2, x1 = cx + w / 2, y0 = cy - h / 2, y1 = cy + h / 2;
    return {Point64(std::llround(x0), std::llround(y0)), Point64(std::llround(x1), std::llround(y0)),
            Point64(std::llround(x1), std::llround(y1)), Point64(std::llround(x0), std::llround(y1))};
}

Path64 regular_polygon(double cx, double cy, double diameter, int n, double rotation_deg) {
    Path64 path;
    const double r = diameter / 2;
    for (int i = 0; i < n; ++i) {
        const double a = (rotation_deg * kPi / 180.0) + 2.0 * kPi * i / n;
        path.emplace_back(std::llround(cx + r * std::cos(a)), std::llround(cy + r * std::sin(a)));
    }
    return path;
}

Paths64 rotate_about_origin(const Paths64& paths, double deg) {
    if (deg == 0.0) return paths;
    Placement p;
    p.rotation_deg = deg;
    return transform(paths, p);
}

// ---------------------------------------------------------------- parser

class Parser {
public:
    explicit Parser(const GeometryOptions& options) : tol_(options.arc_tolerance_nm) {}

    GerberImage run(std::string_view text) {
        std::size_t i = 0;
        const std::size_t n = text.size();
        while (i < n && !done_) {
            const char c = text[i];
            if (std::isspace(static_cast<unsigned char>(c))) {
                ++i;
                continue;
            }
            if (c == '%') {
                const std::size_t close = text.find('%', i + 1);
                if (close == std::string_view::npos) {
                    warn("unterminated extended command");
                    break;
                }
                extended(text.substr(i + 1, close - i - 1));
                i = close + 1;
            } else {
                const std::size_t star = text.find('*', i);
                if (star == std::string_view::npos) break;
                word(text.substr(i, star - i));
                i = star + 1;
            }
        }
        if (region_) warn("file ends inside a region statement");
        if (!blocks_.empty()) warn("file ends inside an aperture block");
        return std::move(img_);
    }

private:
    // ------------------------------------------------------------ helpers

    void warn(const std::string& msg) {
        if (img_.warnings.size() < 50) img_.warnings.push_back(msg);
    }

    std::vector<Feature>& out() { return blocks_.empty() ? img_.features : blocks_.back().features; }

    double length_nm(double value_in_units) const { return value_in_units * unit_nm_; }

    std::int64_t coordinate(const std::string& s) const {
        if (s.find('.') != std::string::npos) return std::llround(std::stod(s) * unit_nm_);
        std::string digits = s;
        bool negative = false;
        if (!digits.empty() && (digits[0] == '+' || digits[0] == '-')) {
            negative = digits[0] == '-';
            digits.erase(0, 1);
        }
        if (omit_trailing_) {
            const std::size_t total = static_cast<std::size_t>(int_digits_ + dec_digits_);
            if (digits.size() < total) digits.append(total - digits.size(), '0');
        }
        std::int64_t v = 0;
        std::from_chars(digits.data(), digits.data() + digits.size(), v);
        if (negative) v = -v;
        return std::llround(static_cast<double>(v) * unit_nm_ / std::pow(10.0, dec_digits_));
    }

    // ------------------------------------------------------------ attributes

    static std::pair<std::string, std::vector<std::string>> parse_attribute(std::string_view body) {
        auto fields = split(body, ',');
        std::string name = fields.empty() ? std::string() : fields[0];
        if (!name.empty() && name[0] == '.') name.erase(0, 1);
        std::vector<std::string> values;
        for (std::size_t i = 1; i < fields.size(); ++i) values.push_back(decode_attribute_field(fields[i]));
        return {name, values};
    }

    void attribute(std::string_view cmd) {
        const std::string code(cmd.substr(0, 2));
        const std::string_view body = cmd.substr(2);
        if (code == "TF") {
            auto [name, values] = parse_attribute(body);
            img_.file_attributes[name] = values;
        } else if (code == "TA") {
            auto [name, values] = parse_attribute(body);
            aperture_attrs_[name] = values;
        } else if (code == "TO") {
            auto [name, values] = parse_attribute(body);
            object_attrs_[name] = values;
            if (name == "N") img_.has_net_attributes = true;
        } else if (code == "TD") {
            std::string name = trim(body);
            if (!name.empty() && name[0] == '.') name.erase(0, 1);
            if (name.empty()) {
                aperture_attrs_.clear();
                object_attrs_.clear();
            } else {
                aperture_attrs_.erase(name);
                object_attrs_.erase(name);
            }
        }
    }

    void apply_object_attributes(Feature& f) const {
        if (auto it = object_attrs_.find("N"); it != object_attrs_.end()) {
            f.has_net = true;
            f.net = it->second.empty() ? std::string() : it->second[0];
        }
        if (auto it = object_attrs_.find("P"); it != object_attrs_.end()) {
            f.pad = true;
            if (!it->second.empty()) f.component = it->second[0];
            if (it->second.size() > 1) f.pin = it->second[1];
            if (it->second.size() > 2) f.pin_function = it->second[2];
        } else if (auto c = object_attrs_.find("C"); c != object_attrs_.end() && !c->second.empty()) {
            f.component = c->second[0];
        }
    }

    // ------------------------------------------------------------ extended commands

    void extended(std::string_view block) {
        const std::string body = strip_newlines(block);
        if (body.rfind("AM", 0) == 0) {
            define_macro(body);
            return;
        }
        for (const auto& w : split(body, '*')) {
            const std::string cmd = trim(w);
            if (cmd.size() < 2) continue;
            extended_command(cmd);
        }
    }

    void extended_command(const std::string& cmd) {
        const std::string code = cmd.substr(0, 2);
        if (code == "FS") {
            format(cmd.substr(2));
        } else if (code == "MO") {
            const std::string u = cmd.substr(2);
            if (u == "MM") unit_nm_ = kNmPerMm;
            else if (u == "IN") unit_nm_ = 25.4 * kNmPerMm;
            else warn("unknown unit " + u);
        } else if (code == "AD") {
            define_aperture(cmd.substr(2));
        } else if (code == "AB") {
            aperture_block(cmd.substr(2));
        } else if (code == "LP") {
            dark_ = cmd.substr(2) != "C";
        } else if (code == "LM") {
            const std::string m = cmd.substr(2);
            mirror_x_ = m.find('X') != std::string::npos;
            mirror_y_ = m.find('Y') != std::string::npos;
        } else if (code == "LR") {
            rotation_ = to_double(cmd.substr(2));
        } else if (code == "LS") {
            scale_ = to_double(cmd.substr(2));
        } else if (code == "SR") {
            step_repeat(cmd.substr(2));
        } else if (code == "TF" || code == "TA" || code == "TO" || code == "TD") {
            attribute(cmd);
        } else if (code == "IP") {
            if (cmd.substr(2) == "NEG") warn("negative image polarity (IPNEG) is not supported");
        } else if (code == "OF" || code == "SF" || code == "MI" || code == "IR") {
            if (cmd.find_first_of("123456789", 2) != std::string::npos)
                warn("deprecated image transformation " + code + " is ignored");
        } else if (code == "IN" || code == "LN" || code == "AS" || code == "IJ") {
            // Image name, level name, axis select: no effect on geometry here.
        } else {
            warn("unknown extended command " + code);
        }
    }

    void format(const std::string& f) {
        omit_trailing_ = f.find('T') != std::string::npos && f.find('L') == std::string::npos;
        incremental_ = f.size() > 1 && f[1] == 'I';
        const std::size_t x = f.find('X');
        if (x != std::string::npos && x + 2 < f.size()) {
            int_digits_ = f[x + 1] - '0';
            dec_digits_ = f[x + 2] - '0';
        }
    }

    void define_macro(const std::string& body) {
        auto words = split(body, '*');
        const std::string name = trim(words[0].substr(2));
        Macro macro;
        for (std::size_t i = 1; i < words.size(); ++i) {
            std::string w = trim(words[i]);
            if (w.empty()) continue;
            // Comment primitive: code 0 followed by text.
            if (w[0] == '0' && (w.size() == 1 || !(std::isdigit(static_cast<unsigned char>(w[1])) ||
                                                   w[1] == '.' || w[1] == ','))) {
                continue;
            }
            std::string compact;
            for (char c : w)
                if (!std::isspace(static_cast<unsigned char>(c))) compact.push_back(c);
            MacroStatement st;
            if (compact[0] == '$') {
                const std::size_t eq = compact.find('=');
                if (eq == std::string::npos) {
                    warn("macro " + name + ": malformed statement " + compact);
                    continue;
                }
                st.assignment = true;
                st.variable = std::stoi(compact.substr(1, eq - 1));
                st.exprs.push_back(compact.substr(eq + 1));
            } else {
                st.exprs = split(compact, ',');
            }
            macro.statements.push_back(std::move(st));
        }
        macros_[name] = std::move(macro);
    }

    Paths64 expand_macro(const std::string& name, const Macro& macro, std::vector<double> vars) {
        Paths64 acc;
        Paths64 pending;
        auto flush = [&]() {
            if (pending.empty()) return;
            acc.insert(acc.end(), pending.begin(), pending.end());
            acc = to_paths(union_polygons(acc));
            pending.clear();
        };
        for (const auto& st : macro.statements) {
            try {
                if (st.assignment) {
                    if (static_cast<int>(vars.size()) <= st.variable) vars.resize(st.variable + 1, 0.0);
                    vars[st.variable] = evaluate_macro_expression(st.exprs[0], vars);
                    continue;
                }
                std::vector<double> v;
                for (std::size_t k = 1; k < st.exprs.size(); ++k)
                    v.push_back(st.exprs[k].empty() ? 0.0 : evaluate_macro_expression(st.exprs[k], vars));
                const int code = static_cast<int>(evaluate_macro_expression(st.exprs[0], vars));
                bool exposure = true;
                Paths64 shape = primitive(code, v, exposure);
                if (shape.empty()) continue;
                for (auto& p : shape) make_positive(p);
                if (exposure) {
                    pending.insert(pending.end(), shape.begin(), shape.end());
                } else {
                    flush();
                    acc = to_paths(difference(union_polygons(acc), union_polygons(shape)));
                }
            } catch (const std::exception& e) {
                warn("macro " + name + ": " + e.what());
            }
        }
        flush();
        return acc;
    }

    // Macro primitive at the macro origin, in nm. Values are in file units.
    Paths64 primitive(int code, const std::vector<double>& v, bool& exposure) {
        auto at = [&](std::size_t i) { return i < v.size() ? v[i] : 0.0; };
        auto L = [&](std::size_t i) { return length_nm(at(i)); };
        switch (code) {
        case 1: {  // circle: exposure, diameter, x, y, [rotation]
            exposure = at(0) != 0.0;
            Paths64 s{circle_path(Point64(std::llround(L(2)), std::llround(L(3))), L(1) / 2, tol_)};
            return rotate_about_origin(s, at(4));
        }
        case 2:
        case 20: {  // vector line: exposure, width, x0, y0, x1, y1, rotation
            exposure = at(0) != 0.0;
            const double w = L(1), x0 = L(2), y0 = L(3), x1 = L(4), y1 = L(5);
            const double len = std::hypot(x1 - x0, y1 - y0);
            if (len == 0.0 || w == 0.0) return {};
            const double nx = -(y1 - y0) / len * w / 2, ny = (x1 - x0) / len * w / 2;
            Path64 p{Point64(std::llround(x0 + nx), std::llround(y0 + ny)),
                     Point64(std::llround(x0 - nx), std::llround(y0 - ny)),
                     Point64(std::llround(x1 - nx), std::llround(y1 - ny)),
                     Point64(std::llround(x1 + nx), std::llround(y1 + ny))};
            return rotate_about_origin({p}, at(6));
        }
        case 21: {  // centre line: exposure, width, height, x, y, rotation
            exposure = at(0) != 0.0;
            return rotate_about_origin({rect_path(L(3), L(4), L(1), L(2))}, at(5));
        }
        case 22: {  // lower-left line (deprecated): exposure, width, height, x, y, rotation
            exposure = at(0) != 0.0;
            return rotate_about_origin({rect_path(L(3) + L(1) / 2, L(4) + L(2) / 2, L(1), L(2))}, at(5));
        }
        case 4: {  // outline: exposure, n, x0, y0, ..., xn, yn, rotation
            exposure = at(0) != 0.0;
            const int n = static_cast<int>(at(1));
            Path64 p;
            for (int k = 0; k <= n; ++k)
                p.emplace_back(std::llround(L(2 + 2 * k)), std::llround(L(3 + 2 * k)));
            if (p.size() > 1 && p.front() == p.back()) p.pop_back();
            return rotate_about_origin({p}, at(2 + 2 * (n + 1)));
        }
        case 5: {  // polygon: exposure, vertices, x, y, diameter, rotation
            exposure = at(0) != 0.0;
            const Path64 p = regular_polygon(L(2), L(3), L(4), static_cast<int>(at(1)), 0.0);
            return rotate_about_origin({p}, at(5));
        }
        case 6: {  // moire (deprecated): x, y, outer diameter, ring thickness, gap, rings,
                   // crosshair thickness, crosshair length, rotation
            exposure = true;
            const Point64 c(std::llround(L(0)), std::llround(L(1)));
            Paths64 rings;
            double od = L(2);
            for (int k = 0; k < static_cast<int>(at(5)) && od > 0; ++k) {
                const double id = std::max(0.0, od - 2 * L(3));
                Paths64 ring{circle_path(c, od / 2, tol_)};
                if (id > 0) {
                    Path64 hole = circle_path(c, id / 2, tol_);
                    std::reverse(hole.begin(), hole.end());
                    ring.push_back(hole);
                }
                rings.insert(rings.end(), ring.begin(), ring.end());
                od = id - 2 * L(4);
            }
            rings.push_back(rect_path(L(0), L(1), L(7), L(6)));
            rings.push_back(rect_path(L(0), L(1), L(6), L(7)));
            return rotate_about_origin(to_paths(union_polygons(rings)), at(8));
        }
        case 7: {  // thermal: x, y, outer diameter, inner diameter, gap, rotation
            exposure = true;
            const Point64 c(std::llround(L(0)), std::llround(L(1)));
            Paths64 ring{circle_path(c, L(2) / 2, tol_)};
            Path64 hole = circle_path(c, L(3) / 2, tol_);
            std::reverse(hole.begin(), hole.end());
            ring.push_back(hole);
            const Paths64 cross{rect_path(L(0), L(1), L(2) * 1.1, L(4)), rect_path(L(0), L(1), L(4), L(2) * 1.1)};
            const PolygonSet result = difference(union_polygons(ring), union_polygons(cross));
            return rotate_about_origin(to_paths(result), at(5));
        }
        default:
            warn("unsupported macro primitive " + std::to_string(code));
            return {};
        }
    }

    void define_aperture(const std::string& def) {
        // D<code><name>[,<p1>X<p2>...]
        std::size_t i = 1;
        while (i < def.size() && std::isdigit(static_cast<unsigned char>(def[i]))) ++i;
        const int dcode = std::stoi(def.substr(1, i - 1));
        const std::size_t comma = def.find(',', i);
        const std::string name = def.substr(i, comma == std::string::npos ? std::string::npos : comma - i);
        std::vector<double> params;
        if (comma != std::string::npos)
            for (const auto& p : split(def.substr(comma + 1), 'X')) params.push_back(to_double(p));
        auto P = [&](std::size_t k) { return k < params.size() ? params[k] : 0.0; };

        Aperture ap;
        if (auto it = aperture_attrs_.find("AperFunction"); it != aperture_attrs_.end() && !it->second.empty())
            ap.function = it->second[0];
        std::size_t hole_index = 0;
        if (name == "C") {
            ap.type = ApertureType::Circle;
            ap.circle_diameter_nm = length_nm(P(0));
            if (ap.circle_diameter_nm > 0) ap.shape = {circle_path(Point64(0, 0), ap.circle_diameter_nm / 2, tol_)};
            hole_index = 1;
        } else if (name == "R") {
            ap.type = ApertureType::Rect;
            ap.shape = {rect_path(0, 0, length_nm(P(0)), length_nm(P(1)))};
            hole_index = 2;
        } else if (name == "O") {
            ap.type = ApertureType::Obround;
            const double w = length_nm(P(0)), h = length_nm(P(1));
            if (w > h) {
                const auto d = std::llround((w - h) / 2);
                ap.shape = {capsule(Point64(-d, 0LL), Point64(d, 0LL), h / 2, tol_)};
            } else if (h > w) {
                const auto d = std::llround((h - w) / 2);
                ap.shape = {capsule(Point64(0LL, -d), Point64(0LL, d), w / 2, tol_)};
            } else {
                ap.shape = {circle_path(Point64(0, 0), w / 2, tol_)};
            }
            hole_index = 2;
        } else if (name == "P") {
            ap.type = ApertureType::Polygon;
            ap.shape = {regular_polygon(0, 0, length_nm(P(0)), static_cast<int>(P(1)), P(2))};
            hole_index = 3;
        } else if (auto m = macros_.find(name); m != macros_.end()) {
            ap.type = ApertureType::Macro;
            std::vector<double> vars(params.size() + 1, 0.0);
            for (std::size_t k = 0; k < params.size(); ++k) vars[k + 1] = params[k];
            ap.shape = expand_macro(name, m->second, vars);
            hole_index = SIZE_MAX;
        } else {
            warn("aperture D" + std::to_string(dcode) + " uses unknown template " + name);
            return;
        }
        if (hole_index < params.size() && params[hole_index] > 0) {
            Paths64 hole{circle_path(Point64(0, 0), length_nm(params[hole_index]) / 2, tol_)};
            ap.shape = to_paths(difference(union_polygons(ap.shape), union_polygons(hole)));
        }
        apertures_[dcode] = std::move(ap);
    }

    void aperture_block(const std::string& arg) {
        if (arg.empty()) {
            if (blocks_.empty()) {
                warn("aperture block end without start");
                return;
            }
            Block b = std::move(blocks_.back());
            blocks_.pop_back();
            Aperture ap;
            ap.type = ApertureType::Block;
            ap.block = std::move(b.features);
            apertures_[b.dcode] = std::move(ap);
        } else {
            blocks_.push_back({std::stoi(arg.substr(1)), {}});
        }
    }

    void step_repeat(const std::string& arg) {
        close_step_repeat();
        if (arg.empty()) return;
        int nx = 1, ny = 1;
        double dx = 0, dy = 0;
        std::size_t i = 0;
        while (i < arg.size()) {
            const char k = arg[i++];
            std::size_t j = i;
            while (j < arg.size() && (std::isdigit(static_cast<unsigned char>(arg[j])) || arg[j] == '.' ||
                                      arg[j] == '-' || arg[j] == '+'))
                ++j;
            const std::string num = arg.substr(i, j - i);
            i = j;
            if (k == 'X') nx = std::stoi(num);
            else if (k == 'Y') ny = std::stoi(num);
            else if (k == 'I') dx = std::stod(num);
            else if (k == 'J') dy = std::stod(num);
        }
        if (nx > 1 || ny > 1) sr_ = StepRepeat{nx, ny, length_nm(dx), length_nm(dy), out().size()};
    }

    void close_step_repeat() {
        if (!sr_) return;
        auto& features = out();
        const std::vector<Feature> block(features.begin() + static_cast<std::ptrdiff_t>(sr_->start), features.end());
        for (int j = 0; j < sr_->ny; ++j) {
            for (int i = 0; i < sr_->nx; ++i) {
                if (i == 0 && j == 0) continue;
                Placement t;
                t.offset = Point64(std::llround(i * sr_->dx), std::llround(j * sr_->dy));
                for (const auto& f : block) features.push_back(placed(f, t, true));
            }
        }
        sr_.reset();
    }

    // ------------------------------------------------------------ operations

    void word(std::string_view raw) {
        const std::string w = strip_newlines(raw);
        if (w.empty()) return;
        if (w.rfind("G04", 0) == 0) {
            // KiCad writes X2 attributes as comments "G04 #@! TF..." in compatibility mode.
            const std::size_t mark = w.find("#@!");
            if (mark != std::string::npos) attribute(trim(std::string_view(w).substr(mark + 3)));
            return;
        }
        if (w == "M02" || w == "M00" || w == "M01") {
            close_step_repeat();
            done_ = true;
            return;
        }
        std::optional<std::string> x, y, xi, yj;
        int d = -1;
        std::size_t i = 0;
        while (i < w.size()) {
            const char letter = w[i++];
            std::size_t j = i;
            while (j < w.size() && (std::isdigit(static_cast<unsigned char>(w[j])) || w[j] == '-' ||
                                    w[j] == '+' || w[j] == '.'))
                ++j;
            const std::string num = w.substr(i, j - i);
            i = j;
            switch (letter) {
            case 'G': g_code(num.empty() ? -1 : std::stoi(num)); break;
            case 'X': x = num; break;
            case 'Y': y = num; break;
            case 'I': xi = num; break;
            case 'J': yj = num; break;
            case 'D': d = num.empty() ? -1 : std::stoi(num); break;
            case 'M': break;
            default: warn(std::string("unknown word ") + w); return;
            }
        }
        if (d >= 10) {
            current_ = d;
            if (!apertures_.count(d)) warn("aperture D" + std::to_string(d) + " is not defined");
            return;
        }
        const bool has_coord = x || y || xi || yj;
        if (d < 0 && !has_coord) return;
        if (d < 0) d = last_d_ > 0 ? last_d_ : 1;  // deprecated modal operation code
        last_d_ = d;

        Point64 target = pos_;
        if (incremental_) {
            if (x) target.x += coordinate(*x);
            if (y) target.y += coordinate(*y);
        } else {
            if (x) target.x = coordinate(*x);
            if (y) target.y = coordinate(*y);
        }
        const Point64 offset(xi ? coordinate(*xi) : 0, yj ? coordinate(*yj) : 0);

        switch (d) {
        case 1: interpolate(target, offset); break;
        case 2:
            if (region_) close_contour();
            if (region_) contour_ = {target};
            break;
        case 3:
            if (region_) warn("flash inside a region statement");
            else flash(target);
            break;
        default: warn("unknown operation D" + std::to_string(d));
        }
        pos_ = target;
    }

    void g_code(int g) {
        switch (g) {
        case 1: interp_ = Interp::Linear; break;
        case 2: interp_ = Interp::Clockwise; break;
        case 3: interp_ = Interp::CounterClockwise; break;
        case 74: multi_quadrant_ = false; break;
        case 75: multi_quadrant_ = true; break;
        case 36:
            region_ = true;
            contour_.clear();
            break;
        case 37:
            close_contour();
            region_ = false;
            contour_.clear();
            break;
        case 70: unit_nm_ = 25.4 * kNmPerMm; break;
        case 71: unit_nm_ = kNmPerMm; break;
        case 90: incremental_ = false; break;
        case 91: incremental_ = true; break;
        case 4:
        case 54:
        case 55:
        case -1: break;
        default: warn("unsupported G code G" + std::to_string(g));
        }
    }

    Point64 arc_centre(Point64 start, Point64 end, Point64 offset, bool ccw) const {
        if (multi_quadrant_) return Point64(start.x + offset.x, start.y + offset.y);
        // Single quadrant: offsets are unsigned; take the centre giving an arc of at most
        // 90 degrees in the given direction with the most equal radii.
        Point64 best = start;
        double best_err = INFINITY;
        for (int sx : {-1, 1}) {
            for (int sy : {-1, 1}) {
                const Point64 c(start.x + sx * std::llabs(offset.x), start.y + sy * std::llabs(offset.y));
                const double a0 = std::atan2(double(start.y - c.y), double(start.x - c.x));
                const double a1 = std::atan2(double(end.y - c.y), double(end.x - c.x));
                double sweep = ccw ? a1 - a0 : a0 - a1;
                while (sweep < 0) sweep += 2 * kPi;
                if (sweep > kPi / 2 + 1e-6) continue;
                const double err = std::abs(std::hypot(double(start.x - c.x), double(start.y - c.y)) -
                                            std::hypot(double(end.x - c.x), double(end.y - c.y)));
                if (err < best_err) {
                    best_err = err;
                    best = c;
                }
            }
        }
        return best;
    }

    void interpolate(Point64 target, Point64 offset) {
        const bool arc = interp_ != Interp::Linear;
        const bool ccw = interp_ == Interp::CounterClockwise;
        const Point64 centre = arc ? arc_centre(pos_, target, offset, ccw) : Point64(0, 0);

        if (region_) {
            if (contour_.empty()) contour_.push_back(pos_);
            if (arc) {
                Path64 pts = arc_points(pos_, target, centre, ccw, tol_, true);
                contour_.insert(contour_.end(), pts.begin() + 1, pts.end());
            } else {
                contour_.push_back(target);
            }
            return;
        }

        auto it = apertures_.find(current_);
        if (it == apertures_.end()) {
            warn("stroke without a defined aperture");
            return;
        }
        const Aperture& ap = it->second;
        Feature f;
        f.kind = FeatureKind::Stroke;
        f.dark = dark_;
        f.at = pos_;
        f.end = target;
        if (arc) {
            f.arc_centre = centre;
            f.arc_ccw = ccw;
        }
        f.function = ap.function;
        apply_object_attributes(f);

        if (ap.type == ApertureType::Circle) {
            const double r = ap.circle_diameter_nm * scale_ / 2;
            f.stroke_width_nm = 2 * r;
            if (r > 0) {
                if (!arc) {
                    f.shape = {capsule(pos_, target, r, tol_)};
                } else {
                    Path64 pts = arc_points(pos_, target, centre, ccw, tol_);
                    if (pos_ == target) {
                        pts.pop_back();
                        f.shape = stroke_ring(pts, r, tol_);
                    } else {
                        f.shape = stroke_polyline(pts, r, tol_);
                    }
                }
            }
        } else if (ap.type != ApertureType::Block && !ap.shape.empty()) {
            // Non-circular apertures: the swept convex hull of the aperture.
            Path64 hull;
            for (const auto& p : flash_placement_shape(ap.shape, Point64(0, 0)))
                hull.insert(hull.end(), p.begin(), p.end());
            hull = convex_hull(hull);
            const Path64 pts = arc ? arc_points(pos_, target, centre, ccw, tol_) : Path64{pos_, target};
            Paths64 parts;
            for (std::size_t k = 0; k + 1 < pts.size(); ++k) parts.push_back(stroke_convex(pts[k], pts[k + 1], hull));
            f.shape = parts.size() == 1 ? parts : to_paths(union_polygons(parts));
            if (ap.type != ApertureType::Rect || arc) warn("stroke with a non-circular aperture approximated by its convex hull");
        }
        out().push_back(std::move(f));
    }

    Paths64 flash_placement_shape(const Paths64& shape, Point64 at) const {
        Placement t;
        t.scale = scale_;
        t.mirror_x = mirror_x_;
        t.mirror_y = mirror_y_;
        t.rotation_deg = rotation_;
        t.offset = at;
        return transform(shape, t);
    }

    Feature placed(const Feature& f, const Placement& t, bool keep_polarity) const {
        Feature g = f;
        g.shape = transform(f.shape, t);
        g.at = transform(f.at, t);
        g.end = transform(f.end, t);
        if (f.arc_centre) {
            g.arc_centre = transform(*f.arc_centre, t);
            if (t.mirror_x != t.mirror_y) g.arc_ccw = !f.arc_ccw;
        }
        g.stroke_width_nm = f.stroke_width_nm * t.scale;
        if (!keep_polarity) g.dark = (f.dark == dark_);
        return g;
    }

    void flash(Point64 at) {
        auto it = apertures_.find(current_);
        if (it == apertures_.end()) {
            warn("flash without a defined aperture");
            return;
        }
        const Aperture& ap = it->second;
        if (ap.type == ApertureType::Block) {
            Placement t;
            t.scale = scale_;
            t.mirror_x = mirror_x_;
            t.mirror_y = mirror_y_;
            t.rotation_deg = rotation_;
            t.offset = at;
            for (const auto& f : ap.block) {
                Feature g = placed(f, t, false);
                if (!g.has_net) apply_object_attributes(g);
                out().push_back(std::move(g));
            }
            return;
        }
        Feature f;
        f.kind = FeatureKind::Flash;
        f.dark = dark_;
        f.at = at;
        f.function = ap.function;
        f.shape = flash_placement_shape(ap.shape, at);
        apply_object_attributes(f);
        out().push_back(std::move(f));
    }

    void close_contour() {
        if (contour_.size() >= 3) {
            Path64 p = contour_;
            if (p.front() == p.back()) p.pop_back();
            if (p.size() >= 3) {
                make_positive(p);
                Feature f;
                f.kind = FeatureKind::Region;
                f.dark = dark_;
                f.at = p.front();
                f.shape = {std::move(p)};
                if (auto a = aperture_attrs_.find("AperFunction"); a != aperture_attrs_.end() && !a->second.empty())
                    f.function = a->second[0];
                apply_object_attributes(f);
                out().push_back(std::move(f));
            }
        }
        contour_.clear();
    }

    // ------------------------------------------------------------ state

    enum class Interp { Linear, Clockwise, CounterClockwise };
    struct Block {
        int dcode;
        std::vector<Feature> features;
    };
    struct StepRepeat {
        int nx, ny;
        double dx, dy;
        std::size_t start;
    };

    double tol_;
    double unit_nm_ = kNmPerMm;
    int int_digits_ = 4;
    int dec_digits_ = 6;
    bool omit_trailing_ = false;
    bool incremental_ = false;
    Point64 pos_{0, 0};
    Interp interp_ = Interp::Linear;
    bool multi_quadrant_ = true;
    int current_ = -1;
    int last_d_ = 0;
    bool region_ = false;
    Path64 contour_;
    bool dark_ = true;
    bool mirror_x_ = false;
    bool mirror_y_ = false;
    double rotation_ = 0.0;
    double scale_ = 1.0;
    bool done_ = false;
    std::map<std::string, std::vector<std::string>> aperture_attrs_;
    std::map<std::string, std::vector<std::string>> object_attrs_;
    std::map<int, Aperture> apertures_;
    std::map<std::string, Macro> macros_;
    std::vector<Block> blocks_;
    std::optional<StepRepeat> sr_;
    GerberImage img_;
};

void append_utf8(std::string& out, unsigned cp) {
    if (cp < 0x80) {
        out.push_back(static_cast<char>(cp));
    } else if (cp < 0x800) {
        out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else {
        out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    }
}

}  // namespace

std::string decode_attribute_field(std::string_view field) {
    std::string out;
    out.reserve(field.size());
    for (std::size_t i = 0; i < field.size(); ++i) {
        if (field[i] == '\\' && i + 5 < field.size() + 0 && (field[i + 1] == 'u' || field[i + 1] == 'U')) {
            unsigned cp = 0;
            const auto r = std::from_chars(field.data() + i + 2, field.data() + i + 6, cp, 16);
            if (r.ec == std::errc() && r.ptr == field.data() + i + 6) {
                append_utf8(out, cp);
                i += 5;
                continue;
            }
        }
        out.push_back(field[i]);
    }
    return out;
}

double evaluate_macro_expression(std::string_view expr, const std::vector<double>& vars) {
    return ExpressionParser(expr, vars).parse();
}

std::vector<std::string> GerberImage::file_function() const {
    auto it = file_attributes.find("FileFunction");
    return it == file_attributes.end() ? std::vector<std::string>{} : it->second;
}

std::string GerberImage::file_polarity() const {
    auto it = file_attributes.find("FilePolarity");
    return it == file_attributes.end() || it->second.empty() ? std::string() : it->second[0];
}

GerberImage parse_gerber(std::string_view text, const GeometryOptions& options) {
    return Parser(options).run(text);
}

bool looks_like_gerber(std::string_view text) {
    const std::string_view head = text.substr(0, 20000);
    return head.find("%FS") != std::string_view::npos || head.find("%MO") != std::string_view::npos;
}

}  // namespace kairoq
