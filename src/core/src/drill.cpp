#include "kairoq/drill.hpp"

#include <cctype>
#include <cmath>
#include <map>
#include <sstream>

namespace kairoq {

namespace {

std::string trim(std::string_view s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return std::string(s.substr(a, b - a));
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, sep)) out.push_back(trim(item));
    return out;
}

bool is_number_char(char c) { return std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+'; }

class ExcellonParser {
public:
    DrillFile run(std::string_view text) {
        std::istringstream in{std::string(text)};
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty()) continue;
            if (line[0] == ';') {
                comment(line.substr(1));
                continue;
            }
            if (header_) header_line(line);
            else body_line(line);
        }
        return std::move(file_);
    }

private:
    void comment(const std::string& c) {
        const std::size_t mark = c.find("#@!");
        if (mark != std::string::npos) {
            const std::string attr = trim(c.substr(mark + 3));
            if (attr.rfind("TF.FileFunction,", 0) == 0) {
                apply_drill_file_function(file_, split(attr.substr(16), ','));
            } else if (attr.rfind("TA.AperFunction,", 0) == 0) {
                const auto fields = split(attr.substr(16), ',');
                pending_function_ = fields.empty() ? std::string() : fields.back();
            }
            return;
        }
        // KiCad: FORMAT={3:3/ absolute / metric / suppress leading zeros}
        const std::size_t f = c.find("FORMAT={");
        if (f != std::string::npos) {
            const std::string spec = c.substr(f + 8);
            if (spec.size() >= 3 && std::isdigit(static_cast<unsigned char>(spec[0])) && spec[1] == ':' &&
                std::isdigit(static_cast<unsigned char>(spec[2]))) {
                int_digits_ = spec[0] - '0';
                dec_digits_ = spec[2] - '0';
                digits_given_ = true;
            }
            if (spec.find("suppress trailing") != std::string::npos) zeros_ = Zeros::Leading;
            else if (spec.find("suppress leading") != std::string::npos) zeros_ = Zeros::Trailing;
            return;
        }
        // Altium and others: FILE_FORMAT=2:5
        const std::size_t ff = c.find("FILE_FORMAT=");
        if (ff != std::string::npos && ff + 15 <= c.size()) {
            int_digits_ = c[ff + 12] - '0';
            dec_digits_ = c[ff + 14] - '0';
            digits_given_ = true;
        }
    }

    void units(const std::string& line) {
        const auto fields = split(line, ',');
        metric_ = fields[0] == "METRIC";
        if (!digits_given_) {
            int_digits_ = metric_ ? 3 : 2;
            dec_digits_ = metric_ ? 3 : 4;
        }
        for (std::size_t i = 1; i < fields.size(); ++i) {
            if (fields[i] == "LZ") zeros_ = Zeros::Leading;
            else if (fields[i] == "TZ") zeros_ = Zeros::Trailing;
            else if (fields[i].find('.') != std::string::npos) {
                const std::size_t dot = fields[i].find('.');
                int_digits_ = static_cast<int>(dot);
                dec_digits_ = static_cast<int>(fields[i].size() - dot - 1);
                digits_given_ = true;
            }
        }
    }

    void header_line(const std::string& line) {
        if (line == "M48") return;
        if (line == "%" || line == "M95") {
            header_ = false;
            return;
        }
        if (line.rfind("METRIC", 0) == 0 || line.rfind("INCH", 0) == 0) {
            units(line);
        } else if (line == "M71") {
            units("METRIC");
        } else if (line == "M72") {
            units("INCH");
        } else if (line[0] == 'T') {
            tool_definition(line);
        }
        // FMAT, VER, ICI and other header words do not affect geometry.
    }

    double value_nm(const std::string& s) const {
        const double unit = metric_ ? kNmPerMm : 25.4 * kNmPerMm;
        if (s.find('.') != std::string::npos) return std::stod(s) * unit;
        std::string digits = s;
        bool negative = false;
        if (!digits.empty() && (digits[0] == '-' || digits[0] == '+')) {
            negative = digits[0] == '-';
            digits.erase(0, 1);
        }
        if (digits.empty()) return 0.0;
        double v = std::stod(digits);
        if (zeros_ == Zeros::Leading) {
            // Leading zeros present, trailing suppressed: digits are left aligned.
            v /= std::pow(10.0, static_cast<double>(digits.size()) - int_digits_);
        } else {
            v /= std::pow(10.0, dec_digits_);
        }
        return (negative ? -v : v) * unit;
    }

    // Letter/value pairs of a word such as T1C0.800 or X1.0Y2.0G85X3.0Y2.0.
    static std::vector<std::pair<char, std::string>> words(const std::string& line) {
        std::vector<std::pair<char, std::string>> out;
        std::size_t i = 0;
        while (i < line.size()) {
            const char letter = line[i++];
            std::size_t j = i;
            while (j < line.size() && is_number_char(line[j])) ++j;
            out.emplace_back(letter, line.substr(i, j - i));
            i = j;
        }
        return out;
    }

    void tool_definition(const std::string& line) {
        int number = -1;
        double diameter = -1;
        for (const auto& [letter, v] : words(line)) {
            if (letter == 'T' && !v.empty()) number = std::stoi(v);
            else if (letter == 'C' && !v.empty()) diameter = std::stod(v) * (metric_ ? kNmPerMm : 25.4 * kNmPerMm);
        }
        if (number < 0) return;
        if (diameter < 0) {
            select_tool(number);
            return;
        }
        DrillTool t;
        t.number = number;
        t.diameter_nm = diameter;
        t.function = pending_function_;
        pending_function_.clear();
        tool_index_[number] = static_cast<int>(file_.tools.size());
        file_.tools.push_back(std::move(t));
        if (!header_) select_tool(number);
    }

    void select_tool(int number) {
        if (number == 0) {
            current_ = -1;
            return;
        }
        auto it = tool_index_.find(number);
        if (it == tool_index_.end()) {
            file_.warnings.push_back("tool T" + std::to_string(number) + " is not defined");
            current_ = -1;
        } else {
            current_ = it->second;
        }
    }

    void body_line(const std::string& line) {
        if (line == "M30" || line == "M00") return;
        if (line[0] == 'T') {
            tool_definition(line);
            return;
        }
        if (line == "M15") {
            tool_down_ = true;
            return;
        }
        if (line == "M16" || line == "M17") {
            tool_down_ = false;
            return;
        }
        auto ws = words(line);
        Point64 p = pos_;
        bool coords = false;
        bool g85 = false;
        Point64 first = pos_;
        int g = -1;
        for (const auto& [letter, v] : ws) {
            if (letter == 'G') {
                const int code = v.empty() ? -1 : std::stoi(v);
                if (code == 85) {
                    g85 = true;
                    first = p;
                } else {
                    g = code;
                }
            } else if (letter == 'X') {
                p.x = std::llround(value_nm(v));
                coords = true;
            } else if (letter == 'Y') {
                p.y = std::llround(value_nm(v));
                coords = true;
            }
        }
        if (g == 0) route_ = true;
        if (g == 5) route_ = false;
        if (g == 90 || g == 91) {
            if (g == 91) file_.warnings.push_back("incremental Excellon coordinates are not supported");
        }
        if (g == 2 || g == 3) file_.warnings.push_back("routed arcs are approximated by straight slots");
        if (!coords) return;
        if (current_ < 0) {
            file_.warnings.push_back("hit without a selected tool");
            pos_ = p;
            return;
        }
        if (g85) {
            file_.hits.push_back({current_, first, p});
        } else if (route_) {
            if (tool_down_ && (g == 1 || g == 2 || g == 3 || g == -1)) file_.hits.push_back({current_, pos_, p});
        } else {
            file_.hits.push_back({current_, p, std::nullopt});
        }
        pos_ = p;
    }

    enum class Zeros { Leading, Trailing };

    DrillFile file_;
    bool header_ = true;
    bool metric_ = true;
    int int_digits_ = 3;
    int dec_digits_ = 3;
    bool digits_given_ = false;
    Zeros zeros_ = Zeros::Trailing;
    std::string pending_function_;
    std::map<int, int> tool_index_;
    int current_ = -1;
    Point64 pos_{0, 0};
    bool route_ = false;
    bool tool_down_ = false;
};

}  // namespace

void apply_drill_file_function(DrillFile& drill, const std::vector<std::string>& fields) {
    if (fields.empty()) return;
    if (fields[0] == "Plated") drill.plated = true;
    else if (fields[0] == "NonPlated") drill.plated = false;
    else return;
    if (fields.size() > 2) {
        try {
            drill.from_layer = std::stoi(fields[1]);
            drill.to_layer = std::stoi(fields[2]);
        } catch (const std::exception&) {
            drill.warnings.push_back("malformed drill span in .FileFunction");
        }
    }
    if (fields.size() > 3) drill.span_type = fields[3];
}

DrillFile parse_excellon(std::string_view text) { return ExcellonParser().run(text); }

bool looks_like_excellon(std::string_view text) {
    const std::string_view head = text.substr(0, 4096);
    return head.find("M48") != std::string_view::npos;
}

DrillFile drill_from_gerber(const GerberImage& image) {
    DrillFile d;
    apply_drill_file_function(d, image.file_function());
    std::map<std::pair<long long, std::string>, int> tools;
    for (const auto& f : image.features) {
        double dia = f.stroke_width_nm;
        if (f.kind == FeatureKind::Flash) {
            // Diameter of the flashed circle from its area.
            dia = 2.0 * std::sqrt(std::abs(area_nm2(union_polygons(f.shape))) / kPi);
        }
        // Round to 0.1 um to group equal tools.
        const long long key = std::llround(dia / 100.0);
        auto [it, inserted] = tools.emplace(std::make_pair(key, f.function), static_cast<int>(d.tools.size()));
        if (inserted) d.tools.push_back({static_cast<int>(d.tools.size()) + 1, key * 100.0, f.function});
        if (f.kind == FeatureKind::Flash) d.hits.push_back({it->second, f.at, std::nullopt});
        else if (f.kind == FeatureKind::Stroke) d.hits.push_back({it->second, f.at, f.end});
    }
    d.warnings = image.warnings;
    return d;
}

}  // namespace kairoq
