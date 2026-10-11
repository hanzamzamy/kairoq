// Gerber job file, IPC-D-356 netlist, placement CSV and stackup YAML readers.
#include "kairoq/fabrication.hpp"

#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

std::string trim(std::string_view s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return std::string(s.substr(a, b - a));
}

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == sep) {
            out.push_back(trim(std::string_view(s).substr(start, i - start)));
            start = i + 1;
        }
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------- job file

namespace {

// Numeric job file field; some generators write numbers as strings.
std::optional<double> json_number(const nlohmann::json& j, const char* key) {
    if (!j.contains(key)) return std::nullopt;
    const auto& v = j[key];
    if (v.is_number()) return v.get<double>();
    if (v.is_string()) {
        try {
            return std::stod(v.get<std::string>());
        } catch (const std::exception&) {
        }
    }
    return std::nullopt;
}

}  // namespace

JobFile parse_job_file(std::string_view text) {
    const auto j = nlohmann::json::parse(text);
    JobFile job;
    if (j.contains("Header") && j["Header"].contains("GenerationSoftware")) {
        const auto& g = j["Header"]["GenerationSoftware"];
        job.generator = g.value("Vendor", "") + " " + g.value("Application", "") + " " + g.value("Version", "");
        job.generator = trim(job.generator);
    }
    if (j.contains("GeneralSpecs")) {
        const auto& g = j["GeneralSpecs"];
        if (g.contains("ProjectId")) job.project = g["ProjectId"].value("Name", "");
        job.layer_count = static_cast<int>(json_number(g, "LayerNumber").value_or(0));
        job.board_thickness_nm = json_number(g, "BoardThickness").value_or(0.0) * kNmPerMm;
    }
    if (j.contains("FilesAttributes")) {
        for (const auto& f : j["FilesAttributes"]) {
            JobFileEntry e;
            e.path = f.value("Path", "");
            e.function = split(f.value("FileFunction", ""), ',');
            e.polarity = f.value("FilePolarity", "");
            job.files.push_back(std::move(e));
        }
    }
    if (j.contains("MaterialStackup")) {
        int copper = 0;
        for (const auto& m : j["MaterialStackup"]) {
            const std::string type = m.value("Type", "");
            StackupLayer l;
            l.name = m.value("Name", "");
            l.thickness_nm = json_number(m, "Thickness").value_or(0.0) * kNmPerMm;
            l.material = m.value("Material", "");
            if (type == "Copper") {
                l.type = LayerType::Copper;
                l.copper_index = ++copper;
                if (l.material.empty()) l.material = "copper";
            } else if (type == "Dielectric") {
                l.type = LayerType::Dielectric;
            } else if (type == "SolderMask") {
                l.type = LayerType::SolderMask;
            } else {
                continue;  // legend and paste have no thickness in the stackup
            }
            l.epsilon_r = json_number(m, "DielectricConstant");
            l.loss_tangent = json_number(m, "LossTangent");
            job.stackup.push_back(std::move(l));
        }
    }
    return job;
}

// ---------------------------------------------------------------- IPC-D-356

bool looks_like_ipc356(std::string_view text) {
    const std::string_view head = text.substr(0, 4096);
    return head.find("UNITS CUST") != std::string_view::npos ||
           head.find("\n317") != std::string_view::npos || head.find("\n327") != std::string_view::npos ||
           head.rfind("P  JOB", 0) == 0 || head.rfind("C  ", 0) == 0;
}

Netlist parse_ipc356(std::string_view text) {
    Netlist nl;
    double unit_nm = 2540.0;  // CUST 0: 0.0001 in
    std::map<std::string, std::string> aliases;
    std::istringstream in{std::string(text)};
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.size() < 3) continue;
        if (line[0] == 'C') continue;
        if (line[0] == 'P') {
            const std::string rest = trim(line.substr(1));
            if (rest.rfind("UNITS", 0) == 0) {
                if (rest.find("CUST 1") != std::string::npos || rest.find("SI") != std::string::npos)
                    unit_nm = 1000.0;  // 0.001 mm
                else
                    unit_nm = 2540.0;
            } else if (rest.rfind("NNAME", 0) == 0) {
                // "NNAME<alias> <name>" (KiCad) or "NNAME <alias> <name>": a net name too long
                // for the 14-character field, used by its alias in the records.
                std::istringstream ss(rest.substr(5));
                std::string alias, name;
                ss >> alias;
                std::getline(ss, name);
                if (!alias.empty()) aliases[alias] = trim(name);
            }
            continue;
        }
        const std::string code = line.substr(0, 3);
        if (code != "317" && code != "327" && code != "367") continue;
        line.resize(std::max<std::size_t>(line.size(), 80), ' ');
        NetlistPoint p;
        p.net = trim(line.substr(3, 14));
        if (auto a = aliases.find(p.net); a != aliases.end()) p.net = a->second;
        p.refdes = trim(line.substr(20, 6));
        p.pin = trim(line.substr(27, 4));
        p.midpoint = line[31] == 'M';
        // Fields from column 33 on are letter-tagged: D hole, P/U plating, A access,
        // X/Y position, X/Y size, R rotation, S solder mask.
        const std::string tail = line.substr(32);
        int xy = 0;
        std::size_t i = 0;
        while (i < tail.size()) {
            const char c = tail[i];
            std::size_t j = i + 1;
            while (j < tail.size() && (std::isdigit(static_cast<unsigned char>(tail[j])) || tail[j] == '+' || tail[j] == '-')) ++j;
            const std::string num = tail.substr(i + 1, j - i - 1);
            try {
                if (c == 'D' && !num.empty()) p.hole_nm = std::stod(num) * unit_nm;
                else if (c == 'P') p.plated = true;
                else if (c == 'U') p.plated = false;
                else if (c == 'A' && !num.empty()) p.access = std::stoi(num);
                else if (c == 'X' && !num.empty() && xy == 0) p.at.x = std::llround(std::stod(num) * unit_nm);
                else if (c == 'Y' && !num.empty() && xy == 0) {
                    p.at.y = std::llround(std::stod(num) * unit_nm);
                    xy = 1;  // later X/Y fields are the pad size
                }
            } catch (const std::exception&) {
                nl.warnings.push_back("malformed IPC-D-356 record: " + trim(line));
                break;
            }
            i = j;
        }
        nl.points.push_back(std::move(p));
    }
    return nl;
}

// ---------------------------------------------------------------- placement CSV

namespace {

std::vector<std::string> csv_fields(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (quoted) {
            if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else if (c == '"') {
                quoted = false;
            } else {
                cur.push_back(c);
            }
        } else if (c == '"') {
            quoted = true;
        } else if (c == ',') {
            out.push_back(trim(cur));
            cur.clear();
        } else if (c != '\r') {
            cur.push_back(c);
        }
    }
    out.push_back(trim(cur));
    return out;
}

int find_column(const std::vector<std::string>& header, std::initializer_list<const char*> names) {
    for (const char* n : names)
        for (std::size_t i = 0; i < header.size(); ++i)
            if (lower(header[i]) == n) return static_cast<int>(i);
    for (const char* n : names)
        for (std::size_t i = 0; i < header.size(); ++i)
            if (lower(header[i]).rfind(n, 0) == 0) return static_cast<int>(i);
    return -1;
}

double length_value(std::string v, double column_unit_nm) {
    v = lower(trim(v));
    double unit = column_unit_nm;
    if (v.size() > 3 && v.substr(v.size() - 3) == "mil") {
        unit = 25400.0;
        v.resize(v.size() - 3);
    } else if (v.size() > 2 && v.substr(v.size() - 2) == "mm") {
        unit = kNmPerMm;
        v.resize(v.size() - 2);
    } else if (v.size() > 2 && v.substr(v.size() - 2) == "in") {
        unit = 25.4 * kNmPerMm;
        v.resize(v.size() - 2);
    }
    return std::stod(v) * unit;
}

}  // namespace

bool looks_like_placement_csv(std::string_view text) {
    const std::string first = lower(std::string(text.substr(0, text.find('\n'))));
    return first.find(',') != std::string::npos &&
           (first.find("ref") != std::string::npos || first.find("designator") != std::string::npos) &&
           (first.find("pos") != std::string::npos || first.find("x") != std::string::npos);
}

std::vector<PlacementRow> parse_placement_csv(std::string_view text) {
    std::istringstream in{std::string(text)};
    std::string line;
    std::vector<std::string> header;
    while (header.empty() && std::getline(in, line))
        if (!trim(line).empty()) header = csv_fields(line);
    const int ref = find_column(header, {"ref", "designator", "refdes", "reference"});
    const int val = find_column(header, {"val", "value", "comment"});
    const int pkg = find_column(header, {"package", "footprint"});
    const int x = find_column(header, {"posx", "mid x", "center-x", "x"});
    const int y = find_column(header, {"posy", "mid y", "center-y", "y"});
    const int rot = find_column(header, {"rot", "rotation"});
    const int side = find_column(header, {"side", "layer"});
    if (ref < 0 || x < 0 || y < 0) throw std::runtime_error("placement CSV: missing reference or position column");
    auto unit_of = [&](int col) {
        const std::string h = lower(header[col]);
        if (h.find("mil") != std::string::npos) return 25400.0;
        if (h.find("(in") != std::string::npos) return 25.4 * kNmPerMm;
        return kNmPerMm;
    };
    std::vector<PlacementRow> rows;
    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        const auto f = csv_fields(line);
        auto at = [&](int col) { return col >= 0 && col < static_cast<int>(f.size()) ? f[col] : std::string(); };
        PlacementRow r;
        r.refdes = at(ref);
        if (r.refdes.empty()) continue;
        r.value = at(val);
        r.package = at(pkg);
        r.at = Point64(std::llround(length_value(at(x), unit_of(x))), std::llround(length_value(at(y), unit_of(y))));
        if (rot >= 0 && !at(rot).empty()) r.rotation_deg = std::stod(at(rot));
        const std::string s = lower(at(side));
        r.side = (s.rfind("b", 0) == 0) ? Side::Bottom : Side::Top;
        rows.push_back(std::move(r));
    }
    return rows;
}

// ---------------------------------------------------------------- stackup YAML

std::vector<StackupLayer> parse_stackup_yaml(std::string_view text) {
    const YAML::Node root = YAML::Load(std::string(text));
    const YAML::Node layers = root["layers"];
    if (!layers || !layers.IsSequence()) throw std::runtime_error("stackup YAML: missing 'layers' list");
    std::vector<StackupLayer> out;
    int copper = 0;
    for (const auto& n : layers) {
        StackupLayer l;
        l.name = n["name"] ? n["name"].as<std::string>() : std::string();
        l.type = layer_type_from_string(n["type"] ? n["type"].as<std::string>() : std::string());
        if (!n["thickness_mm"]) throw std::runtime_error("stackup YAML: layer '" + l.name + "' has no thickness_mm");
        l.thickness_nm = n["thickness_mm"].as<double>() * kNmPerMm;
        if (n["material"]) l.material = n["material"].as<std::string>();
        if (n["epsilon_r"]) l.epsilon_r = n["epsilon_r"].as<double>();
        if (n["loss_tangent"]) l.loss_tangent = n["loss_tangent"].as<double>();
        if (n["conductivity_s_per_m"]) l.conductivity_s_per_m = n["conductivity_s_per_m"].as<double>();
        if (l.type == LayerType::Copper) l.copper_index = ++copper;
        out.push_back(std::move(l));
    }
    return out;
}

}  // namespace kairoq
