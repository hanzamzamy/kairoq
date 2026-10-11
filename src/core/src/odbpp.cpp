// ODB++ importer: the matrix (layer order and types), layer attributes (stackup), feature
// files (copper, solder mask, drills, outline), the EDA data (nets, pads) and the
// component files. Format: ODB++ Design Format Specification 8.1 (Siemens).
#include "kairoq/odbpp.hpp"

#include "kairoq/archive.hpp"
#include "kairoq/feature.hpp"

#include "build.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <future>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace fs = std::filesystem;

std::vector<std::pair<std::string, std::string>> read_tar(const std::string& data) {
    std::vector<std::pair<std::string, std::string>> out;
    std::string long_name;
    std::size_t at = 0;
    auto field = [&](std::size_t off, std::size_t len) {
        std::string s = data.substr(at + off, len);
        const auto z = s.find('\0');
        return z == std::string::npos ? s : s.substr(0, z);
    };
    while (at + 512 <= data.size()) {
        if (data.compare(at, 512, std::string(512, '\0')) == 0) break;
        std::string name = field(0, 100);
        const std::string prefix = field(345, 155);
        if (data.compare(at + 257, 5, "ustar") == 0 && !prefix.empty()) name = prefix + "/" + name;
        const std::string size_text = field(124, 12);
        std::size_t size = 0;
        for (char c : size_text)
            if (c >= '0' && c <= '7') size = size * 8 + static_cast<std::size_t>(c - '0');
        const char type = data[at + 156];
        at += 512;
        if (at + size > data.size()) throw std::runtime_error("truncated tar archive");
        const std::string body = data.substr(at, size);
        at += (size + 511) / 512 * 512;
        if (type == 'L') {  // GNU long name
            long_name = body.substr(0, body.find('\0'));
            continue;
        }
        if (type == 'x') {  // POSIX extended header: records "length key=value\n"
            for (std::size_t r = 0; r < body.size();) {
                const auto sp = body.find(' ', r);
                if (sp == std::string::npos) break;
                const std::size_t len = std::stoul(body.substr(r, sp - r));
                if (len == 0 || r + len > body.size()) break;
                const std::string rec = body.substr(sp + 1, r + len - sp - 2);
                if (rec.rfind("path=", 0) == 0) long_name = rec.substr(5);
                r += len;
            }
            continue;
        }
        if (!long_name.empty()) {
            name = long_name;
            long_name.clear();
        }
        if (type == '0' || type == '\0') out.emplace_back(name, body);
    }
    return out;
}

std::string uncompress_lzw(const std::string& data) {
    // Unix compress: 9- to maxbits-bit codes, least significant bit first, in groups of
    // eight codes; when the code width changes (or the table is cleared) the rest of the
    // current group is skipped.
    const auto* in = reinterpret_cast<const unsigned char*>(data.data());
    const std::size_t n = data.size();
    if (n < 3 || in[0] != 0x1f || in[1] != 0x9d) throw std::runtime_error("not compress (.Z) data");
    const int maxbits = in[2] & 0x1f;
    const bool block = (in[2] & 0x80) != 0;
    if (maxbits < 9 || maxbits > 16) throw std::runtime_error("unsupported compress code size");
    std::string out;
    if (n == 3) return out;
    std::vector<std::uint16_t> prefix(1u << maxbits);
    std::vector<unsigned char> suffix(1u << maxbits);
    std::vector<unsigned char> stack;
    int bits = 9;
    unsigned mask = 0x1ff;
    unsigned end = block ? 256 : 255;
    std::size_t bitpos = 0, mark = 0;  // bits from the start of the code data
    const std::size_t total = (n - 3) * 8;
    auto read = [&](unsigned& code) {
        if (bitpos + static_cast<std::size_t>(bits) > total) return false;
        code = 0;
        for (int k = 0; k < bits; ++k) {
            const std::size_t b = bitpos + static_cast<std::size_t>(k);
            if (in[3 + b / 8] & (1u << (b % 8))) code |= 1u << k;
        }
        bitpos += static_cast<std::size_t>(bits);
        return true;
    };
    auto skip_group = [&]() {
        const std::size_t group = static_cast<std::size_t>(bits) * 8;
        const std::size_t used = (bitpos - mark) % group;
        if (used) bitpos += group - used;
        mark = bitpos;
    };
    unsigned code = 0;
    if (!read(code)) return out;
    if (code > 255) throw std::runtime_error("corrupt compress data");
    unsigned prev = code;
    unsigned char final = static_cast<unsigned char>(code);
    out.push_back(static_cast<char>(final));
    for (;;) {
        if (end >= mask && bits < maxbits) {
            skip_group();
            ++bits;
            mask = (mask << 1) | 1u;
        }
        if (!read(code)) break;
        if (code == 256 && block) {
            skip_group();
            bits = 9;
            mask = 0x1ff;
            end = 255;
            if (!read(code)) break;
            if (code > 255) throw std::runtime_error("corrupt compress data");
            prev = code;
            final = static_cast<unsigned char>(code);
            out.push_back(static_cast<char>(final));
            end = 256;
            continue;
        }
        const unsigned temp = code;
        stack.clear();
        if (code > end) {
            if (code != end + 1 || prev > end) throw std::runtime_error("corrupt compress data");
            stack.push_back(final);
            code = prev;
        }
        while (code >= 256) {
            stack.push_back(suffix[code]);
            code = prefix[code];
        }
        stack.push_back(static_cast<unsigned char>(code));
        final = static_cast<unsigned char>(code);
        if (end < mask) {
            ++end;
            prefix[end] = static_cast<std::uint16_t>(prev);
            suffix[end] = final;
        }
        prev = temp;
        for (auto it = stack.rbegin(); it != stack.rend(); ++it) out.push_back(static_cast<char>(*it));
    }
    return out;
}

namespace {

std::int64_t ll(double v) { return static_cast<std::int64_t>(std::llround(v)); }

Point64 pt(std::int64_t x, std::int64_t y) { return Point64(x, y); }

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r") - a + 1);
}

// The files of a job by their path below the job's root, in lower case.
class JobFiles {
public:
    explicit JobFiles(const std::string& path) {
        if (fs::is_directory(path)) {
            for (const auto& e : fs::recursive_directory_iterator(path))
                if (e.is_regular_file()) disk_[lower(fs::relative(e.path(), path).generic_string())] = e.path().string();
        } else {
            std::vector<std::pair<std::string, std::string>> entries;
            const std::string raw = read_file(path);
            if (is_zip_file(path)) {
                for (auto& e : read_zip(path)) entries.emplace_back(e.name, std::move(e.data));
            } else if (raw.size() > 2 && static_cast<unsigned char>(raw[0]) == 0x1f && static_cast<unsigned char>(raw[1]) == 0x8b) {
                entries = read_tar(gunzip(raw));
            } else {
                entries = read_tar(raw);
            }
            for (auto& [name, data] : entries) mem_[lower(name)] = std::move(data);
        }
        // The root is where matrix/matrix lies.
        auto find_root = [&](const auto& files) {
            for (const auto& [name, x] : files)
                for (const char* m : {"matrix/matrix", "matrix/matrix.z"})
                    if (name.size() >= std::strlen(m) && name.compare(name.size() - std::strlen(m), std::strlen(m), m) == 0) {
                        root_ = name.substr(0, name.size() - std::strlen(m));
                        return true;
                    }
            return false;
        };
        found_ = disk_.empty() ? find_root(mem_) : find_root(disk_);
    }
    bool found() const { return found_; }
    // Content of a file, decompressed when stored as .z.
    std::optional<std::string> get(const std::string& rel) const {
        const std::string key = root_ + lower(rel);
        for (const auto& k : {key, key + ".z"}) {
            std::optional<std::string> data;
            if (auto it = mem_.find(k); it != mem_.end()) data = it->second;
            else if (auto jt = disk_.find(k); jt != disk_.end()) data = read_file(jt->second);
            if (!data) continue;
            if (k.size() > 2 && k.compare(k.size() - 2, 2, ".z") == 0) return uncompress_lzw(*data);
            return data;
        }
        return std::nullopt;
    }
    // Directory names below a path.
    std::vector<std::string> subdirs(const std::string& rel) const {
        std::set<std::string> out;
        const std::string pre = root_ + lower(rel) + "/";
        auto scan = [&](const auto& files) {
            for (const auto& [name, x] : files)
                if (name.compare(0, pre.size(), pre) == 0) {
                    const auto slash = name.find('/', pre.size());
                    if (slash != std::string::npos) out.insert(name.substr(pre.size(), slash - pre.size()));
                }
        };
        scan(mem_);
        scan(disk_);
        return {out.begin(), out.end()};
    }

private:
    std::map<std::string, std::string> mem_, disk_;
    std::string root_;
    bool found_ = false;
};

// Structured text: KEY=VALUE lines and NAME { ... } blocks.
struct Block {
    std::string name;
    std::map<std::string, std::string> values;
};
std::vector<Block> parse_blocks(const std::string& text, std::map<std::string, std::string>* top = nullptr) {
    std::vector<Block> out;
    std::istringstream in(text);
    std::string line;
    Block* cur = nullptr;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line.back() == '{') {
            out.push_back({upper(trim(line.substr(0, line.size() - 1))), {}});
            cur = &out.back();
            continue;
        }
        if (line == "}") {
            cur = nullptr;
            continue;
        }
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string k = upper(trim(line.substr(0, eq))), v = trim(line.substr(eq + 1));
        if (cur) cur->values[k] = v;
        else if (top) (*top)[k] = v;
    }
    return out;
}

std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string w;
    while (in >> w) out.push_back(w);
    return out;
}

// One parsed feature record.
struct Rec {
    char type = 0;  // L P A S
    std::vector<double> v;
    int symbol = -1;
    bool positive = true;
    int orient = 0;
    double rotation = 0;
    bool cw = false;
    std::map<int, std::string> attrs;
    std::vector<std::pair<bool, Paths64>> contours;  // surface: (island, path)
};

struct FeatureFile {
    double unit_nm = kNmPerMm;         // coordinates
    double symbol_unit_nm = 1000.0;    // symbol sizes: micrometres or mils
    std::vector<std::string> symbols;  // $n
    std::map<int, std::string> attr_names, attr_texts;
    std::vector<Rec> recs;
};

class Importer {
public:
    Importer(const ImportOptions& o) : tol_(o.geometry.arc_tolerance_nm) {}

    Board run(const std::string& path) {
        JobFiles files(path);
        if (!files.found()) throw std::runtime_error(path + ": no ODB++ job (matrix/matrix)");
        files_ = &files;
        board_.nets.push_back("");
        board_.source.importer = "odbpp";
        board_.source.files.push_back(fs::path(path).filename().string());
        board_.source.arc_tolerance_nm = tol_;
        std::map<std::string, std::string> info;
        if (auto t = files.get("misc/info")) parse_blocks(*t, &info);
        if (info.count("UNITS") && upper(info["UNITS"]) == "INCH") job_unit_nm_ = 25.4 * kNmPerMm;
        board_.source.generator = info.count("SAVE_APP") ? info["SAVE_APP"] : (info.count("ODB_SOURCE") ? info["ODB_SOURCE"] : "");
        board_.name = fs::path(path).stem().string();  // board.odb.zip, board.tar.gz: "board"
        for (const char* ext : {".odb", ".tar"})
            if (board_.name.size() > 4 && lower(board_.name.substr(board_.name.size() - 4)) == ext) board_.name.resize(board_.name.size() - 4);

        const auto matrix = files.get("matrix/matrix");
        const auto blocks = parse_blocks(*matrix);
        for (const auto& b : blocks)
            if (b.name == "STEP" && step_.empty()) step_ = lower(b.values.count("NAME") ? b.values.at("NAME") : "");
        if (step_.empty()) step_ = files.subdirs("steps").empty() ? "" : files.subdirs("steps").front();
        if (step_.empty()) throw std::runtime_error(path + ": no step");
        for (const auto& b : blocks)
            if (b.name == "LAYER") {
                Layer l;
                l.row = b.values.count("ROW") ? std::stoi(b.values.at("ROW")) : 0;
                l.context = upper(b.values.count("CONTEXT") ? b.values.at("CONTEXT") : "");
                l.type = upper(b.values.count("TYPE") ? b.values.at("TYPE") : "");
                l.name = b.values.count("NAME") ? b.values.at("NAME") : "";
                l.negative = upper(b.values.count("POLARITY") ? b.values.at("POLARITY") : "") == "NEGATIVE";
                l.start = b.values.count("START_NAME") ? b.values.at("START_NAME") : "";
                l.end = b.values.count("END_NAME") ? b.values.at("END_NAME") : "";
                layers_.push_back(l);
            }
        std::sort(layers_.begin(), layers_.end(), [](const Layer& a, const Layer& c) { return a.row < c.row; });
        read_stackup();
        read_eda();
        read_components();
        read_outline();
        read_layers();
        detail::assign_drill_nets(board_);
        detail::check_net_overlaps(board_);
        board_.net_source = "odbpp";
        return std::move(board_);
    }

private:
    struct Layer {
        int row = 0;
        std::string context, type, name, start, end;
        bool negative = false;
    };
    struct FeatureNet {
        int net = -1;
        int side = 0;  // 1 top, 2 bottom component; 0: not a pad
        int comp = -1, pin = -1;
    };
    struct Comp {
        std::string refdes, package, value;
        Point64 at{0, 0};
        double rotation = 0;
        bool mirror = false;
        std::vector<std::string> pins;
    };

    const JobFiles* files_ = nullptr;
    double tol_;
    double job_unit_nm_ = kNmPerMm;
    std::string step_;
    std::vector<Layer> layers_;
    std::map<std::string, int> copper_index_;  // lower-case layer name -> 1..N
    std::vector<std::string> net_names_;       // EDA net index -> name
    std::map<std::pair<std::string, int>, FeatureNet> feature_net_;  // (layer, feature index)
    std::map<int, std::vector<Comp>> comps_;   // side -> components
    Board board_;
    std::set<std::string> warned_;

    void warn(const std::string& m) {
        if (warned_.insert(m).second && board_.warnings.size() < 100) board_.warnings.push_back(m);
    }

    std::string layer_path(const std::string& name, const std::string& file) const {
        return "steps/" + step_ + "/layers/" + lower(name) + "/" + file;
    }

    static bool is_copper(const Layer& l) {
        return l.context == "BOARD" && (l.type == "SIGNAL" || l.type == "POWER_GROUND" || l.type == "MIXED");
    }

    std::map<std::string, std::string> attrlist(const std::string& layer) const {
        std::map<std::string, std::string> out;
        if (auto t = files_->get(layer_path(layer, "attrlist"))) {
            std::istringstream in(*t);
            std::string line;
            while (std::getline(in, line)) {
                line = trim(line);
                const auto eq = line.find('=');
                if (line.empty() || line[0] == '#' || eq == std::string::npos) continue;
                out[lower(trim(line.substr(0, eq)))] = trim(line.substr(eq + 1));
            }
        }
        return out;
    }

    // ---------------------------------------------------------------- stackup

    void read_stackup() {
        int n = 0;
        bool any_dielectric = false;
        for (const auto& l : layers_) {
            if (l.context != "BOARD") continue;
            if (is_copper(l)) copper_index_[lower(l.name)] = ++n;
            any_dielectric = any_dielectric || l.type == "DIELECTRIC";
        }
        if (n == 0) throw std::runtime_error("no copper layer in the matrix");
        bool seen_copper = false;
        for (const auto& l : layers_) {
            if (l.context != "BOARD") continue;
            StackupLayer s;
            s.name = l.name;
            const auto a = attrlist(l.name);
            auto num = [&](const char* k) -> std::optional<double> {
                auto it = a.find(k);
                if (it == a.end()) return std::nullopt;
                try {
                    return std::stod(it->second);
                } catch (...) {
                    return std::nullopt;
                }
            };
            if (auto t = num(".layer_dielectric")) s.thickness_nm = *t * job_unit_nm_;
            if (is_copper(l)) {
                s.type = LayerType::Copper;
                s.copper_index = copper_index_.at(lower(l.name));
                if (s.thickness_nm <= 0) s.thickness_nm = num(".copper_weight").value_or(1.0) * 35000.0;
                s.material = "copper";
                seen_copper = true;
                CopperLayer cl;
                cl.index = s.copper_index;
                cl.name = l.name;
                cl.side = s.copper_index == 1 ? Side::Top : (s.copper_index == n ? Side::Bottom : Side::Inner);
                board_.copper.push_back(cl);
            } else if (l.type == "DIELECTRIC") {
                s.type = LayerType::Dielectric;
                if (auto it = a.find(".material"); it != a.end()) s.material = it->second;
                s.epsilon_r = num(".dielectric_constant");
                s.loss_tangent = num(".loss_tangent");
            } else if (l.type == "SOLDER_MASK") {
                s.type = LayerType::SolderMask;
                s.epsilon_r = num(".dielectric_constant");
                if (s.thickness_nm <= 0) s.thickness_nm = 10000.0;
                (void)seen_copper;
            } else {
                continue;
            }
            board_.stackup.push_back(s);
        }
        bool thick = any_dielectric;
        for (const auto& s : board_.stackup)
            if (s.type == LayerType::Dielectric && s.thickness_nm <= 0) thick = false;
        if (thick) {
            board_.stackup_source = "odbpp";
        } else {
            warn("the job gives no dielectric thicknesses; an even 1.6 mm FR-4 stackup is assumed");
            board_.stackup = detail::default_stackup(n, 1.6e6);
            board_.stackup_source = "default";
        }
    }

    // ---------------------------------------------------------------- nets, components

    void read_eda() {
        const auto t = files_->get("steps/" + step_ + "/eda/data");
        if (!t) {
            warn("no EDA data: copper without nets");
            return;
        }
        std::istringstream in(*t);
        std::string line;
        std::vector<std::string> lyr;
        int net = -1;
        FeatureNet cur;
        while (std::getline(in, line)) {
            const auto w = split(line);
            if (w.empty() || w[0][0] == '#') continue;
            if (w[0] == "LYR") {
                for (std::size_t i = 1; i < w.size(); ++i) lyr.push_back(lower(w[i]));
            } else if (w[0] == "NET") {
                std::string name = w.size() > 1 ? w[1] : "";
                if (const auto semi = name.find(';'); semi != std::string::npos) name = name.substr(0, semi);
                if (name == "$NONE$") name.clear();
                net_names_.push_back(name);
                net = static_cast<int>(net_names_.size()) - 1;
                cur = FeatureNet{net};
            } else if (w[0] == "SNT") {
                cur = FeatureNet{net};
                if (w.size() >= 5 && w[1] == "TOP") {
                    cur.side = w[2] == "B" ? 2 : 1;
                    cur.comp = std::stoi(w[3]);
                    cur.pin = std::stoi(w[4]);
                }
            } else if (w[0] == "FID" && w.size() >= 4) {
                const int li = std::stoi(w[2]), fi = std::stoi(w[3]);
                if (li >= 0 && static_cast<std::size_t>(li) < lyr.size()) feature_net_[{lyr[static_cast<std::size_t>(li)], fi}] = cur;
            }
        }
    }

    void read_components() {
        for (const auto& [side, dir] : {std::pair{1, "comp_+_top"}, std::pair{2, "comp_+_bot"}}) {
            const auto t = files_->get(layer_path(dir, "components"));
            if (!t) continue;
            double unit = job_unit_nm_;
            std::istringstream in(*t);
            std::string line;
            auto& list = comps_[side];
            while (std::getline(in, line)) {
                const std::string tl = trim(line);
                if (tl.rfind("UNITS=", 0) == 0) unit = upper(tl.substr(6)) == "INCH" ? 25.4 * kNmPerMm : kNmPerMm;
                const auto w = split(tl.substr(0, tl.find(';')));
                if (w.empty()) continue;
                if (w[0] == "CMP" && w.size() >= 7) {
                    Comp c;
                    c.at = pt(ll(std::stod(w[2]) * unit), ll(std::stod(w[3]) * unit));
                    c.rotation = std::stod(w[4]);
                    c.mirror = w[5] == "M";
                    c.refdes = w[6];
                    c.package = w.size() > 7 ? w[7] : "";
                    list.push_back(std::move(c));
                } else if (w[0] == "PRP" && w.size() >= 3 && w[1] == "Value" && !list.empty()) {
                    std::string v = tl.substr(tl.find("Value") + 5);
                    v = trim(v);
                    if (v.size() >= 2 && v.front() == '\'' && v.back() == '\'') v = v.substr(1, v.size() - 2);
                    list.back().value = v;
                } else if (w[0] == "TOP" && w.size() >= 9 && !list.empty()) {
                    list.back().pins.push_back(w[8]);
                }
            }
            for (const auto& c : list) {
                Component comp;
                comp.refdes = c.refdes;
                comp.package = c.package;
                comp.value = c.value;
                comp.at = c.at;
                // ODB++ angles are clockwise; the model's are counter-clockwise.
                comp.rotation_deg = c.rotation == 0 ? 0.0 : 360.0 - std::fmod(c.rotation, 360.0);
                comp.side = side == 2 ? Side::Bottom : Side::Top;
                comp.position_source = "odbpp";
                board_.components.push_back(std::move(comp));
            }
        }
    }

    // ---------------------------------------------------------------- feature files

    FeatureFile parse_features(const std::string& text) const {
        FeatureFile f;
        f.unit_nm = job_unit_nm_;
        std::istringstream in(text);
        std::string line;
        Rec* surface = nullptr;
        Path64 contour;
        bool island = true;
        Point64 last{0, 0};
        auto P = [&](double x, double y) { return pt(ll(x * f.unit_nm), ll(y * f.unit_nm)); };
        while (std::getline(in, line)) {
            std::string tl = trim(line);
            if (tl.empty() || tl[0] == '#') continue;
            if (tl.rfind("UNITS=", 0) == 0) {
                const bool inch = upper(tl.substr(6)) == "INCH";
                f.unit_nm = inch ? 25.4 * kNmPerMm : kNmPerMm;
                f.symbol_unit_nm = inch ? 25.4 * 1000.0 : 1000.0;
                continue;
            }
            if (tl.rfind("U ", 0) == 0) {  // older files: U MM / U INCH
                const bool inch = upper(trim(tl.substr(2))) == "INCH";
                f.unit_nm = inch ? 25.4 * kNmPerMm : kNmPerMm;
                f.symbol_unit_nm = inch ? 25.4 * 1000.0 : 1000.0;
                continue;
            }
            if (tl[0] == '$') {
                const auto sp = tl.find(' ');
                const int idx = std::stoi(tl.substr(1, sp - 1));
                if (static_cast<int>(f.symbols.size()) <= idx) f.symbols.resize(static_cast<std::size_t>(idx + 1));
                f.symbols[static_cast<std::size_t>(idx)] = trim(tl.substr(sp + 1));
                continue;
            }
            if (tl[0] == '@' || tl[0] == '&') {
                const auto sp = tl.find(' ');
                const int idx = std::stoi(tl.substr(1, sp - 1));
                (tl[0] == '@' ? f.attr_names : f.attr_texts)[idx] = trim(tl.substr(sp + 1));
                continue;
            }
            std::map<int, std::string> attrs;
            if (const auto semi = tl.find(';'); semi != std::string::npos) {
                std::string a = tl.substr(semi + 1);
                tl = trim(tl.substr(0, semi));
                std::stringstream as(a);
                std::string item;
                while (std::getline(as, item, ',')) {
                    item = trim(item);
                    if (item.empty()) continue;
                    const auto eq = item.find('=');
                    try {
                        attrs[std::stoi(item.substr(0, eq))] = eq == std::string::npos ? "" : item.substr(eq + 1);
                    } catch (...) {
                    }
                }
            }
            const auto w = split(tl);
            const std::string& k = w[0];
            if (surface) {
                if (k == "OB" && w.size() >= 4) {
                    contour = {P(std::stod(w[1]), std::stod(w[2]))};
                    last = contour.back();
                    island = w[3] != "H";
                } else if (k == "OS" && w.size() >= 3) {
                    last = P(std::stod(w[1]), std::stod(w[2]));
                    contour.push_back(last);
                } else if (k == "OC" && w.size() >= 6) {
                    const Point64 e = P(std::stod(w[1]), std::stod(w[2])), c = P(std::stod(w[3]), std::stod(w[4]));
                    const bool cw = w[5] == "Y";
                    const Path64 arc = arc_points(last, e, c, !cw, tol_, true);
                    contour.insert(contour.end(), arc.begin() + 1, arc.end());
                    last = e;
                } else if (k == "OE") {
                    if (contour.size() >= 3) surface->contours.push_back({island, Paths64{contour}});
                    contour.clear();
                } else if (k == "SE") {
                    surface = nullptr;
                }
                continue;
            }
            Rec r;
            r.type = k[0];
            r.attrs = std::move(attrs);
            auto num = [&](std::size_t i) { return i < w.size() ? std::stod(w[i]) : 0.0; };
            if (k == "L" && w.size() >= 7) {
                r.v = {num(1), num(2), num(3), num(4)};
                r.symbol = std::stoi(w[5]);
                r.positive = w[6] != "N";
            } else if (k == "P" && w.size() >= 6) {
                r.v = {num(1), num(2)};
                r.symbol = std::stoi(w[3]);
                r.positive = w[4] != "N";
                r.orient = w.size() > 6 ? std::stoi(w[6]) : 0;
                if ((r.orient == 8 || r.orient == 9) && w.size() > 7) r.rotation = num(7);
            } else if (k == "A" && w.size() >= 10) {
                r.v = {num(1), num(2), num(3), num(4), num(5), num(6)};
                r.symbol = std::stoi(w[7]);
                r.positive = w[8] != "N";
                r.cw = w.size() > 10 && w[10] == "Y";
            } else if (k == "S" && w.size() >= 2) {
                r.positive = w[1] != "N";
                f.recs.push_back(std::move(r));
                surface = &f.recs.back();
                continue;
            } else if (k == "T" || k == "B") {
                // Text and barcodes carry no copper the analyses use; they keep their index.
                r.type = 'X';
            } else {
                continue;  // not a feature record (the count line F, identifiers)
            }
            f.recs.push_back(std::move(r));
        }
        return f;
    }

    // Shape of a standard symbol at the origin (nm).
    Paths64 symbol_shape(const std::string& name, double unit) {
        auto nums = [&](const std::string& s, char sep = 'x') {
            std::vector<double> v;
            std::stringstream ss(s);
            std::string part;
            while (std::getline(ss, part, sep)) {
                try {
                    v.push_back(std::stod(part) * unit);
                } catch (...) {
                    v.push_back(0);
                }
            }
            return v;
        };
        auto rect = [](double w, double h) {
            const auto x = ll(w / 2), y = ll(h / 2);
            return Path64{pt(-x, -y), pt(x, -y), pt(x, y), pt(-x, y)};
        };
        auto starts = [&](const char* p) { return name.rfind(p, 0) == 0; };
        if (starts("rect")) {
            std::string rest = name.substr(4);
            double r = 0, c = 0;
            if (auto p = rest.find("xr"); p != std::string::npos) {
                r = nums(rest.substr(p + 2))[0];
                rest = rest.substr(0, p);
            } else if (auto q = rest.find("xc"); q != std::string::npos) {
                c = nums(rest.substr(q + 2))[0];
                rest = rest.substr(0, q);
            }
            const auto wh = nums(rest);
            if (wh.size() < 2) return {};
            const double w = wh[0], h = wh[1];
            if (r > 0) {
                r = std::min({r, w / 2, h / 2});
                return to_paths(union_polygons(Clipper2Lib::InflatePaths(Paths64{rect(w - 2 * r, h - 2 * r)}, r, Clipper2Lib::JoinType::Round,
                                                                         Clipper2Lib::EndType::Polygon, 2.0, tol_)));
            }
            if (c > 0) {
                const auto x = w / 2, y = h / 2;
                return {Path64{pt(ll(-x + c), ll(-y)), pt(ll(x - c), ll(-y)),
                               pt(ll(x), ll(-y + c)), pt(ll(x), ll(y - c)),
                               pt(ll(x - c), ll(y)), pt(ll(-x + c), ll(y)),
                               pt(ll(-x), ll(y - c)), pt(ll(-x), ll(-y + c))}};
            }
            return {rect(w, h)};
        }
        if (starts("r") && name.size() > 1 && (std::isdigit(static_cast<unsigned char>(name[1])) || name[1] == '.'))
            return {circle_path(pt(std::int64_t{0}, 0), nums(name.substr(1))[0] / 2, tol_)};
        if (starts("s") && name.size() > 1 && (std::isdigit(static_cast<unsigned char>(name[1])) || name[1] == '.')) {
            const double s = nums(name.substr(1))[0];
            return {rect(s, s)};
        }
        if (starts("oval")) {
            const auto wh = nums(name.substr(4));
            if (wh.size() < 2) return {};
            const double w = wh[0], h = wh[1], d = std::min(w, h), half = (std::max(w, h) - d) / 2;
            const Point64 a = w >= h ? pt(ll(-half), std::int64_t{0}) : pt(std::int64_t{0}, ll(-half));
            const Point64 b = w >= h ? pt(ll(half), std::int64_t{0}) : pt(std::int64_t{0}, ll(half));
            return {capsule(a, b, d / 2, tol_)};
        }
        if (starts("di")) {
            const auto wh = nums(name.substr(2));
            if (wh.size() < 2) return {};
            const auto x = ll(wh[0] / 2), y = ll(wh[1] / 2);
            return {Path64{pt(-x, 0), pt(std::int64_t{0}, -y), pt(x, 0), pt(std::int64_t{0}, y)}};
        }
        if (starts("oct")) {
            const auto v = nums(name.substr(3));
            if (v.size() < 3) return {};
            const double x = v[0] / 2, y = v[1] / 2, c = v[2];
            return {Path64{pt(ll(-x + c), ll(-y)), pt(ll(x - c), ll(-y)),
                           pt(ll(x), ll(-y + c)), pt(ll(x), ll(y - c)),
                           pt(ll(x - c), ll(y)), pt(ll(-x + c), ll(y)),
                           pt(ll(-x), ll(y - c)), pt(ll(-x), ll(-y + c))}};
        }
        if (starts("donut_r")) {
            const auto v = nums(name.substr(7));
            if (v.size() < 2) return {};
            Path64 inner = circle_path(pt(std::int64_t{0}, 0), v[1] / 2, tol_);
            std::reverse(inner.begin(), inner.end());
            return {circle_path(pt(std::int64_t{0}, 0), v[0] / 2, tol_), inner};
        }
        if (starts("donut_s")) {
            const auto v = nums(name.substr(7));
            if (v.size() < 2) return {};
            Path64 inner = rect(v[1], v[1]);
            std::reverse(inner.begin(), inner.end());
            return {rect(v[0], v[0]), inner};
        }
        if (starts("el")) {
            const auto wh = nums(name.substr(2));
            if (wh.size() < 2) return {};
            const int k = std::max(16, circle_segments(std::max(wh[0], wh[1]) / 2, tol_));
            Path64 p;
            for (int i = 0; i < k; ++i)
                p.emplace_back(ll(wh[0] / 2 * std::cos(2 * kPi * i / k)), ll(wh[1] / 2 * std::sin(2 * kPi * i / k)));
            return {p};
        }
        if (starts("hex_l")) {
            const auto v = nums(name.substr(5));
            if (v.size() < 3) return {};
            const double x = v[0] / 2, y = v[1] / 2, c = v[2];
            return {Path64{pt(ll(-x), std::int64_t{0}), pt(ll(-x + c), ll(-y)), pt(ll(x - c), ll(-y)),
                           pt(ll(x), std::int64_t{0}), pt(ll(x - c), ll(y)), pt(ll(-x + c), ll(y))}};
        }
        if (starts("hex_s")) {
            const auto v = nums(name.substr(5));
            if (v.size() < 3) return {};
            const double x = v[0] / 2, y = v[1] / 2, c = v[2];
            return {Path64{pt(std::int64_t{0}, ll(-y)), pt(ll(x), ll(-y + c)), pt(ll(x), ll(y - c)),
                           pt(std::int64_t{0}, ll(y)), pt(ll(-x), ll(y - c)), pt(ll(-x), ll(-y + c))}};
        }
        if (starts("tri")) {
            const auto v = nums(name.substr(3));
            if (v.size() < 2) return {};
            const auto b = ll(v[0] / 2), h = ll(v[1] / 2);
            return {Path64{pt(-b, -h), pt(b, -h), pt(std::int64_t{0}, h)}};
        }
        // A user symbol: its own feature file, drawn positive.
        if (auto t = files_->get("symbols/" + lower(name) + "/features")) {
            const FeatureFile uf = parse_features(*t);
            std::vector<Feature> feats;
            for (std::size_t i = 0; i < uf.recs.size(); ++i) to_features(uf, uf.recs[i], "", static_cast<int>(i), feats);
            return to_paths(render_all(feats));
        }
        warn("symbol " + name + " is not supported; drawn as its nominal circle");
        const auto v = nums(name.substr(name.find_first_of("0123456789")));
        return v.empty() ? Paths64{} : Paths64{circle_path(pt(std::int64_t{0}, 0), v[0] / 2, tol_)};
    }

    std::map<std::pair<std::string, double>, Paths64> symbol_cache_;  // (name, unit) -> shape
    const Paths64& symbol(const FeatureFile& f, int idx) {
        static const Paths64 none;
        if (idx < 0 || static_cast<std::size_t>(idx) >= f.symbols.size()) return none;
        auto key = std::make_pair(lower(f.symbols[static_cast<std::size_t>(idx)]), f.symbol_unit_nm);
        auto it = symbol_cache_.find(key);
        if (it != symbol_cache_.end()) return it->second;
        Paths64 s = symbol_shape(key.first, key.second);
        return symbol_cache_[key] = std::move(s);
    }

    // Features of one record, appended to out; the layer name selects the nets.
    void to_features(const FeatureFile& f, const Rec& r, const std::string& layer, int index, std::vector<Feature>& out) {
        auto P = [&](double x, double y) { return pt(ll(x * f.unit_nm), ll(y * f.unit_nm)); };
        Feature ft;
        ft.dark = r.positive;
        if (!layer.empty()) {
            auto it = feature_net_.find({layer, index});
            if (it != feature_net_.end()) {
                const FeatureNet& fn = it->second;
                ft.has_net = true;
                ft.net = fn.net >= 0 && static_cast<std::size_t>(fn.net) < net_names_.size() ? net_names_[static_cast<std::size_t>(fn.net)] : "";
                if (fn.side) {
                    auto& list = comps_[fn.side];
                    if (fn.comp >= 0 && static_cast<std::size_t>(fn.comp) < list.size()) {
                        const auto& c = list[static_cast<std::size_t>(fn.comp)];
                        ft.pad = true;
                        ft.component = c.refdes;
                        ft.pin = fn.pin >= 0 && static_cast<std::size_t>(fn.pin) < c.pins.size() ? c.pins[static_cast<std::size_t>(fn.pin)] : "";
                    }
                }
            }
        }
        switch (r.type) {
        case 'L': {
            ft.kind = FeatureKind::Stroke;
            ft.at = P(r.v[0], r.v[1]);
            ft.end = P(r.v[2], r.v[3]);
            const std::string sym = r.symbol >= 0 && static_cast<std::size_t>(r.symbol) < f.symbols.size() ? lower(f.symbols[static_cast<std::size_t>(r.symbol)]) : "";
            if (!sym.empty() && sym[0] == 'r') {
                ft.stroke_width_nm = std::stod(sym.substr(1)) * f.symbol_unit_nm;
                ft.shape = {capsule(ft.at, ft.end, ft.stroke_width_nm / 2, tol_)};
            } else {
                const Paths64& s = symbol(f, r.symbol);
                if (s.empty()) return;
                ft.shape = {stroke_convex(ft.at, ft.end, s.front())};
            }
            break;
        }
        case 'A': {
            ft.kind = FeatureKind::Stroke;
            ft.at = P(r.v[0], r.v[1]);
            ft.end = P(r.v[2], r.v[3]);
            ft.arc_centre = P(r.v[4], r.v[5]);
            ft.arc_ccw = !r.cw;
            const std::string sym = r.symbol >= 0 && static_cast<std::size_t>(r.symbol) < f.symbols.size() ? lower(f.symbols[static_cast<std::size_t>(r.symbol)]) : "";
            const double d = !sym.empty() && sym[0] == 'r' ? std::stod(sym.substr(1)) * f.symbol_unit_nm : 0.0;
            if (!(d > 0)) warn("arc drawn with " + sym + ", taken as a round line");
            ft.stroke_width_nm = d;
            const Path64 c = arc_points(ft.at, ft.end, *ft.arc_centre, ft.arc_ccw, tol_);
            ft.shape = d > 0 ? stroke_polyline(c, d / 2, tol_) : Paths64{};
            break;
        }
        case 'P': {
            ft.kind = FeatureKind::Flash;
            ft.at = P(r.v[0], r.v[1]);
            Placement t;
            t.offset = ft.at;
            const int o = r.orient;
            if (o >= 0 && o <= 7) {
                t.rotation_deg = -90.0 * (o % 4);
                t.mirror_x = o >= 4;
            } else {
                t.rotation_deg = -r.rotation;  // clockwise in ODB++
                t.mirror_x = o == 9;
            }
            ft.shape = transform(symbol(f, r.symbol), t);
            break;
        }
        case 'S': {
            ft.kind = FeatureKind::Region;
            Paths64 islands, holes;
            for (const auto& [isl, p] : r.contours) (isl ? islands : holes).insert((isl ? islands : holes).end(), p.begin(), p.end());
            ft.shape = to_paths(holes.empty() ? union_polygons(islands) : difference(union_polygons(islands), union_polygons(holes)));
            if (!ft.shape.empty() && !ft.shape.front().empty()) ft.at = ft.shape.front().front();
            break;
        }
        default:
            return;
        }
        out.push_back(std::move(ft));
    }

    std::vector<Feature> layer_features(const Layer& l) {
        std::vector<Feature> out;
        const auto t = files_->get(layer_path(l.name, "features"));
        if (!t) return out;
        const FeatureFile f = parse_features(*t);
        const std::string key = lower(l.name);
        for (std::size_t i = 0; i < f.recs.size(); ++i) to_features(f, f.recs[i], key, static_cast<int>(i), out);
        return out;
    }

    // ---------------------------------------------------------------- outline, layers

    void read_outline() {
        const auto t = files_->get("steps/" + step_ + "/profile");
        if (!t) return;
        const FeatureFile f = parse_features(*t);
        std::vector<Feature> feats;
        for (std::size_t i = 0; i < f.recs.size(); ++i) to_features(f, f.recs[i], "", static_cast<int>(i), feats);
        board_.outline = render_all(feats);
    }

    void read_layers() {
        const int n = static_cast<int>(board_.copper.size());
        std::vector<std::vector<Feature>> copper(n);
        std::map<Side, std::vector<Feature>> masks;
        int first_copper_row = INT32_MAX;
        for (const auto& l : layers_)
            if (is_copper(l)) first_copper_row = std::min(first_copper_row, l.row);
        for (const auto& l : layers_) {
            if (l.context != "BOARD") continue;
            if (is_copper(l)) {
                if (l.negative) warn("negative copper layer " + l.name + " is not supported");
                copper[copper_index_.at(lower(l.name)) - 1] = layer_features(l);
            } else if (l.type == "SOLDER_MASK") {
                auto& dst = masks[l.row < first_copper_row ? Side::Top : Side::Bottom];
                auto f = layer_features(l);
                dst.insert(dst.end(), std::make_move_iterator(f.begin()), std::make_move_iterator(f.end()));
            } else if (l.type == "DRILL") {
                read_drills(l);
            } else if (l.type == "ROUT") {
                warn("rout layer " + l.name + " is not read");
            }
        }
        std::vector<std::future<std::map<std::string, PolygonSet>>> rendered;
        for (int i = 0; i < n; ++i)
            rendered.push_back(std::async(std::launch::async, [&copper, i]() { return render_by_net(copper[i]); }));
        for (int i = 0; i < n; ++i) detail::set_layer_copper(board_, board_.copper[i], rendered[i].get());
        {
            std::vector<const std::vector<Feature>*> layers;
            for (int i = 0; i < n; ++i) layers.push_back(&copper[i]);
            detail::collect_tracks(board_, layers);
        }
        for (int i = 0; i < n; ++i)
            for (const auto& f : copper[i]) {
                if (!f.is_pad() || !f.dark || f.shape.empty()) continue;
                Pad p;
                p.component = f.component;
                p.pin = f.pin;
                p.function = "ComponentPad";
                p.net = board_.net_id(f.net);
                p.layer = i + 1;
                p.at = f.at;
                p.shape = union_polygons(f.shape);
                board_.pads.push_back(std::move(p));
            }
        for (auto& [side, features] : masks) board_.masks.push_back({side, render_all(features)});
        if (board_.outline.empty()) {
            warn("no profile; the outline is the bounding box of the copper");
            Rect64 r(INT64_MAX, INT64_MAX, INT64_MIN, INT64_MIN);
            for (const auto& l : board_.copper)
                for (const auto& nc : l.nets) {
                    const Rect64 b = bounds(nc.polygons);
                    r.left = std::min(r.left, b.left);
                    r.top = std::min(r.top, b.top);
                    r.right = std::max(r.right, b.right);
                    r.bottom = std::max(r.bottom, b.bottom);
                }
            board_.outline = {Polygon{{pt(r.left, r.top), pt(r.right, r.top), pt(r.right, r.bottom), pt(r.left, r.bottom)}, {}}};
        }
    }

    void read_drills(const Layer& l) {
        const int n = static_cast<int>(board_.copper.size());
        int from = 1, to = n;
        auto a = copper_index_.find(lower(l.start)), b = copper_index_.find(lower(l.end));
        if (a != copper_index_.end() && b != copper_index_.end()) {
            from = std::min(a->second, b->second);
            to = std::max(a->second, b->second);
        }
        const auto t = files_->get(layer_path(l.name, "features"));
        if (!t) return;
        const FeatureFile f = parse_features(*t);
        int drill_attr = -1;
        for (const auto& [i, name] : f.attr_names)
            if (lower(name) == ".drill") drill_attr = i;
        // The layer's tools: type (0 plated, 1 non-plated, 2 via, as .drill) and finished size.
        std::vector<std::pair<int, double>> tools;
        if (const auto tt = files_->get(layer_path(l.name, "tools"))) {
            std::istringstream in(*tt);
            std::string line;
            double unit = 1000.0;  // micrometres (MM) or mils (INCH), in nm
            int type = -1;
            while (std::getline(in, line)) {
                const auto eq = line.find('=');
                if (eq == std::string::npos) continue;
                std::string key = line.substr(0, eq), val = line.substr(eq + 1);
                key.erase(0, key.find_first_not_of(" \t"));
                val.erase(val.find_last_not_of(" \t\r") + 1);
                if (key == "UNITS") unit = lower(val) == "inch" ? 25400.0 : 1000.0;
                else if (key == "TYPE") type = val == "PLATED" ? 0 : val == "NON_PLATED" ? 1 : val == "VIA" ? 2 : -1;
                else if (key == "FINISH_SIZE" && type >= 0) {
                    try {
                        tools.push_back({type, std::stod(val) * unit});
                    } catch (...) {
                    }
                    type = -1;
                }
            }
        }
        const bool one_type = !tools.empty() && std::all_of(tools.begin(), tools.end(), [&](const auto& t) { return t.first == tools.front().first; });
        const bool layer_nonplated = lower(l.name).find("non") != std::string::npos;
        bool guessed = false;
        for (const auto& r : f.recs) {
            Drill d;
            d.from_layer = from;
            d.to_layer = to;
            int kind = -1;
            if (auto it = r.attrs.find(drill_attr); drill_attr >= 0 && it != r.attrs.end()) {
                try {
                    kind = std::stoi(it->second);
                } catch (...) {
                }
            }
            const std::string sym = r.symbol >= 0 && static_cast<std::size_t>(r.symbol) < f.symbols.size() ? lower(f.symbols[static_cast<std::size_t>(r.symbol)]) : "";
            auto P = [&](double x, double y) { return pt(ll(x * f.unit_nm), ll(y * f.unit_nm)); };
            if (r.type == 'P' && !sym.empty() && sym[0] == 'r') {
                d.at = P(r.v[0], r.v[1]);
                d.diameter_nm = std::stod(sym.substr(1)) * f.symbol_unit_nm;
            } else if (r.type == 'P' && sym.rfind("oval", 0) == 0) {
                std::stringstream ss(sym.substr(4));
                std::string a1, b1;
                std::getline(ss, a1, 'x');
                std::getline(ss, b1, 'x');
                const double w = std::stod(a1) * f.symbol_unit_nm, h = std::stod(b1) * f.symbol_unit_nm;
                const double dia = std::min(w, h), half = (std::max(w, h) - dia) / 2;
                Placement t2;
                t2.offset = P(r.v[0], r.v[1]);
                t2.rotation_deg = (r.orient == 8 || r.orient == 9) ? -r.rotation : -90.0 * (r.orient % 4);
                const Point64 p0 = w >= h ? pt(ll(-half), std::int64_t{0}) : pt(std::int64_t{0}, ll(-half));
                const Point64 p1 = w >= h ? pt(ll(half), std::int64_t{0}) : pt(std::int64_t{0}, ll(half));
                d.at = transform(p0, t2);
                d.slot_end = transform(p1, t2);
                d.diameter_nm = dia;
            } else if (r.type == 'L' && !sym.empty() && sym[0] == 'r') {
                d.at = P(r.v[0], r.v[1]);
                d.slot_end = P(r.v[2], r.v[3]);
                d.diameter_nm = std::stod(sym.substr(1)) * f.symbol_unit_nm;
            } else {
                if (r.type != 'X') warn("drill feature with symbol " + sym + " is not read");
                continue;
            }
            // Without .drill: the tool of the finished size, the layer's only tool type, or
            // the layer's name.
            if (kind < 0)
                for (const auto& [type, size] : tools)
                    if (std::abs(size - d.diameter_nm) <= 1000.0) {
                        kind = type;
                        break;
                    }
            if (kind < 0 && one_type) kind = tools.front().first;
            if (kind < 0) {
                kind = layer_nonplated ? 1 : 0;
                guessed = true;
            }
            d.plated = kind != 1;
            d.function = kind == 2 ? DrillFunction::Via : (kind == 0 ? DrillFunction::Component : DrillFunction::Unknown);
            board_.drills.push_back(d);
        }
        if (guessed)
            warn("drill layer " + l.name + ": plating of some holes taken from the layer's name (no .drill attribute or matching tool)");
    }
};

}  // namespace

Board import_odbpp(const std::string& path, const ImportOptions& options) { return Importer(options).run(path); }

bool looks_like_odbpp(const std::string& path) {
    auto matrix = [](std::string p) {
        p = lower(p);
        for (const std::string m : {"matrix/matrix", "matrix/matrix.z"})
            if (p.size() >= m.size() && p.compare(p.size() - m.size(), m.size(), m) == 0) return true;
        return false;
    };
    try {
        if (fs::is_directory(path)) {
            for (const auto& e : fs::recursive_directory_iterator(path))
                if (matrix(e.path().generic_string())) return true;
            return false;
        }
        if (is_zip_file(path)) {
            for (const auto& e : read_zip(path))
                if (matrix(e.name)) return true;
            return false;
        }
        const std::string ext = lower(fs::path(path).extension().string());
        const std::string name = lower(fs::path(path).filename().string());
        if (ext == ".tgz" || ext == ".tar" || (name.size() > 7 && name.compare(name.size() - 7, 7, ".tar.gz") == 0)) {
            const std::string raw = read_file(path);
            const std::string tar = raw.size() > 2 && static_cast<unsigned char>(raw[0]) == 0x1f ? gunzip(raw) : raw;
            for (const auto& [n, d] : read_tar(tar))
                if (matrix(n)) return true;
        }
        return false;
    } catch (const std::exception&) {
        return false;
    }
}

}  // namespace kairoq
