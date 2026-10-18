#include "kairoq/open.hpp"

#include "kairoq/archive.hpp"
#include "kairoq/fabrication.hpp"
#include "kairoq/ipc2581.hpp"
#include "kairoq/model_io.hpp"
#include "kairoq/odbpp.hpp"
#include "kairoq/provenance.hpp"

#include <cstdio>
#include <filesystem>
#include <map>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

// A stackup file replaces the stackup of a board read from IPC-2581, ODB++ or a model; the
// fabrication importer applies it itself, in place of the job file's.
void override_stackup(Board& b, const std::string& path) {
    std::vector<StackupLayer> layers = parse_stackup_yaml(read_file(path));
    int copper = 0;
    for (const auto& l : layers)
        if (l.type == LayerType::Copper) ++copper;
    if (copper != int(b.copper.size()))
        throw std::runtime_error("the stackup file has " + std::to_string(copper) + " copper layers, the board " +
                                 std::to_string(b.copper.size()));
    b.stackup = std::move(layers);
    b.stackup_source = "yaml";
    b.source.files.push_back(std::filesystem::path(path).filename().string());
}

Board read_board(const std::vector<std::string>& inputs, const ImportOptions& options) {
    if (inputs.empty()) throw std::runtime_error("no input given");
    if (inputs.size() == 1) {
        Board b;
        if (is_board_model(inputs[0])) b = load_board(inputs[0]);
        else if (looks_like_odbpp(inputs[0])) b = import_odbpp(inputs[0], options);
        else if (looks_like_ipc2581(inputs[0])) b = import_ipc2581(inputs[0], options);
        else return import_fabrication(inputs, options);
        if (!options.stackup_yaml_path.empty()) override_stackup(b, options.stackup_yaml_path);
        return b;
    }
    return import_fabrication(inputs, options);
}

}  // namespace

Board open_board(const std::vector<std::string>& inputs, const ImportOptions& options) {
    Board b = read_board(inputs, options);
    // A board read from a canonical model keeps the digests of the files it was imported from.
    if (!(inputs.size() == 1 && is_board_model(inputs[0]))) {
        b.source.digests.clear();
        // By file name, as the board names the files it read: a model does not record where
        // they lay on the machine that imported them.
        auto add = [&](const std::string& path) {
            FileDigest d = digest(path);
            d.path = std::filesystem::path(path).lexically_normal().filename().string();
            if (d.path.empty()) d.path = std::filesystem::path(path).lexically_normal().parent_path().filename().string();
            b.source.digests.push_back(d);
        };
        for (const auto& path : inputs) add(path);
        if (!options.stackup_yaml_path.empty()) add(options.stackup_yaml_path);
    }
    return options.crop ? crop_board(b, *options.crop) : b;
}

namespace {

std::string fmt(const char* f, double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, f, v);
    return buf;
}

}  // namespace

std::string board_summary(const Board& b) {
    std::ostringstream s;
    s << "board        " << b.name << "\n";
    s << "source       " << b.source.importer;
    if (!b.source.generator.empty()) s << " (" << b.source.generator << ")";
    s << "\n";
    const Rect64 r = bounds(b.outline);
    s << "outline      " << fmt("%.3f", to_mm(double(r.right - r.left))) << " x " << fmt("%.3f", to_mm(double(r.bottom - r.top)))
      << " mm, area " << fmt("%.2f", area_nm2(b.outline) / 1e12) << " mm^2\n";
    s << "stackup      " << b.stackup_source << ", " << fmt("%.4f", to_mm(b.thickness_nm())) << " mm\n";
    for (const auto& l : b.stackup) {
        char line[160];
        std::snprintf(line, sizeof line, "  %-11s %-24s %9.4f mm  %s", to_string(l.type), l.name.c_str(), to_mm(l.thickness_nm),
                      l.material.c_str());
        s << line << "\n";
    }
    s << "copper       " << b.copper.size() << " layers, nets from " << b.net_source << "\n";
    for (const auto& l : b.copper) {
        double area = 0;
        for (const auto& nc : l.nets) area += area_nm2(nc.polygons);
        char line[160];
        std::snprintf(line, sizeof line, "  L%-3d %-8s %-24s %5zu nets %12.3f mm^2", l.index, to_string(l.side), l.name.c_str(),
                      l.nets.size(), area / 1e12);
        s << line << "\n";
    }
    s << "nets         " << (b.nets.empty() ? 0 : b.nets.size() - 1) << "\n";
    s << "pads         " << b.pads.size() << "\n";
    std::map<std::string, int> drills;
    for (const auto& d : b.drills) drills[std::string(d.plated ? "plated " : "non-plated ") + to_string(d.function)]++;
    s << "drills       " << b.drills.size() << "\n";
    for (const auto& [k, n] : drills) s << "  " << k << ": " << n << "\n";
    s << "components   " << b.components.size() << "\n";
    s << "masks        " << b.masks.size() << "\n";
    if (!b.warnings.empty()) {
        s << "warnings     " << b.warnings.size() << "\n";
        for (const auto& w : b.warnings) s << "  " << w << "\n";
    }
    return s.str();
}

}  // namespace kairoq
