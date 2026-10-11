#include "kairoq/open.hpp"

#include "kairoq/ipc2581.hpp"
#include "kairoq/model_io.hpp"
#include "kairoq/odbpp.hpp"

#include <cstdio>
#include <map>
#include <sstream>
#include <stdexcept>

namespace kairoq {

namespace {

Board read_board(const std::vector<std::string>& inputs, const ImportOptions& options) {
    if (inputs.empty()) throw std::runtime_error("no input given");
    if (inputs.size() == 1) {
        if (is_board_model(inputs[0])) return load_board(inputs[0]);
        if (looks_like_odbpp(inputs[0])) return import_odbpp(inputs[0], options);
        if (looks_like_ipc2581(inputs[0])) return import_ipc2581(inputs[0], options);
    }
    return import_fabrication(inputs, options);
}

}  // namespace

Board open_board(const std::vector<std::string>& inputs, const ImportOptions& options) {
    Board b = read_board(inputs, options);
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
