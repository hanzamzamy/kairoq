// KairoQ command-line tool.
#include "kairoq/archive.hpp"
#include "kairoq/dc.hpp"
#include "kairoq/deembed.hpp"
#include "kairoq/model_io.hpp"
#include "kairoq/network.hpp"
#include "kairoq/open.hpp"
#include "kairoq/openems.hpp"
#include "kairoq/pdn.hpp"
#include "kairoq/server.hpp"
#include "kairoq/si.hpp"
#include "kairoq/emi.hpp"
#include "kairoq/elmer.hpp"
#include "kairoq/eye.hpp"
#include "kairoq/thermal.hpp"
#include "kairoq/timedomain.hpp"

#include <nlohmann/json.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

const char* kUsage = R"(usage:
  kairoq import <input>... -o <model> [options]   write the canonical model
  kairoq info <input>... [--json] [options]       summarise a board
  kairoq view <input>... [--port N] [--host H] [--no-browser] [options]
  kairoq dc <project.yaml> [<input>...] [--view] [options]   DC current solution
  kairoq thermal <project.yaml> [<input>...] [--view] [options]   steady or transient temperature
  kairoq emi <project.yaml> [<input>...] [--view] [options]   loops, return paths, switch-node coupling
  kairoq openems <project.yaml> [<input>...] -o <dir> [--loop NAME | --channel NAME] [options]   export a loop or a channel to openEMS
  kairoq si <project.yaml> [<input>...] [--view] [options]   impedance, delay and loss of nets; pairs
  kairoq sparams <project.yaml> [<input>...] [-o <dir>] [options]   S-parameters and TDR of channels
  kairoq waveforms <project.yaml> [<input>...] [-o <dir>] [options]   time-domain response with drivers
  kairoq elmer <project.yaml> [<input>...] -o <dir> [--run] [--thermal] [options]   the DC or thermal problem for Elmer (3D reference)
  kairoq deembed <spec.yaml>   de-embed a Touchstone file: reference, lines, 2x-thru
  kairoq eye <project.yaml> [<input>...] [-o <dir>] [options]   eye diagrams of channels at a bit error ratio
  kairoq pdn <project.yaml> [<input>...] [-o <dir>] [options]   impedance of a power distribution network

An input is a canonical model (directory or .zip), an IPC-2581 file (.xml or .zip),
an ODB++ job (directory, .zip or .tgz), or fabrication files: Gerber X2, Excellon or
Gerber drill, .gbrjob, IPC-D-356, placement CSV, given as files, directories or a zip
archive.

options:
  --stackup <file.yaml>     stackup, overriding the job file
  --arc-tolerance <um>      largest chord deviation from arcs (default 1 um)
  --crop <x0,y0,x1,y1>      the part of the board inside this rectangle (mm); its edges are open
)";

struct Args {
    std::string command;
    std::vector<std::string> inputs;
    std::string output;
    std::string loop;
    std::string channel;
    bool json = false;
    bool view = false;
    bool run = false;
    bool thermal = false;
    kairoq::ImportOptions import;
    kairoq::ServeOptions serve;
};

Args parse(int argc, char** argv) {
    if (argc < 2) throw std::invalid_argument("missing command");
    Args a;
    a.command = argv[1];
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
            return argv[++i];
        };
        if (arg == "-o" || arg == "--output") a.output = value();
        else if (arg == "--json") a.json = true;
        else if (arg == "--loop") a.loop = value();
        else if (arg == "--run") a.run = true;
        else if (arg == "--channel") a.channel = value();
        else if (arg == "--stackup") a.import.stackup_yaml_path = value();
        else if (arg == "--arc-tolerance") a.import.geometry.arc_tolerance_nm = std::stod(value()) * 1000.0;
        else if (arg == "--crop") {
            std::stringstream text(value());
            std::array<double, 4> c{};
            char comma = ',';
            text >> c[0] >> comma >> c[1] >> comma >> c[2] >> comma >> c[3];
            if (!text) throw std::invalid_argument("--crop needs x0,y0,x1,y1 in mm");
            a.import.crop = Clipper2Lib::Rect64(kairoq::nm_from_mm(c[0]), kairoq::nm_from_mm(c[1]), kairoq::nm_from_mm(c[2]), kairoq::nm_from_mm(c[3]));
        }
        else if (arg == "--thermal") a.thermal = true;
        else if (arg == "--port") a.serve.port = std::stoi(value());
        else if (arg == "--host") a.serve.host = value();
        else if (arg == "--no-browser") a.serve.open_browser = false;
        else if (arg == "--view") a.view = true;
        else if (arg.rfind("--", 0) == 0) throw std::invalid_argument("unknown option " + arg);
        else a.inputs.push_back(arg);
    }
    if (a.inputs.empty()) throw std::invalid_argument("no input given");
    return a;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc >= 2 && (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help"))) {
        std::cout << kUsage;
        return 0;
    }
    Args args;
    try {
        args = parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "kairoq: " << e.what() << "\n\n" << kUsage;
        return 2;
    }
    try {
        if (args.command == "dc") {
            const kairoq::Project project = kairoq::load_project(args.inputs[0]);
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const kairoq::DcResult result = kairoq::solve_dc(board, project.dc);
            std::cout << kairoq::dc_summary(board, result);
            if (args.view) kairoq::serve(board, args.serve, &result);
            return 0;
        }
        if (args.command == "elmer") {
            const kairoq::Project project = kairoq::load_project(args.inputs[0]);
            if (args.output.empty()) throw std::invalid_argument("elmer needs -o <dir>");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            if (args.thermal) {
                const kairoq::ThermalProject tp = kairoq::load_thermal_project(args.inputs[0]);
                if (!tp.thermal) throw std::invalid_argument("the project has no thermal section");
                std::vector<std::string> warnings;
                const kairoq::ThermalSettings ts = kairoq::resolve_thermal_project(board, tp, warnings);
                const kairoq::ElmerSettings es = kairoq::parse_elmer_settings(kairoq::read_file(args.inputs[0]), "thermal");
                const auto e = kairoq::export_elmer_thermal(board, ts, es, args.output);
                if (!args.run) {
                    std::cout << kairoq::elmer_thermal_summary(e, nullptr, nullptr) << "wrote " << args.output << "; run ElmerSolver on thermal.sif there\n";
                    return 0;
                }
                const auto r = kairoq::run_elmer_thermal(std::filesystem::absolute(args.output).string(), es.solver);
                const kairoq::ThermalResult native = kairoq::solve_thermal(board, ts);
                std::cout << kairoq::elmer_thermal_summary(e, &r, &native);
                return 0;
            }
            const kairoq::ElmerSettings es = kairoq::parse_elmer_settings(kairoq::read_file(args.inputs[0]));
            const auto e = kairoq::export_elmer_dc(board, project.dc, es, args.output);
            if (!args.run) {
                std::cout << kairoq::elmer_dc_summary(e, nullptr, nullptr) << "wrote " << args.output << "; run ElmerSolver on each port_<k>.sif there\n";
                return 0;
            }
            const auto r = kairoq::run_elmer_dc(project.dc, std::filesystem::absolute(args.output).string(), es.solver);
            const kairoq::DcResult native = kairoq::solve_dc(board, project.dc);
            std::cout << kairoq::elmer_dc_summary(e, &r, &native);
            return 0;
        }
        if (args.command == "emi") {
            const kairoq::EmiProject project = kairoq::load_emi_project(args.inputs[0]);
            if (!project.emi) throw std::invalid_argument("the project has no emi section");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const kairoq::EmiResult result = kairoq::analyse_emi(board, *project.emi);
            std::cout << kairoq::emi_summary(board, result);
            if (args.view) kairoq::serve(board, args.serve, nullptr, nullptr, &result);
            return 0;
        }
        if (args.command == "si") {
            const kairoq::SiProject project = kairoq::load_si_project(args.inputs[0]);
            if (!project.si) throw std::invalid_argument("the project has no si section");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const kairoq::SiResult result = kairoq::analyse_si(board, *project.si);
            std::cout << kairoq::si_summary(board, result);
            if (args.view) kairoq::serve(board, args.serve, nullptr, nullptr, nullptr, &result);
            return 0;
        }
        if (args.command == "waveforms") {
            const kairoq::SiProject project = kairoq::load_si_project(args.inputs[0]);
            if (!project.si) throw std::invalid_argument("the project has no si section");
            std::ifstream pf(args.inputs[0]);
            std::stringstream text;
            text << pf.rdbuf();
            const kairoq::ChannelProject cp = kairoq::parse_channel_project(text.str());
            const auto sims = kairoq::parse_simulations(text.str());
            if (sims.empty()) throw std::invalid_argument("the project has no si.simulations");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const std::string base = std::filesystem::absolute(args.inputs[0]).parent_path().string();
            std::map<std::string, kairoq::ChannelResult> channels;
            std::vector<kairoq::SimulationResult> results;
            for (const auto& sp : sims) {
                if (!channels.count(sp.channel)) {
                    const kairoq::ChannelSpec* c = nullptr;
                    for (const auto& x : cp.channels)
                        if (x.name == sp.channel) c = &x;
                    if (!c) throw std::invalid_argument("simulation " + sp.name + ": no channel " + sp.channel);
                    channels[sp.channel] = kairoq::analyse_channels(board, *project.si, {*c}, cp.network, base).front();
                }
                results.push_back(kairoq::simulate(channels[sp.channel], sp, base));
            }
            std::cout << kairoq::simulations_summary(results);
            if (!args.output.empty()) {
                std::filesystem::create_directories(args.output);
                for (const auto& r : results) {
                    std::ofstream f(std::filesystem::path(args.output) / (r.name + ".csv"));
                    f << "time_s";
                    for (const auto& p : r.ports) f << "," << p;
                    f << "\n";
                    for (std::size_t k = 0; k < r.time_s.size(); ++k) {
                        f << r.time_s[k];
                        for (const auto& v : r.voltage) f << "," << v[k];
                        f << "\n";
                    }
                }
                std::cout << "wrote " << args.output << "\n";
            }
            return 0;
        }
        if (args.command == "deembed") {
            const kairoq::DeembedResult r = kairoq::deembed(kairoq::read_file(args.inputs[0]),
                                                                std::filesystem::absolute(args.inputs[0]).parent_path().string());
            for (const auto& s : r.steps) std::cout << s << "\n";
            if (!r.output.empty()) std::cout << "wrote " << r.output << "\n";
            return 0;
        }
        if (args.command == "eye") {
            const kairoq::SiProject project = kairoq::load_si_project(args.inputs[0]);
            if (!project.si) throw std::invalid_argument("the project has no si section");
            std::ifstream pf(args.inputs[0]);
            std::stringstream text;
            text << pf.rdbuf();
            const kairoq::ChannelProject cp = kairoq::parse_channel_project(text.str());
            const auto eyes = kairoq::parse_eyes(text.str());
            if (eyes.empty()) throw std::invalid_argument("the project has no si.eyes");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const std::string base = std::filesystem::absolute(args.inputs[0]).parent_path().string();
            std::map<std::string, kairoq::ChannelResult> channels;
            std::vector<kairoq::EyeResult> results;
            for (const auto& e : eyes) {
                if (!channels.count(e.channel)) {
                    const kairoq::ChannelSpec* c = nullptr;
                    for (const auto& x : cp.channels)
                        if (x.name == e.channel) c = &x;
                    if (!c) throw std::invalid_argument("eye " + e.name + ": no channel " + e.channel);
                    channels[e.channel] = kairoq::analyse_channels(board, *project.si, {*c}, cp.network, base).front();
                }
                results.push_back(kairoq::eye(channels[e.channel], e));
            }
            std::cout << kairoq::eyes_summary(results);
            if (!args.output.empty()) {
                std::filesystem::create_directories(args.output);
                for (const auto& r : results) {
                    const auto dir = std::filesystem::path(args.output);
                    std::ofstream(dir / (r.name + ".svg")) << kairoq::eye_svg(r);
                    std::ofstream f(dir / (r.name + "_pulse.csv"));
                    f << "time_s,voltage_v\n";
                    for (std::size_t k = 0; k < r.pulse_t.size(); ++k) f << r.pulse_t[k] << "," << r.pulse_v[k] << "\n";
                    std::ofstream g(dir / (r.name + "_contour.csv"));
                    g << "phase_s,upper_v,lower_v\n";
                    for (std::size_t k = 0; k < r.phase_t.size(); ++k) g << r.phase_t[k] << "," << r.upper_v[k] << "," << r.lower_v[k] << "\n";
                }
                std::cout << "wrote " << args.output << "\n";
            }
            return 0;
        }
        if (args.command == "pdn") {
            const kairoq::Project project = kairoq::load_project(args.inputs[0]);
            const kairoq::PdnSettings settings = kairoq::parse_pdn(kairoq::read_file(args.inputs[0]));
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const kairoq::PdnResult r = kairoq::analyse_pdn(board, settings);
            std::cout << kairoq::pdn_summary(r, settings);
            if (!args.output.empty()) {
                std::filesystem::create_directories(args.output);
                for (std::size_t l = 0; l < r.loads.size(); ++l) {
                    std::ofstream f(std::filesystem::path(args.output) / (r.loads[l] + "_z.csv"));
                    f << "frequency_hz,z_re_ohm,z_im_ohm,z_abs_ohm\n";
                    for (std::size_t k = 0; k < r.frequency_hz.size(); ++k)
                        f << r.frequency_hz[k] << "," << r.z_ohm[l][k].real() << "," << r.z_ohm[l][k].imag() << "," << std::abs(r.z_ohm[l][k]) << "\n";
                }
                std::cout << "wrote " << args.output << "\n";
            }
            return 0;
        }
        if (args.command == "sparams") {
            const kairoq::SiProject project = kairoq::load_si_project(args.inputs[0]);
            if (!project.si) throw std::invalid_argument("the project has no si section");
            std::ifstream pf(args.inputs[0]);
            std::stringstream text;
            text << pf.rdbuf();
            const kairoq::ChannelProject cp = kairoq::parse_channel_project(text.str());
            if (cp.channels.empty()) throw std::invalid_argument("the project has no si.channels");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const std::string base = std::filesystem::absolute(args.inputs[0]).parent_path().string();
            const auto results = kairoq::analyse_channels(board, *project.si, cp.channels, cp.network, base);
            std::cout << kairoq::channels_summary(results);
            if (!args.output.empty()) {
                std::filesystem::create_directories(args.output);
                for (const auto& r : results) {
                    const std::string n = std::to_string(r.ports.size());
                    std::ofstream(std::filesystem::path(args.output) / (r.name + ".s" + n + "p")) << kairoq::touchstone(r, false);
                    if (!r.mixed.empty()) {
                        const auto m = static_cast<std::size_t>(std::llround(std::sqrt(double(r.mixed.front().size()))));
                        std::ofstream(std::filesystem::path(args.output) / (r.name + "_mixed.s" + std::to_string(m) + "p")) << kairoq::touchstone(r, true);
                    }
                }
                std::cout << "wrote " << args.output << "\n";
            }
            return 0;
        }
        if (args.command == "openems" && !args.channel.empty()) {
            const kairoq::SiProject project = kairoq::load_si_project(args.inputs[0]);
            if (!project.si) throw std::invalid_argument("the project has no si section");
            if (args.output.empty()) throw std::invalid_argument("openems needs -o <dir>");
            const std::string text = kairoq::read_file(args.inputs[0]);
            const kairoq::ChannelProject cp = kairoq::parse_channel_project(text);
            const kairoq::ChannelSpec* c = nullptr;
            for (const auto& x : cp.channels)
                if (x.name == args.channel) c = &x;
            if (!c) throw std::invalid_argument("the project has no channel " + args.channel);
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const auto model = kairoq::export_openems_channel(board, *project.si, *c, kairoq::parse_openems_channel(text), args.output);
            std::cout << kairoq::openems_channel_summary(model) << "wrote " << args.output << "/model.json and simulate.py\n";
            return 0;
        }
        if (args.command == "openems") {
            const kairoq::EmiProject project = kairoq::load_emi_project(args.inputs[0]);
            if (!project.emi || project.emi->loops.empty()) throw std::invalid_argument("the project has no emi loops");
            if (args.output.empty()) throw std::invalid_argument("openems needs -o <dir>");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::LoopSpec* loop = &project.emi->loops.front();
            if (!args.loop.empty()) {
                loop = nullptr;
                for (const auto& l : project.emi->loops)
                    if (l.name == args.loop) loop = &l;
                if (!loop) throw std::invalid_argument("the project has no loop " + args.loop);
            }
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            const auto model = kairoq::export_openems(board, *loop, *project.emi, args.output);
            std::cout << kairoq::openems_summary(model) << "wrote " << args.output << "/model.json and simulate.py\n";
            return 0;
        }
        if (args.command == "thermal") {
            const kairoq::ThermalProject project = kairoq::load_thermal_project(args.inputs[0]);
            if (!project.thermal) throw std::invalid_argument("the project has no thermal section");
            std::vector<std::string> inputs(args.inputs.begin() + 1, args.inputs.end());
            if (inputs.empty()) inputs = project.board;
            if (inputs.empty()) throw std::invalid_argument("the project names no board and none is given");
            const kairoq::Board board = kairoq::open_board(inputs, args.import);
            std::vector<std::string> warnings;
            const kairoq::ThermalSettings settings = kairoq::resolve_thermal_project(board, project, warnings);
            kairoq::ThermalResult result = kairoq::solve_thermal(board, settings);
            result.warnings.insert(result.warnings.begin(), warnings.begin(), warnings.end());
            std::cout << kairoq::thermal_summary(board, result);
            if (args.view) kairoq::serve(board, args.serve, result.dc ? &*result.dc : nullptr, &result);
            return 0;
        }
        const auto t0 = std::chrono::steady_clock::now();
        const kairoq::Board board = kairoq::open_board(args.inputs, args.import);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (args.command == "import") {
            if (args.output.empty()) throw std::invalid_argument("import needs -o <model>");
            kairoq::save_board(board, args.output);
            std::cerr << "read in " << seconds << " s; wrote " << args.output << " (" << board.warnings.size() << " warnings)\n";
        } else if (args.command == "info") {
            if (args.json) std::cout << kairoq::board_json(board, false).dump(1) << "\n";
            else std::cout << kairoq::board_summary(board) << "read in " << seconds << " s\n";
        } else if (args.command == "view") {
            kairoq::serve(board, args.serve);
        } else {
            std::cerr << "kairoq: unknown command " << args.command << "\n\n" << kUsage;
            return 2;
        }
    } catch (const std::exception& e) {
        std::cerr << "kairoq: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
