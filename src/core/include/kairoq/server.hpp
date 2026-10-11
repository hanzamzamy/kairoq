// Local HTTP server for the 3D viewer.
#pragma once

#include "kairoq/board.hpp"

#include <functional>
#include <string>

namespace kairoq {

struct DcResult;
struct EmiResult;
struct SiResult;
struct ThermalResult;

struct ServeOptions {
    std::string host = "127.0.0.1";
    int port = 8150;  // 0: any free port
    bool open_browser = true;
    // Polled a few times per second; the server stops when it returns true.
    std::function<bool()> should_stop;
};

// Serve the viewer for the board, with the fields of a DC and a thermal solution, the EMI
// estimates and the signal-integrity results when given, until interrupted or should_stop returns true. Prints the URL.
void serve(const Board& board, const ServeOptions& options = {}, const DcResult* dc = nullptr,
           const ThermalResult* thermal = nullptr, const EmiResult* emi = nullptr, const SiResult* si = nullptr);

}  // namespace kairoq
