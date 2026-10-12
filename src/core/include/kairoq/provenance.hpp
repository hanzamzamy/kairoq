// Provenance of a result: the build that computed it and the digests of everything it was
// computed from, so that a result can be traced to its exact inputs (docs/provenance.md).
#pragma once

#include "kairoq/board.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace kairoq {

struct BuildInfo {
    std::string version;
    std::string commit;  // git commit of the source; "+modified" when the tree had changes; "unknown" outside git
    std::string compiler;
    std::string platform;
};
BuildInfo build_info();

// The SHA-256 of a file, or of a directory (see FileDigest). Throws if it cannot be read.
FileDigest digest(const std::string& path);

struct Provenance {
    BuildInfo build;
    std::string time_utc;               // ISO 8601
    std::vector<std::string> command;   // the command line, when there is one
    std::vector<FileDigest> inputs;     // the files given, and those a project file names
    std::vector<FileDigest> board_sources;  // the files the board was imported from
};

// The provenance of a result computed now from the given input paths. A project file
// (YAML) among them also contributes every file its values name, relative to its
// directory (boards, stackups, Touchstone and IBIS files), except a value under an `output`
// key. Paths that do not exist are left out.
Provenance make_provenance(const std::vector<std::string>& command, const std::vector<std::string>& inputs, const std::string& project,
                           const Board* board);
nlohmann::json provenance_json(const Provenance& p);
// The provenance as comment lines, each starting with the prefix ("! " for Touchstone).
std::string provenance_comment(const Provenance& p, const std::string& prefix);

}  // namespace kairoq
