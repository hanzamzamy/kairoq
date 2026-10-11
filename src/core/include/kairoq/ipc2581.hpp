// IPC-2581 (revisions B and C) importer.
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/fabrication.hpp"

#include <string>

namespace kairoq {

// Import an IPC-2581 file (XML, or a zip archive holding one).
Board import_ipc2581(const std::string& path, const ImportOptions& options = {});

// True if the file is IPC-2581 XML or a zip archive whose first XML entry is.
bool looks_like_ipc2581(const std::string& path);

}  // namespace kairoq
