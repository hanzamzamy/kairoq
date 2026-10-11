// ODB++ importer (job directory, or a .zip or .tgz archive of one).
#pragma once

#include "kairoq/board.hpp"
#include "kairoq/fabrication.hpp"

#include <string>
#include <utility>
#include <vector>

namespace kairoq {

// Import an ODB++ job: the first step's copper, solder mask, drills, outline, nets,
// components and stackup.
Board import_odbpp(const std::string& path, const ImportOptions& options = {});

// True if the path is a directory, .zip or .tgz holding an ODB++ job (matrix/matrix).
bool looks_like_odbpp(const std::string& path);

// Archive formats ODB++ uses besides zip and gzip: tar (POSIX ustar and GNU long names)
// and Unix compress (.Z, LZW).
std::vector<std::pair<std::string, std::string>> read_tar(const std::string& data);
std::string uncompress_lzw(const std::string& data);

}  // namespace kairoq
