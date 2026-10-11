// File and zip archive helpers.
#pragma once

#include <string>
#include <utility>
#include <vector>

namespace kairoq {

struct ArchiveEntry {
    std::string name;
    std::string data;
};

std::string read_file(const std::string& path);
void write_file(const std::string& path, const std::string& data);

bool is_zip_file(const std::string& path);
std::vector<ArchiveEntry> read_zip(const std::string& path);
void write_zip(const std::string& path, const std::vector<ArchiveEntry>& entries);

// Content of a gzip stream (RFC 1952, the first member).
std::string gunzip(const std::string& data);

}  // namespace kairoq
