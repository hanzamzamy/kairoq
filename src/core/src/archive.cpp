#include "kairoq/archive.hpp"

#include <miniz.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace kairoq {

std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void write_file(const std::string& path, const std::string& data) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);
    out.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!out) throw std::runtime_error("cannot write " + path);
}

bool is_zip_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    char magic[4] = {};
    in.read(magic, 4);
    return in.gcount() == 4 && magic[0] == 'P' && magic[1] == 'K' && magic[2] == 3 && magic[3] == 4;
}

std::vector<ArchiveEntry> read_zip(const std::string& path) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_file(&zip, path.c_str(), 0)) throw std::runtime_error("cannot read zip archive " + path);
    std::vector<ArchiveEntry> out;
    const mz_uint n = mz_zip_reader_get_num_files(&zip);
    for (mz_uint i = 0; i < n; ++i) {
        mz_zip_archive_file_stat st;
        if (!mz_zip_reader_file_stat(&zip, i, &st) || st.m_is_directory) continue;
        std::size_t size = 0;
        void* data = mz_zip_reader_extract_to_heap(&zip, i, &size, 0);
        if (!data) {
            mz_zip_reader_end(&zip);
            throw std::runtime_error(std::string("cannot extract ") + st.m_filename + " from " + path);
        }
        out.push_back({st.m_filename, std::string(static_cast<const char*>(data), size)});
        mz_free(data);
    }
    mz_zip_reader_end(&zip);
    return out;
}

void write_zip(const std::string& path, const std::vector<ArchiveEntry>& entries) {
    mz_zip_archive zip{};
    if (!mz_zip_writer_init_file(&zip, path.c_str(), 0)) throw std::runtime_error("cannot create " + path);
    for (const auto& e : entries) {
        if (!mz_zip_writer_add_mem(&zip, e.name.c_str(), e.data.data(), e.data.size(), MZ_DEFAULT_COMPRESSION)) {
            mz_zip_writer_end(&zip);
            throw std::runtime_error("cannot add " + e.name + " to " + path);
        }
    }
    const bool ok = mz_zip_writer_finalize_archive(&zip);
    mz_zip_writer_end(&zip);
    if (!ok) throw std::runtime_error("cannot finalise " + path);
}

std::string gunzip(const std::string& data) {
    const auto* p = reinterpret_cast<const unsigned char*>(data.data());
    const std::size_t n = data.size();
    if (n < 18 || p[0] != 0x1f || p[1] != 0x8b || p[2] != 8) throw std::runtime_error("not gzip data");
    const unsigned flags = p[3];
    std::size_t at = 10;
    if (flags & 4) {
        if (at + 2 > n) throw std::runtime_error("truncated gzip header");
        at += 2 + (p[at] | (p[at + 1] << 8));
    }
    if (flags & 8)
        while (at < n && p[at++] != 0) {
        }
    if (flags & 16)
        while (at < n && p[at++] != 0) {
        }
    if (flags & 2) at += 2;
    if (at >= n) throw std::runtime_error("truncated gzip data");
    std::size_t out_len = 0;
    void* out = tinfl_decompress_mem_to_heap(p + at, n - at, &out_len, 0);
    if (!out) throw std::runtime_error("corrupt gzip data");
    std::string result(static_cast<const char*>(out), out_len);
    mz_free(out);
    return result;
}

}  // namespace kairoq
