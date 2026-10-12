// SHA-256 (FIPS 180-4) of byte strings and files, for the digests of a result's inputs.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace kairoq::detail {

class Sha256 {
public:
    Sha256();
    void update(const void* data, std::size_t size);
    void update(std::string_view s) { update(s.data(), s.size()); }
    std::string hex();  // finishes the digest; lower-case hexadecimal

private:
    void block(const std::uint8_t* p);
    std::array<std::uint32_t, 8> h_;
    std::array<std::uint8_t, 64> buf_{};
    std::size_t used_ = 0;
    std::uint64_t bytes_ = 0;
};

std::string sha256_hex(std::string_view data);
std::string sha256_file(const std::string& path);  // throws if the file cannot be read

}  // namespace kairoq::detail
