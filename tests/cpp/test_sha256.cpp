#include "sha256.hpp"

#include <doctest/doctest.h>

#include <string>

using namespace kairoq::detail;

// The examples of FIPS 180-4 (NIST), and messages of lengths about the 56- and 64-byte
// boundaries of the padding.
TEST_CASE("SHA-256 test vectors") {
    CHECK(sha256_hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    CHECK(sha256_hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    CHECK(sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
          "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    CHECK(sha256_hex(std::string(1000000, 'a')) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
    const std::pair<int, const char*> lengths[] = {
        {55, "d5e285683cd4efc02d021a5c62014694958901005d6f71e89e0989fac77e4072"},
        {56, "04c26261370ee7541549d16dee320c723e3fd14671e66a099afe0a377c16888e"},
        {63, "75220b47218278e656f2013bb8f0c455a25eaf01e86c64924e9d48d89776d6f2"},
        {64, "7ce100971f64e7001e8fe5a51973ecdfe1ced42befe7ee8d5fd6219506b5393c"},
        {65, "9537c5fdf120482f7d58d25e9ed583f52c02b4e304ea814db1633ad565aed7e9"},
        {119, "000b48d4edf0fa7bee3c6236ecd2785baa5db4eeb8bb54341b029e0d9fa5fb0c"},
        {120, "13f05a0b594787f5ecd315edc96141bd3243203d1b7d4f0836f37308b276ba98"}};
    for (const auto& [n, hex] : lengths) CHECK(sha256_hex(std::string(static_cast<std::size_t>(n), 'x')) == hex);
}

TEST_CASE("SHA-256 in pieces") {
    const std::string m(1000, 'q');
    for (std::size_t piece : {1u, 7u, 63u, 64u, 65u, 333u}) {
        Sha256 s;
        for (std::size_t i = 0; i < m.size(); i += piece) s.update(std::string_view(m).substr(i, piece));
        CHECK(s.hex() == sha256_hex(m));
    }
}
