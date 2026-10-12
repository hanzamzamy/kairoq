// Provenance of results (provenance.hpp, docs/provenance.md).
#include "kairoq/provenance.hpp"

#include "sha256.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <ctime>
#include <filesystem>
#include <functional>
#include <set>
#include <sstream>

namespace kairoq {

namespace fs = std::filesystem;
using json = nlohmann::json;

FileDigest digest(const std::string& path) {
    FileDigest d;
    d.path = path;
    if (!fs::is_directory(path)) {
        d.sha256 = detail::sha256_file(path);
        d.bytes = fs::file_size(path);
        return d;
    }
    // A directory: its files in sorted order of their relative paths, one line each.
    d.directory = true;
    std::vector<std::pair<std::string, fs::path>> files;
    for (const auto& e : fs::recursive_directory_iterator(path))
        if (e.is_regular_file()) files.push_back({fs::relative(e.path(), path).generic_string(), e.path()});
    std::sort(files.begin(), files.end());
    detail::Sha256 s;
    for (const auto& [rel, full] : files) {
        s.update(rel + " " + detail::sha256_file(full.string()) + "\n");
        d.bytes += fs::file_size(full);
    }
    d.sha256 = s.hex();
    return d;
}

namespace {

std::string utc_now() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

// The files a project's values name, relative to its directory, except those under an
// `output` key.
void named_files(const std::string& project, std::vector<std::pair<std::string, std::string>>& out) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(project);
    } catch (const std::exception&) {
        return;  // not YAML: only the file itself is digested
    }
    const fs::path base = fs::absolute(project).parent_path();
    std::function<void(const YAML::Node&)> walk = [&](const YAML::Node& n) {
        if (n.IsScalar()) {
            const std::string v = n.Scalar();
            if (v.empty() || v.size() > 4096 || v.find('\n') != std::string::npos) return;
            std::error_code ec;
            const fs::path p = fs::path(v).is_absolute() ? fs::path(v) : base / v;
            if (fs::is_regular_file(p, ec) || fs::is_directory(p, ec)) out.push_back({v, p.string()});
        } else if (n.IsMap()) {
            for (const auto& c : n)
                if (c.first.Scalar() != "output") walk(c.second);  // a file the run writes, not an input
        } else if (n.IsSequence()) {
            for (const auto& c : n) walk(c);
        }
    };
    walk(root);
}

json digest_json(const FileDigest& d) {
    json j = {{"path", d.path}, {"sha256", d.sha256}, {"bytes", d.bytes}};
    if (d.directory) j["directory"] = true;
    return j;
}

}  // namespace

Provenance make_provenance(const std::vector<std::string>& command, const std::vector<std::string>& inputs, const std::string& project,
                           const Board* board) {
    Provenance p;
    p.build = build_info();
    p.time_utc = utc_now();
    p.command = command;
    // (name as written, path to read): the project's files are named as the project names them.
    std::vector<std::pair<std::string, std::string>> paths;
    if (!project.empty()) {
        paths.push_back({project, project});
        named_files(project, paths);
    }
    for (const auto& i : inputs) paths.push_back({i, i});
    std::set<std::string> seen;
    for (const auto& [name, path] : paths) {
        std::error_code ec;
        if (!fs::exists(path, ec)) continue;
        const std::string key = fs::weakly_canonical(path, ec).string();
        if (!seen.insert(key).second) continue;
        FileDigest d = digest(path);
        d.path = name;
        p.inputs.push_back(d);
    }
    if (board) p.board_sources = board->source.digests;
    return p;
}

json provenance_json(const Provenance& p) {
    json inputs = json::array(), sources = json::array();
    for (const auto& d : p.inputs) inputs.push_back(digest_json(d));
    for (const auto& d : p.board_sources) sources.push_back(digest_json(d));
    return {{"kairoq", {{"version", p.build.version}, {"commit", p.build.commit}, {"compiler", p.build.compiler}, {"platform", p.build.platform}}},
            {"time_utc", p.time_utc},
            {"command", p.command},
            {"inputs", inputs},
            {"board_sources", sources}};
}

std::string provenance_comment(const Provenance& p, const std::string& prefix) {
    std::ostringstream o;
    o << prefix << "kairoq " << p.build.version << ", commit " << p.build.commit << ", " << p.time_utc << "\n";
    if (!p.command.empty()) {
        o << prefix << "command:";
        for (const auto& c : p.command) o << " " << c;
        o << "\n";
    }
    for (const auto& d : p.inputs) o << prefix << "input " << d.path << " sha256 " << d.sha256 << "\n";
    for (const auto& d : p.board_sources) o << prefix << "board source " << d.path << " sha256 " << d.sha256 << "\n";
    return o.str();
}

}  // namespace kairoq
