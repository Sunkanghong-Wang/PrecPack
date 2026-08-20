#pragma once

#include "precpack/solver_profile.hpp"

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string_view>

namespace precpack {

struct CommandLineOptions {
    std::filesystem::path instance_path;
    std::optional<std::filesystem::path> graph_path;
    std::filesystem::path output_directory = "results";
    ProblemKind problem = ProblemKind::kBppGp;
    double time_limit_seconds = 300.0;
    std::uint64_t memory_limit_mb = 24ULL * 1024ULL;
    int threads = 1;
    bool show_help = false;
};

[[nodiscard]] CommandLineOptions parse_command_line(
    int argc,
    char* const argv[]);

void print_help(std::ostream& output, std::string_view executable);

}
