#pragma once

#include "precpack/cli.hpp"

#include <filesystem>
#include <optional>
#include <vector>

namespace precpack {

struct BatchCase {
    std::filesystem::path instance_path;
    std::optional<std::filesystem::path> graph_path;
};

[[nodiscard]] std::vector<BatchCase> collect_batch_cases(
    ProblemKind problem,
    const std::optional<std::filesystem::path>& input_path,
    const std::optional<std::filesystem::path>& graph_path,
    const std::filesystem::path& instance_root = "data/instances",
    const std::filesystem::path& graph_root = "data/bpp-gp-graphs");

[[nodiscard]] int run_batch(const CommandLineOptions& options);

}
