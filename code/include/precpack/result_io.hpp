#pragma once

#include "precpack/types.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace precpack {

struct ResultReference {
    std::string instance_key;
    std::string problem;
    double time_limit_seconds = 0.0;
    int threads = 0;
    std::uint64_t memory_limit_mb = 0U;
    bool gurobi_enabled = false;
    bool gurobi_required = false;
    std::filesystem::path solution_file;
};

[[nodiscard]] std::string make_instance_key(
    const std::filesystem::path& instance_path,
    const std::optional<std::filesystem::path>& graph_path);

void append_result_csv(const std::filesystem::path& path,
                       std::string_view instance_key,
                       const std::filesystem::path& instance_path,
                       const std::optional<std::filesystem::path>& graph_path,
                       const Instance& instance,
                       const Solution& solution,
                       const std::filesystem::path& assignment_path);

void write_assignment(const std::filesystem::path& path,
                      const Instance& instance,
                      const Solution& solution);

[[nodiscard]] std::vector<ResultReference> read_result_references(
    const std::filesystem::path& path);

void require_unused_instance_key(const std::filesystem::path& path,
                                 std::string_view instance_key);

}
