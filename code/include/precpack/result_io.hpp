#pragma once

#include "precpack/types.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

namespace precpack {

[[nodiscard]] std::string make_instance_set(
    const std::filesystem::path& instance_path,
    const std::optional<std::filesystem::path>& graph_path);

void append_result_csv(const std::filesystem::path& path,
                       std::string_view instance_set,
                       std::string_view instance_name,
                       const Instance& instance,
                       const Solution& solution);

void validate_result_csv(const std::filesystem::path& path);

void write_assignment(const std::filesystem::path& path,
                      const Instance& instance,
                      const Solution& solution);

[[nodiscard]] bool has_nonempty_solution(
    const std::filesystem::path& path);

}
