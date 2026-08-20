#pragma once

#include "precpack/types.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

namespace precpack {

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

}
