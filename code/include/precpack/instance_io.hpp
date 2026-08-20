#pragma once

#include "precpack/types.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace precpack {

[[nodiscard]] Instance read_instance(
    const std::filesystem::path& alb_path,
    const std::optional<std::filesystem::path>& graph_path,
    const std::string& problem_type,
    int instance_id = -1);

}
