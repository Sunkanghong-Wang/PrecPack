#pragma once

#include "precpack/batch.hpp"

#include <filesystem>

namespace precpack::internal {

[[nodiscard]] double effective_batch_time_limit(
    const CommandLineOptions& options,
    const BatchCase& batch_case,
    const std::filesystem::path& instance_root);

}
