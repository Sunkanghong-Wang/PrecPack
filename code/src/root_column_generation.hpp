#pragma once

#include "precpack/types.hpp"

class GRBEnv;

namespace precpack {

inline constexpr int kMaximumPositionFreeRootItems = 100;

[[nodiscard]] RootStatistics run_position_free_root_column_generation(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics);

}
