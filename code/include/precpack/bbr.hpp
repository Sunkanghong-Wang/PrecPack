#pragma once

#include "precpack/initial_bounds.hpp"
#include "precpack/types.hpp"

namespace precpack {

struct BbrResult {
    bool attempted = false;
    bool optimal = false;
    bool timed_out = false;
    bool state_limited = false;
    bool memory_limited = false;
    int certified_lower_bound = 0;
    Assignment incumbent;
    BbrStatistics statistics;
};

[[nodiscard]] BbrResult run_branch_bound_remember(
    const PreparedInstance& prepared,
    int initial_lower_bound,
    const Config& config,
    Deadline& deadline);

}
