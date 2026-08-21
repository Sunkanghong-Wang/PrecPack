#pragma once

#include "precpack/types.hpp"

#include <cstdint>
#include <limits>

class GRBEnv;

namespace precpack {

struct BinIndexedRootBoundResult {
    bool attempted = false;
    bool completed = false;
    bool timed_out = false;
    int certified_lower_bound = 0;
    double lp_value = std::numeric_limits<double>::quiet_NaN();
    std::uint64_t column_count = 0;
};

[[nodiscard]] BinIndexedRootBoundResult run_bin_indexed_root_bound(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics);

}
