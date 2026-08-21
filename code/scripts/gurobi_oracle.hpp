#pragma once

#include "precpack/types.hpp"

class GRBEnv;

namespace precpack::test {

enum class OracleStatus {
    kNotSolved,
    kOptimal,
    kFeasible,
    kTimeLimit,
    kError,
};

struct MipResult {
    OracleStatus status = OracleStatus::kNotSolved;
    bool has_incumbent = false;
    bool optimal = false;
    int certified_lower_bound = 0;
    double best_bound = 0.0;
    Assignment assignment;
    std::uint64_t explored_nodes = 0;
};

[[nodiscard]] MipResult solve_compact_mip(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline);

}
