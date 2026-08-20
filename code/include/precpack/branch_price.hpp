#pragma once

#include "precpack/types.hpp"

class GRBEnv;

namespace precpack {

struct BranchPriceResult {
    bool attempted = false;
    bool optimal = false;
    bool timed_out = false;
    bool root_only_completed = false;
    int certified_lower_bound = 0;
    int root_integer_lower_bound = 0;
    double root_lp_value = 0.0;
    std::uint64_t root_column_count = 0;
    Assignment incumbent;
};

struct BppcBoundResult {
    bool attempted = false;
    bool optimal = false;
    bool timed_out = false;
    int certified_lower_bound = 0;
    int incumbent_value = 0;  // BPPC value; not feasible for BPP-GP in general.
    int root_integer_lower_bound = 0;
    double root_lp_value = 0.0;
    std::uint64_t root_column_count = 0;
};

[[nodiscard]] BranchPriceResult run_branch_price_and_cut(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics);

[[nodiscard]] BranchPriceResult run_m_branch_price_and_cut(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics);

[[nodiscard]] BppcBoundResult run_bppc_branch_price_bound(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& original_feasible_incumbent,
    int bppc_lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics);

}
