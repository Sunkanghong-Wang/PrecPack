#pragma once

#include "precpack/types.hpp"

class GRBEnv;

namespace precpack {

struct ColumnGenerationResult {
    bool attempted = false;
    bool converged = false;
    bool pricing_proven = false;
    bool timed_out = false;
    int integer_lower_bound = 0;
    double lp_value = 0.0;
};

[[nodiscard]] ColumnGenerationResult run_initial_bpp_column_generation(
    GRBEnv& environment,
    const Instance& strengthened_instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& global_deadline,
    Statistics& statistics);

}
