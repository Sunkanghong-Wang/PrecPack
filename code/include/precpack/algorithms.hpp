#pragma once

#include "precpack/types.hpp"

#include <random>
#include <string>

namespace precpack {

[[nodiscard]] LowerBounds compute_simple_lower_bounds(const Instance& instance);

[[nodiscard]] Assignment construct_initial_assignment(
    const Instance& instance,
    int known_lower_bound,
    std::mt19937& random,
    int randomized_trials = 20,
    const Deadline* deadline = nullptr);

[[nodiscard]] bool check_assignment(
    const Instance& instance,
    const Assignment& assignment,
    std::string* diagnostic = nullptr);

}
