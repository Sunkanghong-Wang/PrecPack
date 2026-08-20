#pragma once

#include "precpack/types.hpp"

#include <functional>
#include <vector>

class GRBEnv;

namespace precpack {

struct FullCapacityRemoval {
    int original_item = -1;
    std::vector<int> predecessor_original_items;
    std::vector<int> successor_original_items;
};

struct PreparedInstance {
    // Search-space transformations retain a mapping to original item indices.
    Instance search_instance;
    Assignment search_incumbent;
    std::vector<int> search_to_original;
    std::vector<int> original_to_search;
    int original_item_count = 0;
    int fixed_bin_offset = 0;
    std::vector<FullCapacityRemoval> full_capacity_removals;
    std::vector<std::vector<int>> fixed_prefix_bins;
    std::vector<std::vector<int>> fixed_suffix_bins;
    bool structured_preprocessing_enabled = false;
    bool reversed = false;
};

struct InitialBoundsResult {
    int lb1 = 0;
    int lb2 = 0;
    int lb3 = 0;
    int lb4 = 0;
    int lower_bound = 0;
    Assignment incumbent;
    bool preprocessing_reversed = false;
    bool initial_column_generation_attempted = false;
    bool initial_column_generation_converged = false;
    bool early_bbr_attempted = false;
    bool early_bbr_optimal = false;
    BbrStatistics early_bbr_statistics;
    PreparedInstance prepared;
};

[[nodiscard]] Assignment map_prepared_assignment_to_original(
    const PreparedInstance& prepared,
    const Assignment& search_assignment);

[[nodiscard]] int compute_initial_dff_lower_bound(
    const std::vector<int>& weights,
    int capacity);

[[nodiscard]] InitialBoundsResult compute_initial_bounds(
    const Instance& original,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics,
    const std::function<GRBEnv&()>& environment_provider);

}
