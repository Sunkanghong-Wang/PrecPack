#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace precpack {

inline constexpr double kEpsilon = 1e-9;

struct Item {
    int original_index = -1;
    int weight = 0;
};

struct Arc {
    int from = -1;
    int to = -1;
    int separation = 0;
};

struct Instance {
    std::string problem_type = "BPP-GP";
    std::string order_strength;
    int id = -1;
    int capacity = 0;
    std::vector<Item> items;
    std::vector<Arc> arcs;

    std::int64_t total_weight = 0;
    std::vector<std::vector<int>> predecessors;
    std::vector<std::vector<int>> successors;
    std::vector<std::vector<std::pair<int, int>>> predecessor_arcs;
    std::vector<std::vector<std::pair<int, int>>> successor_arcs;
    std::vector<int> topological_order;
    std::vector<int> front;
    std::vector<int> back;
    std::vector<int> longest_separation;

    [[nodiscard]] int size() const noexcept {
        return static_cast<int>(items.size());
    }

    [[nodiscard]] int separation(int i, int j) const noexcept {
        return longest_separation[static_cast<std::size_t>(i) * items.size() +
                                  static_cast<std::size_t>(j)];
    }

    void initialize();
};

struct Assignment {
    std::vector<int> bin_of_item;
    int bin_count = 0;

    [[nodiscard]] bool complete() const noexcept;
};

enum class SolveStatus {
    kNotSolved,
    kOptimal,
    kTimeLimit,
    kStateLimit,
    kMemoryLimit,
};

struct Statistics {
    std::uint64_t explored_nodes = 0;
    std::uint64_t infeasible_nodes = 0;
    std::uint64_t phase_one_count = 0;
    std::uint64_t pricing_search_nodes = 0;
    std::uint64_t rmp_count = 0;
    std::uint64_t pricing_count = 0;
    std::uint64_t cg_count = 0;
    std::uint64_t cg_iterations = 0;
    std::uint64_t generated_columns = 0;
    std::uint64_t generated_precedence_rows = 0;
    std::uint64_t initial_bdp_states = 0;
    std::uint64_t initial_bdp_transitions = 0;

    double preprocessing_seconds = 0.0;
    double lower_bound_seconds = 0.0;
    double upper_bound_seconds = 0.0;
    double rmp_seconds = 0.0;
    double pricing_seconds = 0.0;
    double cg_seconds = 0.0;
    double initial_bdp_seconds = 0.0;
    double total_seconds = 0.0;
};

struct RootStatistics {
    bool attempted = false;
    bool completed = false;
    bool timed_out = false;
    bool numerical_failure = false;
    int certified_lower_bound = 0;
    double lp_value = std::numeric_limits<double>::quiet_NaN();
    std::uint64_t column_count = 0;
    std::uint64_t generated_columns = 0;
    std::uint64_t iterations = 0;
    std::uint64_t pricing_count = 0;
    std::uint64_t pricing_search_nodes = 0;
    std::uint64_t rmp_count = 0;
    std::uint64_t generated_precedence_rows = 0;
    std::uint64_t phase_one_count = 0;
    std::uint64_t fixed_point_scale_min = 0;
    std::uint64_t fixed_point_scale_max = 0;
    std::uint64_t certificate_scale_min = 0;
    double rmp_objective_multiplier = 1.0;
    double total_seconds = 0.0;
    double pricing_seconds = 0.0;
    double rmp_seconds = 0.0;
};

struct BbrStatistics {
    bool attempted = false;
    bool timed_out = false;
    bool state_limited = false;
    bool memory_limited = false;
    bool parallel = false;
    bool shared_memory_saturated = false;
    bool reverse_direction = false;
    bool exact_phase_attempted = false;
    bool item_dominance_enabled = false;
    bool generalized_item_dominance_enabled = false;
    bool paper_queue_order_enabled = false;
    bool complete_dff_enabled = false;
    bool binlb_enabled = false;
    bool structured_preprocessing_enabled = false;
    std::uint64_t states_created = 0;
    std::uint64_t states_expanded = 0;
    std::uint64_t states_reopened = 0;
    std::uint64_t peak_open_states = 0;
    std::uint64_t loads_generated = 0;
    std::uint64_t load_search_nodes = 0;
    std::uint64_t forced_empty_transitions = 0;
    std::uint64_t hash_lookups = 0;
    std::uint64_t hash_probes = 0;
    std::uint64_t exact_memory_prunes = 0;
    std::uint64_t profile_dominance_prunes = 0;
    std::uint64_t superset_memory_prunes = 0;
    std::uint64_t bound_prunes = 0;
    std::uint64_t nonmaximal_load_prunes = 0;
    std::uint64_t jackson_prunes = 0;
    std::uint64_t generalized_item_dominance_pairs = 0;
    std::uint64_t generalized_item_dominance_search_nodes = 0;
    std::uint64_t generalized_item_dominance_checks = 0;
    std::uint64_t generalized_item_dominance_prunes = 0;
    std::uint64_t no_successor_prunes = 0;
    std::uint64_t machine_bound_calls = 0;
    std::uint64_t machine_bound_prunes = 0;
    std::uint64_t closure_bound_calls = 0;
    std::uint64_t closure_bound_improvements = 0;
    std::uint64_t closure_bound_prunes = 0;
    std::uint64_t binlb_calls = 0;
    std::uint64_t binlb_completed = 0;
    std::uint64_t binlb_aborted = 0;
    std::uint64_t binlb_timeouts = 0;
    std::uint64_t binlb_node_limits = 0;
    std::uint64_t binlb_load_limits = 0;
    std::uint64_t binlb_item_skips = 0;
    std::uint64_t binlb_disabled_skips = 0;
    std::uint64_t binlb_budget_skips = 0;
    std::uint64_t binlb_bound_improvements = 0;
    std::uint64_t binlb_prunes = 0;
    std::uint64_t binlb_search_nodes = 0;
    std::uint64_t binlb_loads = 0;
    std::uint64_t binlb_memo_hits = 0;
    std::uint64_t binlb_memo_entries = 0;
    std::uint64_t binlb_node_limit = 0;
    std::uint64_t binlb_load_limit = 0;
    std::uint64_t binlb_memo_limit = 0;
    int binlb_max_items = 0;
    std::uint64_t incumbent_updates = 0;
    std::uint64_t dff_transform_count = 0;
    std::uint64_t structured_items_removed = 0;
    std::uint64_t structured_full_capacity_items = 0;
    std::uint64_t structured_prefix_bins = 0;
    std::uint64_t structured_suffix_bins = 0;
    std::uint64_t structured_fixed_bins = 0;
    std::uint64_t structured_search_items = 0;
    std::uint64_t parallel_tasks_generated = 0;
    std::uint64_t parallel_tasks_completed = 0;
    std::uint64_t parallel_tasks_stolen = 0;
    std::uint64_t shared_exact_memory_prunes = 0;
    std::uint64_t shared_profile_dominance_prunes = 0;
    std::uint64_t shared_superset_memory_prunes = 0;
    std::uint64_t parallel_shared_peak_memory_bytes = 0;
    std::uint64_t parallel_worker_peak_memory_bytes = 0;
    std::uint64_t parallel_task_memory_bytes = 0;
    std::uint64_t parallel_worker_tasks_min = 0;
    std::uint64_t parallel_worker_tasks_max = 0;
    std::uint64_t parallel_worker_expanded_min = 0;
    std::uint64_t parallel_worker_expanded_max = 0;
    std::uint64_t parallel_initial_work_min = 0;
    std::uint64_t parallel_initial_work_max = 0;

    std::uint64_t peak_memory_bytes = 0;
    std::uint64_t state_limit = 0;
    std::uint64_t configured_state_limit = 0;
    std::uint64_t memory_limit_bytes = 0;
    bool initial_bdp_enabled = true;
    int seed = 1;
    int requested_threads = 1;
    int threads = 1;
    double time_limit_seconds = 0.0;
    double initialization_time_limit_seconds = 0.0;
    double search_seconds = 0.0;
    double exact_search_seconds = 0.0;
    double parallel_split_seconds = 0.0;
    double parallel_worker_busy_seconds_sum = 0.0;
    double parallel_worker_busy_seconds_min = 0.0;
    double parallel_worker_busy_seconds_max = 0.0;
    double generalized_item_dominance_seconds = 0.0;
    double binlb_seconds = 0.0;
    double binlb_call_time_limit_seconds = 0.0;
    double binlb_total_time_limit_seconds = 0.0;
};

struct Solution {
    SolveStatus status = SolveStatus::kNotSolved;
    bool optimal = false;
    bool gurobi_runtime_required = false;
    int lower_bound = 0;
    int upper_bound = std::numeric_limits<int>::max();
    double relative_gap = std::numeric_limits<double>::infinity();
    Assignment assignment;
    Statistics stats;
    BbrStatistics bbr_stats;
    RootStatistics root_stats;
    int threads = 1;
};

struct Config {
    double time_limit_seconds = 300.0;
    double initialization_time_limit_seconds = 0.0;
    int seed = 1;
    // -1 uses available hardware concurrency; 1 preserves serial BBR.
    int threads = 1;
    bool require_gurobi_runtime = false;
    // State and memory limits are global across all BBR workers.
    std::uint64_t bbr_state_limit = 60'000'000ULL;
    std::uint64_t bbr_memory_limit_mb = 24ULL * 1024ULL;
    bool bbr_enable_early_exact_probe = true;
    bool bbr_enable_initial_bdp = true;
    bool bbr_enable_bbr12_mhh = false;
    bool bbr_enable_bbr12_mhh_portfolio = false;
    int bbr12_mhh_portfolio_max_items = 200;
    int bbr12_mhh_full_load_limit = 5'000;
    bool bbr_enable_jackson = true;
    bool bbr_enable_bbr12_jackson = false;
    bool bbr_enable_no_successor = true;
    bool bbr_enable_superset_memory = true;
    bool bbr_enable_profile_dominance = true;
    bool bbr_enable_generalized_item_dominance = false;
    bool bbr_enable_paper_queue_order = true;
    bool bbr_enable_complete_dff = true;
    int bbr_dff_transform_limit = 0;
    bool bbr_enable_closure_bound = true;
    bool bbr_enable_structured_preprocessing = true;
    bool bbr_enable_binlb = false;
    double bbr_binlb_call_time_limit_seconds = 1.0;
    double bbr_binlb_total_time_limit_seconds = 1.0;
    std::uint64_t bbr_binlb_node_limit = 250'000U;
    std::uint64_t bbr_binlb_load_limit = 50U;
    std::uint64_t bbr_binlb_memo_limit = 200'000U;
    int bbr_binlb_max_items = 400;
    bool bbr_enable_root_strengthening = true;
    double bbr_root_cg_time_limit_seconds = 5.0;
    bool bbr_initialization_mode = false;
    int max_precedence_rows_per_round = 64;
    int max_cg_iterations = 10'000;
    int max_columns_per_pricing = 8;
    double reduced_cost_tolerance = 1e-8;
    double row_violation_tolerance = 1e-8;
};

class Deadline {
public:
    explicit Deadline(double seconds)
        : begin_(Clock::now()),
          end_(begin_ + std::chrono::duration_cast<Clock::duration>(
                            std::chrono::duration<double>(seconds))) {}

    [[nodiscard]] bool expired() const noexcept {
        return Clock::now() >= end_;
    }

    [[nodiscard]] double elapsed_seconds() const noexcept {
        return std::chrono::duration<double>(Clock::now() - begin_).count();
    }

    [[nodiscard]] double remaining_seconds() const noexcept {
        return std::max(0.0,
                        std::chrono::duration<double>(end_ - Clock::now()).count());
    }

    [[nodiscard]] std::chrono::steady_clock::time_point end_time() const noexcept {
        return end_;
    }

private:
    using Clock = std::chrono::steady_clock;
    Clock::time_point begin_;
    Clock::time_point end_;
};

[[nodiscard]] const char* to_string(SolveStatus status) noexcept;

}
