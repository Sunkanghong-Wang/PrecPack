#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
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
    std::string precedence_graph_structure;
    std::string order_strength;
    std::string item_weight_distribution;
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

struct LowerBounds {
    int capacity = 0;
    int precedence = 0;
    int combined = 0;
};

enum class SolveStatus {
    kNotSolved,
    kOptimal,
    kFeasible,
    kInfeasible,
    kTimeLimit,
    kStateLimit,
    kMemoryLimit,
    kError,
};

enum class ExactMethod {
    kBbr,
    // BPPC supplies a lower bound, never a BPP-GP incumbent.
    kBppc,
    kBpc,
    kMip,
};

enum class RootModelKind {
    kDirectPrecedence,
    kFixedK,
    kM,
};

enum class BbrRootCgMode {
    kNone,
    kDirectPrecedence,
    kFixedK,
    kM,
    // Switches from the position-free to the position-indexed root bound.
    kPriceAndSwitch,
};

struct Statistics {
    std::uint64_t explored_nodes = 0;
    std::uint64_t infeasible_nodes = 0;
    std::uint64_t integral_nodes = 0;
    std::uint64_t rf_branches = 0;
    std::uint64_t position_branches = 0;
    std::uint64_t phase_one_count = 0;
    std::uint64_t pricing_search_nodes = 0;
    std::uint64_t rmp_count = 0;
    std::uint64_t pricing_count = 0;
    std::uint64_t cg_count = 0;
    std::uint64_t cg_iterations = 0;
    std::uint64_t generated_columns = 0;
    std::uint64_t generated_precedence_rows = 0;
    std::uint64_t generated_sr_rows = 0;
    std::uint64_t precedence_check_count = 0;
    std::uint64_t precedence_check_cache_hits = 0;
    std::uint64_t initial_bdp_states = 0;
    std::uint64_t initial_bdp_transitions = 0;

    double preprocessing_seconds = 0.0;
    double lower_bound_seconds = 0.0;
    double upper_bound_seconds = 0.0;
    double rmp_seconds = 0.0;
    double pricing_seconds = 0.0;
    double cg_seconds = 0.0;
    double precedence_check_seconds = 0.0;
    double mip_seconds = 0.0;
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
    std::uint64_t certificate_scale_max = 0;
    double rmp_objective_multiplier = 1.0;
    double total_seconds = 0.0;
    double pricing_seconds = 0.0;
    double rmp_seconds = 0.0;
    double time_limit_seconds = 0.0;
    double initialization_seconds = 0.0;
    double initialization_time_limit_seconds = 0.0;
    int seed = 1;
    int run_order = 1;
};

struct BbrStatistics {
    bool attempted = false;
    bool timed_out = false;
    bool state_limited = false;
    bool memory_limited = false;
    bool parallel = false;
    bool shared_memory_saturated = false;
    bool reverse_direction = false;
    bool heuristic_phase_attempted = false;
    bool exact_phase_attempted = false;
    bool load_generation_truncated = false;
    bool item_dominance_enabled = false;
    bool generalized_item_dominance_enabled = false;
    bool paper_queue_order_enabled = false;
    bool bbr12_load_order_enabled = false;
    bool complete_dff_enabled = false;
    bool binlb_enabled = false;
    bool conflict_binlb_enabled = false;
    bool structured_preprocessing_enabled = false;
    std::string stop_reason;
    BbrRootCgMode root_cg_mode = BbrRootCgMode::kNone;

    std::uint64_t states_created = 0;
    std::uint64_t states_expanded = 0;
    std::uint64_t states_reopened = 0;
    std::uint64_t peak_open_states = 0;
    std::uint64_t loads_generated = 0;
    std::uint64_t load_search_nodes = 0;
    std::uint64_t truncated_load_states = 0;
    std::uint64_t heuristic_states_created = 0;
    std::uint64_t heuristic_states_expanded = 0;
    std::uint64_t heuristic_loads_generated = 0;
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
    std::uint64_t conflict_binlb_calls = 0;
    std::uint64_t conflict_binlb_completed = 0;
    std::uint64_t conflict_binlb_aborted = 0;
    std::uint64_t conflict_binlb_timeouts = 0;
    std::uint64_t conflict_binlb_node_limits = 0;
    std::uint64_t conflict_binlb_load_limits = 0;
    std::uint64_t conflict_binlb_item_skips = 0;
    std::uint64_t conflict_binlb_backoff_skips = 0;
    std::uint64_t conflict_binlb_budget_skips = 0;
    std::uint64_t conflict_binlb_bound_improvements = 0;
    std::uint64_t conflict_binlb_prunes = 0;
    std::uint64_t conflict_binlb_search_nodes = 0;
    std::uint64_t conflict_binlb_loads = 0;
    std::uint64_t conflict_binlb_memo_hits = 0;
    std::uint64_t conflict_binlb_ordinary_memo_hits = 0;
    std::uint64_t conflict_binlb_memo_entries = 0;
    std::uint64_t conflict_binlb_conflict_edges = 0;
    std::uint64_t conflict_binlb_node_limit = 0;
    std::uint64_t conflict_binlb_load_limit = 0;
    std::uint64_t conflict_binlb_memo_limit = 0;
    int conflict_binlb_max_items = 0;
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

    std::uint64_t peak_memory_bytes = 0;
    std::uint64_t state_limit = 0;
    std::uint64_t configured_state_limit = 0;
    std::uint64_t memory_limit_bytes = 0;
    std::uint64_t heuristic_load_limit = 0;
    bool initial_alns_enabled = true;
    bool initial_bdp_enabled = true;
    int seed = 1;
    int requested_threads = 1;
    int threads = 1;
    double time_limit_seconds = 0.0;
    double initialization_time_limit_seconds = 0.0;
    double search_seconds = 0.0;
    double heuristic_search_seconds = 0.0;
    double exact_search_seconds = 0.0;
    double parallel_split_seconds = 0.0;
    double generalized_item_dominance_seconds = 0.0;
    double binlb_seconds = 0.0;
    double binlb_call_time_limit_seconds = 0.0;
    double binlb_total_time_limit_seconds = 0.0;
    double conflict_binlb_seconds = 0.0;
    double conflict_binlb_call_time_limit_seconds = 0.0;
    double conflict_binlb_total_time_limit_seconds = 0.0;
};

struct BppcBoundStatistics {
    bool attempted = false;
    bool completed = false;
    bool timed_out = false;
    int certified_lower_bound = 0;
    int incumbent_value = std::numeric_limits<int>::max();
    int root_integer_lower_bound = 0;
    double root_lp_value = std::numeric_limits<double>::quiet_NaN();
    std::uint64_t root_column_count = 0;
    double total_seconds = 0.0;
};

struct Solution {
    SolveStatus status = SolveStatus::kNotSolved;
    std::string status_detail;
    bool optimal = false;
    int lb1 = 0;
    int lb2 = 0;
    int lb3 = 0;
    int lb4 = 0;
    int root_lp_lower_bound = 0;
    int initial_lower_bound = 0;
    int lower_bound = 0;
    int initial_upper_bound = std::numeric_limits<int>::max();
    int upper_bound = std::numeric_limits<int>::max();
    double root_lp_value = std::numeric_limits<double>::quiet_NaN();
    double mip_best_bound = std::numeric_limits<double>::quiet_NaN();
    double relative_gap = std::numeric_limits<double>::infinity();
    Assignment assignment;
    Statistics stats;
    ExactMethod exact_method = ExactMethod::kBbr;
    BbrStatistics bbr_stats;
    BppcBoundStatistics bppc_stats;
    RootModelKind root_model = RootModelKind::kDirectPrecedence;
    RootStatistics root_stats;
    int threads = 1;
};

struct Config {
    double time_limit_seconds = 300.0;
    double cg_time_limit_seconds = 0.0;
    double initialization_time_limit_seconds = 0.0;
    int seed = 1;
    // -1 uses available hardware concurrency; 1 preserves serial BBR.
    int threads = 1;
    bool enable_initial_alns = true;
    bool enable_initial_column_generation = true;
    bool run_column_generation = true;
    bool initialization_only = false;
    bool root_node_only = false;
    bool root_comparison_mode = false;
    bool set_covering_master = false;
    bool enable_sr_cuts = true;
    bool gurobi_log = false;
    ExactMethod exact_method = ExactMethod::kBbr;
    // State and memory limits are global across all BBR workers.
    std::uint64_t bbr_state_limit = 60'000'000ULL;
    std::uint64_t bbr_memory_limit_mb = 24ULL * 1024ULL;
    std::uint64_t bbr_heuristic_load_limit = 0U;
    bool bbr_enable_early_exact_probe = true;
    bool bbr_enable_initial_bdp = true;
    bool bbr_enable_bidirectional_bdp = false;
    bool bbr_enable_bbr12_mhh = false;
    bool bbr_enable_bbr12_mhh_portfolio = false;
    int bbr12_mhh_portfolio_max_items = 200;
    int bbr12_mhh_full_load_limit = 5'000;
    bool bbr_enable_initial_alns = false;
    bool bbr_enable_jackson = true;
    bool bbr_enable_bbr12_jackson = false;
    bool bbr_enable_no_successor = true;
    bool bbr_enable_superset_memory = true;
    bool bbr_enable_profile_dominance = true;
    bool bbr_enable_generalized_item_dominance = false;
    bool bbr_enable_paper_queue_order = true;
    bool bbr_enable_bbr12_load_order = true;
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
    bool bbr_enable_conflict_binlb = false;
    double bbr_conflict_binlb_call_time_limit_seconds = 1.0;
    double bbr_conflict_binlb_total_time_limit_seconds = 1.0;
    std::uint64_t bbr_conflict_binlb_node_limit = 250'000U;
    std::uint64_t bbr_conflict_binlb_load_limit = 50U;
    std::uint64_t bbr_conflict_binlb_memo_limit = 200'000U;
    int bbr_conflict_binlb_max_items = 400;
    BbrRootCgMode bbr_root_cg_mode = BbrRootCgMode::kPriceAndSwitch;
    double bbr_root_cg_time_limit_seconds = 5.0;
    bool bbr_initialization_mode = false;
    int max_precedence_rows_per_round = 64;
    int max_sr_cuts = 100;
    int max_cg_iterations = 10'000;
    int max_sr_rounds_per_node = 10;
    int max_columns_per_pricing = 8;
    int sr_enumeration_item_limit = 80;
    double reduced_cost_tolerance = 1e-8;
    double row_violation_tolerance = 1e-8;
    RootModelKind root_model = RootModelKind::kDirectPrecedence;
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
[[nodiscard]] const char* to_string(ExactMethod method) noexcept;
[[nodiscard]] const char* to_string(RootModelKind model) noexcept;
[[nodiscard]] const char* to_string(BbrRootCgMode mode) noexcept;

}
