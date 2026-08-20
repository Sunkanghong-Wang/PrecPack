#include "precpack/solver.hpp"

#include "precpack/algorithms.hpp"
#include "precpack/bbr.hpp"
#include "precpack/build_config.hpp"
#if PRECPACK_HAS_GUROBI
#include "precpack/branch_price.hpp"
#include "precpack/gurobi_solver.hpp"
#endif
#include "precpack/initial_bounds.hpp"
#if PRECPACK_HAS_GUROBI
#include "precpack/root_column_generation.hpp"
#include <gurobi_c++.h>
#endif
#include "precpack/solver_profile.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

namespace precpack {
namespace {

#if PRECPACK_HAS_GUROBI
void accumulate_statistics(Statistics& target, const Statistics& source) {
    target.explored_nodes += source.explored_nodes;
    target.infeasible_nodes += source.infeasible_nodes;
    target.integral_nodes += source.integral_nodes;
    target.rf_branches += source.rf_branches;
    target.position_branches += source.position_branches;
    target.phase_one_count += source.phase_one_count;
    target.pricing_search_nodes += source.pricing_search_nodes;
    target.rmp_count += source.rmp_count;
    target.pricing_count += source.pricing_count;
    target.cg_count += source.cg_count;
    target.cg_iterations += source.cg_iterations;
    target.generated_columns += source.generated_columns;
    target.generated_precedence_rows += source.generated_precedence_rows;
    target.generated_sr_rows += source.generated_sr_rows;
    target.precedence_check_count += source.precedence_check_count;
    target.precedence_check_cache_hits +=
        source.precedence_check_cache_hits;
    target.initial_bdp_states += source.initial_bdp_states;
    target.initial_bdp_transitions += source.initial_bdp_transitions;
    target.preprocessing_seconds += source.preprocessing_seconds;
    target.lower_bound_seconds += source.lower_bound_seconds;
    target.upper_bound_seconds += source.upper_bound_seconds;
    target.rmp_seconds += source.rmp_seconds;
    target.pricing_seconds += source.pricing_seconds;
    target.cg_seconds += source.cg_seconds;
    target.precedence_check_seconds += source.precedence_check_seconds;
    target.mip_seconds += source.mip_seconds;
    target.initial_bdp_seconds += source.initial_bdp_seconds;
}

[[nodiscard]] RootStatistics direct_root_statistics(
    const BranchPriceResult& result,
    const Statistics& statistics,
    double elapsed_seconds) {
    RootStatistics root;
    root.attempted = result.attempted;
    root.completed = result.root_only_completed;
    root.timed_out = result.timed_out;
    root.certified_lower_bound = result.root_integer_lower_bound;
    root.lp_value = result.root_lp_value;
    root.column_count = result.root_column_count;
    root.generated_columns = statistics.generated_columns;
    root.iterations = statistics.cg_iterations;
    root.pricing_count = statistics.pricing_count;
    root.pricing_search_nodes = statistics.pricing_search_nodes;
    root.rmp_count = statistics.rmp_count;
    root.generated_precedence_rows =
        statistics.generated_precedence_rows;
    root.phase_one_count = statistics.phase_one_count;
    root.total_seconds = elapsed_seconds;
    root.pricing_seconds = statistics.pricing_seconds;
    root.rmp_seconds = statistics.rmp_seconds;
    return root;
}

void append_root_statistics(RootStatistics& target,
                            const RootStatistics& source) {
    const bool target_was_attempted = target.attempted;
    target.attempted = target.attempted || source.attempted;
    target.completed = source.attempted ? source.completed : target.completed;
    target.timed_out = source.attempted ? source.timed_out : target.timed_out;
    target.numerical_failure =
        target.numerical_failure || source.numerical_failure;
    target.certified_lower_bound = std::max(
        target.certified_lower_bound, source.certified_lower_bound);
    if (std::isfinite(source.lp_value)) {
        target.lp_value = std::isfinite(target.lp_value)
                              ? std::max(target.lp_value, source.lp_value)
                              : source.lp_value;
    } else if (!target_was_attempted) {
        target.lp_value = source.lp_value;
    }
    target.column_count += source.column_count;
    target.generated_columns += source.generated_columns;
    target.iterations += source.iterations;
    target.pricing_count += source.pricing_count;
    target.pricing_search_nodes += source.pricing_search_nodes;
    target.rmp_count += source.rmp_count;
    target.generated_precedence_rows += source.generated_precedence_rows;
    target.phase_one_count += source.phase_one_count;
    target.total_seconds += source.total_seconds;
    target.pricing_seconds += source.pricing_seconds;
    target.rmp_seconds += source.rmp_seconds;
}
#endif

void accumulate_preliminary_bbr_statistics(
    BbrStatistics& target,
    const BbrStatistics& source) noexcept {
    target.attempted = target.attempted || source.attempted;
    target.heuristic_phase_attempted =
        target.heuristic_phase_attempted || source.heuristic_phase_attempted;
    target.exact_phase_attempted =
        target.exact_phase_attempted || source.exact_phase_attempted;
    target.load_generation_truncated =
        target.load_generation_truncated || source.load_generation_truncated;
    target.item_dominance_enabled =
        target.item_dominance_enabled || source.item_dominance_enabled;
    target.generalized_item_dominance_enabled =
        target.generalized_item_dominance_enabled ||
        source.generalized_item_dominance_enabled;
    target.paper_queue_order_enabled =
        target.paper_queue_order_enabled || source.paper_queue_order_enabled;
    target.bbr12_load_order_enabled =
        target.bbr12_load_order_enabled || source.bbr12_load_order_enabled;
    target.complete_dff_enabled =
        target.complete_dff_enabled || source.complete_dff_enabled;
    target.binlb_enabled = target.binlb_enabled || source.binlb_enabled;
    target.conflict_binlb_enabled =
        target.conflict_binlb_enabled || source.conflict_binlb_enabled;
    target.structured_preprocessing_enabled =
        target.structured_preprocessing_enabled ||
        source.structured_preprocessing_enabled;

    target.states_created += source.states_created;
    target.states_expanded += source.states_expanded;
    target.states_reopened += source.states_reopened;
    target.loads_generated += source.loads_generated;
    target.load_search_nodes += source.load_search_nodes;
    target.truncated_load_states += source.truncated_load_states;
    target.heuristic_states_created += source.heuristic_states_created;
    target.heuristic_states_expanded += source.heuristic_states_expanded;
    target.heuristic_loads_generated += source.heuristic_loads_generated;
    target.forced_empty_transitions += source.forced_empty_transitions;
    target.hash_lookups += source.hash_lookups;
    target.hash_probes += source.hash_probes;
    target.exact_memory_prunes += source.exact_memory_prunes;
    target.profile_dominance_prunes += source.profile_dominance_prunes;
    target.superset_memory_prunes += source.superset_memory_prunes;
    target.bound_prunes += source.bound_prunes;
    target.nonmaximal_load_prunes += source.nonmaximal_load_prunes;
    target.jackson_prunes += source.jackson_prunes;
    target.generalized_item_dominance_search_nodes +=
        source.generalized_item_dominance_search_nodes;
    target.generalized_item_dominance_checks +=
        source.generalized_item_dominance_checks;
    target.generalized_item_dominance_prunes +=
        source.generalized_item_dominance_prunes;
    target.no_successor_prunes += source.no_successor_prunes;
    target.machine_bound_calls += source.machine_bound_calls;
    target.machine_bound_prunes += source.machine_bound_prunes;
    target.closure_bound_calls += source.closure_bound_calls;
    target.closure_bound_improvements += source.closure_bound_improvements;
    target.closure_bound_prunes += source.closure_bound_prunes;
    target.binlb_calls += source.binlb_calls;
    target.binlb_completed += source.binlb_completed;
    target.binlb_aborted += source.binlb_aborted;
    target.binlb_timeouts += source.binlb_timeouts;
    target.binlb_node_limits += source.binlb_node_limits;
    target.binlb_load_limits += source.binlb_load_limits;
    target.binlb_item_skips += source.binlb_item_skips;
    target.binlb_disabled_skips += source.binlb_disabled_skips;
    target.binlb_budget_skips += source.binlb_budget_skips;
    target.binlb_bound_improvements += source.binlb_bound_improvements;
    target.binlb_prunes += source.binlb_prunes;
    target.binlb_search_nodes += source.binlb_search_nodes;
    target.binlb_loads += source.binlb_loads;
    target.binlb_memo_hits += source.binlb_memo_hits;
    target.binlb_memo_entries = std::max(
        target.binlb_memo_entries, source.binlb_memo_entries);
    target.conflict_binlb_calls += source.conflict_binlb_calls;
    target.conflict_binlb_completed += source.conflict_binlb_completed;
    target.conflict_binlb_aborted += source.conflict_binlb_aborted;
    target.conflict_binlb_timeouts += source.conflict_binlb_timeouts;
    target.conflict_binlb_node_limits += source.conflict_binlb_node_limits;
    target.conflict_binlb_load_limits += source.conflict_binlb_load_limits;
    target.conflict_binlb_item_skips += source.conflict_binlb_item_skips;
    target.conflict_binlb_backoff_skips +=
        source.conflict_binlb_backoff_skips;
    target.conflict_binlb_budget_skips += source.conflict_binlb_budget_skips;
    target.conflict_binlb_bound_improvements +=
        source.conflict_binlb_bound_improvements;
    target.conflict_binlb_prunes += source.conflict_binlb_prunes;
    target.conflict_binlb_search_nodes += source.conflict_binlb_search_nodes;
    target.conflict_binlb_loads += source.conflict_binlb_loads;
    target.conflict_binlb_memo_hits += source.conflict_binlb_memo_hits;
    target.conflict_binlb_ordinary_memo_hits +=
        source.conflict_binlb_ordinary_memo_hits;
    target.conflict_binlb_memo_entries = std::max(
        target.conflict_binlb_memo_entries,
        source.conflict_binlb_memo_entries);
    target.conflict_binlb_conflict_edges = std::max(
        target.conflict_binlb_conflict_edges,
        source.conflict_binlb_conflict_edges);
    target.incumbent_updates += source.incumbent_updates;
    target.peak_open_states =
        std::max(target.peak_open_states, source.peak_open_states);
    target.peak_memory_bytes =
        std::max(target.peak_memory_bytes, source.peak_memory_bytes);
    target.structured_search_items =
        std::max(target.structured_search_items,
                 source.structured_search_items);
    target.structured_items_removed =
        std::max(target.structured_items_removed,
                 source.structured_items_removed);
    target.structured_full_capacity_items =
        std::max(target.structured_full_capacity_items,
                 source.structured_full_capacity_items);
    target.structured_prefix_bins =
        std::max(target.structured_prefix_bins,
                 source.structured_prefix_bins);
    target.structured_suffix_bins =
        std::max(target.structured_suffix_bins,
                 source.structured_suffix_bins);
    target.structured_fixed_bins =
        std::max(target.structured_fixed_bins,
                 source.structured_fixed_bins);
    target.generalized_item_dominance_pairs =
        std::max(target.generalized_item_dominance_pairs,
                 source.generalized_item_dominance_pairs);
    target.dff_transform_count =
        std::max(target.dff_transform_count, source.dff_transform_count);
    target.search_seconds += source.search_seconds;
    target.heuristic_search_seconds += source.heuristic_search_seconds;
    target.exact_search_seconds += source.exact_search_seconds;
    target.generalized_item_dominance_seconds +=
        source.generalized_item_dominance_seconds;
    target.binlb_seconds += source.binlb_seconds;
    target.conflict_binlb_seconds += source.conflict_binlb_seconds;
}

}

Solution solve(const Instance& instance, const Config& requested_config) {
    Config config = requested_config;
#if !PRECPACK_HAS_GUROBI
    if (config.exact_method != ExactMethod::kBbr ||
        config.root_comparison_mode) {
        throw std::invalid_argument(
            "this build excludes Gurobi; only the exact BBR solver is available");
    }
    config.enable_initial_column_generation = false;
    config.run_column_generation = false;
    config.bbr_root_cg_mode = BbrRootCgMode::kNone;
#endif
    if (config.time_limit_seconds <= 0.0) {
        throw std::invalid_argument("time limit must be positive");
    }
    if (config.threads == 0 || config.threads < -1) {
        throw std::invalid_argument("threads must be -1 or a positive integer");
    }
    const auto solve_start = std::chrono::steady_clock::now();
    Deadline deadline(config.time_limit_seconds);
    Solution solution;
    solution.threads = resolve_thread_count(config.threads);
    solution.exact_method = config.exact_method;
    solution.root_model = config.root_model;

#if PRECPACK_HAS_GUROBI
    bool optional_root_backend_failed = false;
    std::unique_ptr<GRBEnv> environment;
    const auto environment_provider = [&]() -> GRBEnv& {
        if (!environment) {
            environment = std::make_unique<GRBEnv>(true);
            environment->set(GRB_IntParam_OutputFlag, config.gurobi_log ? 1 : 0);
            environment->set(GRB_IntParam_Threads, 1);
            environment->set(GRB_IntParam_Seed, static_cast<int>(config.seed));
            environment->start();
        }
        return *environment;
    };
#else
    const std::function<GRBEnv&()> environment_provider;
#endif
    Config initialization_config = config;
    if (config.exact_method == ExactMethod::kBppc) {
        initialization_config.enable_initial_column_generation = false;
        initialization_config.enable_initial_alns = false;
    }
    if (config.exact_method == ExactMethod::kBbr &&
        !config.root_comparison_mode) {
        initialization_config.enable_initial_column_generation = false;
        initialization_config.enable_initial_alns =
            config.initialization_only ? config.enable_initial_alns
                                       : config.bbr_enable_initial_alns;
        initialization_config.bbr_initialization_mode =
            !config.initialization_only;
        if (initialization_config.initialization_time_limit_seconds <= 0.0) {
            initialization_config.initialization_time_limit_seconds =
                std::min(3.0, 0.1 * config.time_limit_seconds);
        }
    }
    if (config.exact_method == ExactMethod::kBbr &&
        !config.root_comparison_mode) {
        constexpr std::uint64_t kMegabyte = 1024U * 1024U;
        solution.bbr_stats.configured_state_limit = config.bbr_state_limit;
        solution.bbr_stats.state_limit = config.bbr_state_limit;
        solution.bbr_stats.memory_limit_bytes =
            config.bbr_memory_limit_mb * kMegabyte;
        solution.bbr_stats.heuristic_load_limit =
            config.bbr_heuristic_load_limit;
        solution.bbr_stats.initial_alns_enabled =
            initialization_config.enable_initial_alns;
        solution.bbr_stats.initial_bdp_enabled =
            initialization_config.bbr_enable_initial_bdp;
        solution.bbr_stats.requested_threads = config.threads;
        solution.bbr_stats.threads = solution.threads;
        solution.bbr_stats.parallel = solution.threads > 1;
        solution.bbr_stats.item_dominance_enabled =
            config.bbr_enable_jackson;
        solution.bbr_stats.generalized_item_dominance_enabled =
            config.bbr_enable_generalized_item_dominance;
        solution.bbr_stats.paper_queue_order_enabled =
            config.bbr_enable_paper_queue_order;
        solution.bbr_stats.bbr12_load_order_enabled =
            config.bbr_enable_bbr12_load_order;
        solution.bbr_stats.complete_dff_enabled =
            config.bbr_enable_complete_dff;
        solution.bbr_stats.binlb_enabled =
            config.bbr_enable_binlb || config.bbr_enable_conflict_binlb;
        solution.bbr_stats.binlb_call_time_limit_seconds =
            config.bbr_binlb_call_time_limit_seconds;
        solution.bbr_stats.binlb_total_time_limit_seconds =
            config.bbr_binlb_total_time_limit_seconds;
        solution.bbr_stats.binlb_node_limit = config.bbr_binlb_node_limit;
        solution.bbr_stats.binlb_load_limit = config.bbr_binlb_load_limit;
        solution.bbr_stats.binlb_memo_limit = config.bbr_binlb_memo_limit;
        solution.bbr_stats.binlb_max_items = config.bbr_binlb_max_items;
        solution.bbr_stats.conflict_binlb_enabled =
            config.bbr_enable_conflict_binlb;
        solution.bbr_stats.conflict_binlb_call_time_limit_seconds =
            config.bbr_conflict_binlb_call_time_limit_seconds;
        solution.bbr_stats.conflict_binlb_total_time_limit_seconds =
            config.bbr_conflict_binlb_total_time_limit_seconds;
        solution.bbr_stats.conflict_binlb_node_limit =
            config.bbr_conflict_binlb_node_limit;
        solution.bbr_stats.conflict_binlb_load_limit =
            config.bbr_conflict_binlb_load_limit;
        solution.bbr_stats.conflict_binlb_memo_limit =
            config.bbr_conflict_binlb_memo_limit;
        solution.bbr_stats.conflict_binlb_max_items =
            config.bbr_conflict_binlb_max_items;
        solution.bbr_stats.root_cg_mode = config.bbr_root_cg_mode;
        solution.bbr_stats.seed = config.seed;
        solution.bbr_stats.time_limit_seconds = config.time_limit_seconds;
        solution.bbr_stats.initialization_time_limit_seconds =
            initialization_config.initialization_time_limit_seconds;
    }
    const InitialBoundsResult initial = compute_initial_bounds(
        instance, initialization_config, deadline, solution.stats,
        environment_provider);
    if (initial.early_bbr_attempted) {
        accumulate_preliminary_bbr_statistics(
            solution.bbr_stats, initial.early_bbr_statistics);
        solution.stats.explored_nodes +=
            initial.early_bbr_statistics.states_expanded;
        solution.stats.infeasible_nodes +=
            initial.early_bbr_statistics.bound_prunes;
    }
    const double initialization_elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      solve_start)
            .count();
    solution.lb1 = initial.lb1;
    solution.lb2 = initial.lb2;
    solution.lb3 = initial.lb3;
    solution.lb4 = initial.lb4;
    solution.initial_lower_bound = initial.lower_bound;
    solution.lower_bound = initial.lower_bound;
    solution.assignment = initial.incumbent;
    solution.initial_upper_bound = solution.assignment.bin_count;
    solution.upper_bound = solution.assignment.bin_count;
    if (config.exact_method == ExactMethod::kBbr &&
        !config.root_comparison_mode) {
        solution.bbr_stats.reverse_direction = initial.prepared.reversed;
        solution.bbr_stats.structured_preprocessing_enabled =
            initial.prepared.structured_preprocessing_enabled;
        solution.bbr_stats.structured_search_items =
            static_cast<std::uint64_t>(
                initial.prepared.search_instance.size());
        solution.bbr_stats.structured_items_removed =
            initial.prepared.original_item_count >
                    initial.prepared.search_instance.size()
                ? static_cast<std::uint64_t>(
                      initial.prepared.original_item_count -
                      initial.prepared.search_instance.size())
                : 0U;
        solution.bbr_stats.structured_full_capacity_items =
            initial.prepared.full_capacity_removals.size();
        solution.bbr_stats.structured_prefix_bins =
            initial.prepared.fixed_prefix_bins.size();
        solution.bbr_stats.structured_suffix_bins =
            initial.prepared.fixed_suffix_bins.size();
        solution.bbr_stats.structured_fixed_bins =
            static_cast<std::uint64_t>(initial.prepared.fixed_bin_offset);
    }

    std::string diagnostic;
    if (!check_assignment(instance, solution.assignment, &diagnostic)) {
        throw std::logic_error("initial assignment check failed: " + diagnostic);
    }

#if PRECPACK_HAS_GUROBI
    const int lower_bound_before_optional_root = solution.lower_bound;
    const int root_lower_bound_before_optional_root =
        solution.root_lp_lower_bound;
    const double root_value_before_optional_root = solution.root_lp_value;
    const RootModelKind root_model_before_optional_root = solution.root_model;
    const RootStatistics root_statistics_before_optional_root =
        solution.root_stats;
    const Statistics statistics_before_optional_root = solution.stats;
    try {
        if (config.exact_method == ExactMethod::kBbr &&
            !config.root_comparison_mode && !config.initialization_only &&
            config.bbr_root_cg_mode != BbrRootCgMode::kNone &&
            solution.lower_bound < solution.upper_bound &&
            !deadline.expired()) {
            const double configured_budget =
                config.bbr_root_cg_time_limit_seconds > 0.0
                    ? config.bbr_root_cg_time_limit_seconds
                    : deadline.remaining_seconds();
            const double root_budget =
                std::min(configured_budget, deadline.remaining_seconds());
            solution.root_stats.time_limit_seconds = root_budget;
            solution.root_stats.initialization_seconds = initialization_elapsed;
            solution.root_stats.initialization_time_limit_seconds =
                initialization_config.initialization_time_limit_seconds;
            solution.root_stats.seed = config.seed;

            const auto run_root_model = [&](RootModelKind kind, double budget,
                                            int current_lower_bound) {
                RootStatistics root;
                const int residual_lower_bound = std::max(
                    1, current_lower_bound - initial.prepared.fixed_bin_offset);
                if (budget <= 0.0 || deadline.expired()) {
                    root.attempted = true;
                    root.timed_out = true;
                    root.certified_lower_bound = residual_lower_bound;
                    root.certified_lower_bound +=
                        initial.prepared.fixed_bin_offset;
                    return root;
                }
                Deadline root_deadline(
                    std::min(budget, deadline.remaining_seconds()));
                Statistics root_algorithm_statistics;
                if (kind == RootModelKind::kDirectPrecedence) {
                    Config root_config = config;
                    root_config.root_node_only = true;
                    root_config.root_comparison_mode = false;
                    root_config.set_covering_master = true;
                    root_config.enable_sr_cuts = false;
                    const auto root_start = std::chrono::steady_clock::now();
                    const BranchPriceResult result = run_branch_price_and_cut(
                        environment_provider(),
                        initial.prepared.search_instance,
                        initial.prepared.search_incumbent, residual_lower_bound,
                        root_config, root_deadline, root_algorithm_statistics);
                    root = direct_root_statistics(
                        result, root_algorithm_statistics,
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - root_start)
                            .count());
                } else {
                    root = run_position_free_root_column_generation(
                        environment_provider(),
                        initial.prepared.search_instance,
                        initial.prepared.search_incumbent, residual_lower_bound,
                        kind, config, root_deadline, root_algorithm_statistics);
                }
                if (root.attempted) {
                    root.certified_lower_bound +=
                        initial.prepared.fixed_bin_offset;
                    if (std::isfinite(root.lp_value)) {
                        root.lp_value += initial.prepared.fixed_bin_offset;
                    }
                }
                accumulate_statistics(solution.stats,
                                      root_algorithm_statistics);
                return root;
            };

            if (config.bbr_root_cg_mode == BbrRootCgMode::kPriceAndSwitch) {
                constexpr int kAdaptiveRootItemLimit = 128;
                constexpr double kMProbeSeconds = 0.15;
                constexpr double kDirectCapSeconds = 3.0;
                const Instance& search_instance =
                    initial.prepared.search_instance;
                const int maximum_separation = std::accumulate(
                    search_instance.arcs.begin(), search_instance.arcs.end(), 0,
                    [](int value, const Arc& arc) {
                        return std::max(value, arc.separation);
                    });
                solution.root_model = RootModelKind::kM;
                RootStatistics first;
                if (search_instance.size() <= kAdaptiveRootItemLimit) {
                    first =
                        run_root_model(RootModelKind::kM,
                                       std::min(root_budget, kMProbeSeconds),
                                       solution.lower_bound);
                    append_root_statistics(solution.root_stats, first);
                    solution.lower_bound =
                        std::min(solution.upper_bound,
                                 std::max(solution.lower_bound,
                                          first.certified_lower_bound));
                }
                const bool use_direct =
                    first.completed && maximum_separation >= 2 &&
                    solution.upper_bound - solution.lower_bound >= 2 &&
                    !deadline.expired();
                if (use_direct) {
                    const double remaining_root_budget = std::min(
                        {std::max(0.0, root_budget - first.total_seconds),
                         kDirectCapSeconds, deadline.remaining_seconds()});
                    RootStatistics second = run_root_model(
                        RootModelKind::kDirectPrecedence, remaining_root_budget,
                        solution.lower_bound);
                    append_root_statistics(solution.root_stats, second);
                    solution.lower_bound =
                        std::min(solution.upper_bound,
                                 std::max(solution.lower_bound,
                                          second.certified_lower_bound));
                    solution.root_model = RootModelKind::kDirectPrecedence;
                }
            } else {
                RootModelKind kind = RootModelKind::kDirectPrecedence;
                if (config.bbr_root_cg_mode == BbrRootCgMode::kFixedK) {
                    kind = RootModelKind::kFixedK;
                } else if (config.bbr_root_cg_mode == BbrRootCgMode::kM) {
                    kind = RootModelKind::kM;
                }
                solution.root_model = kind;
                RootStatistics root =
                    run_root_model(kind, root_budget, solution.lower_bound);
                append_root_statistics(solution.root_stats, root);
                solution.lower_bound = std::min(
                    solution.upper_bound,
                    std::max(solution.lower_bound, root.certified_lower_bound));
            }
            solution.root_stats.time_limit_seconds = root_budget;
            solution.root_stats.initialization_seconds = initialization_elapsed;
            solution.root_stats.initialization_time_limit_seconds =
                initialization_config.initialization_time_limit_seconds;
            solution.root_stats.seed = config.seed;
            solution.root_lp_lower_bound = solution.lower_bound;
            solution.root_lp_value = solution.root_stats.lp_value;
        }
    } catch (const GRBException&) {
        optional_root_backend_failed = true;
        environment.reset();
        config.bbr_root_cg_mode = BbrRootCgMode::kNone;
        solution.lower_bound = lower_bound_before_optional_root;
        solution.root_lp_lower_bound = root_lower_bound_before_optional_root;
        solution.root_lp_value = root_value_before_optional_root;
        solution.root_model = root_model_before_optional_root;
        solution.root_stats = root_statistics_before_optional_root;
        solution.stats = statistics_before_optional_root;
    }
#endif

    if (solution.lower_bound >= solution.upper_bound) {
        solution.lower_bound = solution.upper_bound;
        solution.root_lp_lower_bound = solution.lower_bound;
        solution.root_lp_value = static_cast<double>(solution.lower_bound);
        solution.optimal = true;
        solution.status = SolveStatus::kOptimal;
        solution.status_detail =
            initial.early_bbr_optimal
                ? "early exact BBR probe proved optimality"
                : (solution.root_stats.attempted
                       ? "root column generation closed the initial gap"
                       : "initial lower and upper bounds coincide");
        if (config.exact_method == ExactMethod::kBbr &&
            !config.root_comparison_mode) {
            solution.bbr_stats.stop_reason =
                initial.early_bbr_optimal
                    ? "EARLY_BBR_OPTIMAL"
                    : (solution.root_stats.attempted ? "ROOT_CG_OPTIMAL"
                                                     : "INITIAL_BOUNDS_OPTIMAL");
        }
        if (config.root_comparison_mode) {
            solution.root_stats.completed = true;
            solution.root_stats.certified_lower_bound = solution.lower_bound;
            solution.root_stats.lp_value =
                static_cast<double>(solution.lower_bound);
            solution.root_stats.time_limit_seconds =
                config.cg_time_limit_seconds > 0.0
                    ? config.cg_time_limit_seconds
                    : deadline.remaining_seconds();
            solution.root_stats.initialization_seconds =
                initialization_elapsed;
            solution.root_stats.initialization_time_limit_seconds =
                config.initialization_time_limit_seconds;
            solution.root_stats.seed = config.seed;
        }
    } else if (config.initialization_only) {
        solution.root_lp_lower_bound = solution.lower_bound;
        solution.status = deadline.expired() ? SolveStatus::kTimeLimit
                                             : SolveStatus::kFeasible;
        solution.status_detail = deadline.expired()
                                     ? "global time limit reached during initialization; "
                                       "exact BPC and compact MIP were skipped"
                                     : "initialization-only mode; exact BPC and compact "
                                       "MIP were skipped";
    } else if (config.exact_method == ExactMethod::kBbr &&
               !config.root_comparison_mode) {
        if (config.root_node_only) {
            throw std::invalid_argument(
                "root-only mode is defined only for branch-price-and-cut");
        }
        if (!solution.root_stats.attempted) {
            solution.root_lp_lower_bound = solution.initial_lower_bound;
        }
        if (deadline.expired()) {
            solution.status = SolveStatus::kTimeLimit;
            solution.status_detail =
                "global time limit reached during BBR initialization";
            solution.bbr_stats.timed_out = true;
            solution.bbr_stats.stop_reason = "TIME_LIMIT_INITIALIZATION";
        } else {
            Config search_config = config;
            search_config.enable_initial_alns =
                initialization_config.enable_initial_alns;
            const BbrResult bbr = run_branch_bound_remember(
                initial.prepared, solution.lower_bound, search_config,
                deadline);
            BbrStatistics combined_statistics = bbr.statistics;
            combined_statistics.initial_alns_enabled =
                initialization_config.enable_initial_alns;
            combined_statistics.initial_bdp_enabled =
                initialization_config.bbr_enable_initial_bdp;
            if (initial.early_bbr_attempted) {
                accumulate_preliminary_bbr_statistics(
                    combined_statistics, initial.early_bbr_statistics);
            }
            solution.bbr_stats = std::move(combined_statistics);
            solution.assignment = map_prepared_assignment_to_original(
                initial.prepared, bbr.incumbent);
            solution.upper_bound = solution.assignment.bin_count;
            solution.lower_bound = std::min(
                solution.upper_bound,
                std::max(solution.lower_bound, bbr.certified_lower_bound));
            solution.mip_best_bound =
                static_cast<double>(bbr.certified_lower_bound);
            solution.optimal = bbr.optimal;
            solution.stats.explored_nodes += bbr.statistics.states_expanded;
            solution.stats.infeasible_nodes += bbr.statistics.bound_prunes;
            if (bbr.optimal) {
                solution.status = SolveStatus::kOptimal;
                solution.status_detail =
                    "branch, bound, and remember proved optimality";
            } else if (bbr.state_limited) {
                solution.status = SolveStatus::kStateLimit;
                solution.status_detail =
                    "branch, bound, and remember reached its global state limit";
            } else if (bbr.memory_limited) {
                solution.status = SolveStatus::kMemoryLimit;
                solution.status_detail =
                    "branch, bound, and remember reached its safe memory limit";
            } else {
                solution.status = SolveStatus::kTimeLimit;
                solution.status_detail =
                    "branch, bound, and remember reached the global time limit";
            }
        }
#if PRECPACK_HAS_GUROBI
    } else {
        GRBEnv& solver_environment = environment_provider();

        if (config.root_comparison_mode && !deadline.expired()) {
            const double budget = config.cg_time_limit_seconds > 0.0
                                      ? std::min(config.cg_time_limit_seconds,
                                                 deadline.remaining_seconds())
                                      : deadline.remaining_seconds();
            Deadline root_deadline(budget);
            Statistics root_algorithm_statistics;
            if (config.root_model == RootModelKind::kDirectPrecedence) {
                Config root_config = config;
                root_config.root_node_only = true;
                root_config.set_covering_master = true;
                root_config.enable_sr_cuts = false;
                const auto root_start = std::chrono::steady_clock::now();
                const BranchPriceResult root = run_branch_price_and_cut(
                    solver_environment, instance, solution.assignment,
                    solution.lower_bound, root_config, root_deadline,
                    root_algorithm_statistics);
                solution.root_stats = direct_root_statistics(
                    root, root_algorithm_statistics,
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - root_start)
                        .count());
            } else {
                solution.root_stats =
                    run_position_free_root_column_generation(
                        solver_environment, instance, solution.assignment,
                        solution.lower_bound, config.root_model, config,
                        root_deadline, root_algorithm_statistics);
            }
            solution.root_stats.time_limit_seconds = budget;
            solution.root_stats.initialization_seconds =
                initialization_elapsed;
            solution.root_stats.initialization_time_limit_seconds =
                config.initialization_time_limit_seconds;
            solution.root_stats.seed = config.seed;
            accumulate_statistics(solution.stats, root_algorithm_statistics);
            solution.root_lp_value = solution.root_stats.lp_value;
            solution.root_lp_lower_bound =
                solution.root_stats.certified_lower_bound;
            solution.lower_bound = std::min(
                solution.upper_bound,
                std::max(solution.lower_bound,
                         solution.root_stats.certified_lower_bound));
            solution.mip_best_bound =
                static_cast<double>(solution.root_stats.certified_lower_bound);
            solution.optimal = solution.lower_bound >= solution.upper_bound;
            if (solution.optimal) {
                solution.status = SolveStatus::kOptimal;
                solution.status_detail =
                    "set-covering root relaxation closed the initial gap";
            } else if (solution.root_stats.completed) {
                solution.status = SolveStatus::kFeasible;
                solution.status_detail =
                    "set-covering root relaxation completed";
            } else if (solution.root_stats.numerical_failure) {
                solution.status = SolveStatus::kFeasible;
                solution.status_detail =
                    "root relaxation stopped without a numerical certificate";
            } else {
                solution.status = SolveStatus::kTimeLimit;
                solution.status_detail =
                    "root relaxation reached its time or iteration limit";
            }
        } else if (config.exact_method == ExactMethod::kBppc &&
                   !deadline.expired()) {
            const auto bppc_start = std::chrono::steady_clock::now();
            const int capacity_lower_bound = static_cast<int>(
                instance.total_weight / instance.capacity +
                (instance.total_weight % instance.capacity != 0 ? 1 : 0));
            const BppcBoundResult bppc = run_bppc_branch_price_bound(
                solver_environment, instance, solution.assignment,
                capacity_lower_bound, config, deadline, solution.stats);
            solution.bppc_stats.attempted = bppc.attempted;
            solution.bppc_stats.completed = bppc.optimal;
            solution.bppc_stats.timed_out = bppc.timed_out;
            solution.bppc_stats.certified_lower_bound =
                bppc.certified_lower_bound;
            solution.bppc_stats.incumbent_value = bppc.incumbent_value;
            solution.bppc_stats.root_integer_lower_bound =
                bppc.root_integer_lower_bound;
            solution.bppc_stats.root_lp_value = bppc.root_lp_value;
            solution.bppc_stats.root_column_count = bppc.root_column_count;
            solution.bppc_stats.total_seconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - bppc_start)
                    .count();
            solution.root_model = RootModelKind::kM;
            solution.root_lp_value = bppc.root_lp_value;
            solution.root_lp_lower_bound = std::max(
                solution.initial_lower_bound,
                bppc.root_integer_lower_bound);
            solution.lower_bound = std::min(
                solution.upper_bound,
                std::max(solution.lower_bound,
                         bppc.certified_lower_bound));
            solution.mip_best_bound =
                static_cast<double>(solution.lower_bound);
            solution.optimal =
                solution.lower_bound >= solution.upper_bound;
            if (solution.optimal) {
                solution.status = SolveStatus::kOptimal;
                solution.status_detail =
                    "the exact BPPC relaxation matched the validated "
                    "BPP-GP incumbent";
            } else if (bppc.optimal) {
                solution.status = SolveStatus::kFeasible;
                solution.status_detail =
                    "the exact BPPC relaxation completed; its integer "
                    "optimum is only a BPP-GP lower bound";
            } else {
                solution.status = SolveStatus::kTimeLimit;
                solution.status_detail =
                    "the BPPC branch-price-and-cut tree stopped before "
                    "proving its relaxation optimum";
            }
        } else if (config.exact_method == ExactMethod::kBpc &&
                   config.run_column_generation && !deadline.expired()) {
            const BranchPriceResult bpc =
                config.root_model == RootModelKind::kM
                    ? run_m_branch_price_and_cut(
                          solver_environment, instance, solution.assignment,
                          solution.lower_bound, config, deadline,
                          solution.stats)
                    : run_branch_price_and_cut(
                          solver_environment, instance, solution.assignment,
                          solution.lower_bound, config, deadline,
                          solution.stats);
            solution.root_lp_value = bpc.root_lp_value;
            solution.root_lp_lower_bound = bpc.root_integer_lower_bound;
            solution.assignment = bpc.incumbent;
            solution.upper_bound = bpc.incumbent.bin_count;
            solution.lower_bound = std::min(
                solution.upper_bound,
                std::max(solution.lower_bound, bpc.certified_lower_bound));
            solution.mip_best_bound =
                static_cast<double>(bpc.certified_lower_bound);
            solution.optimal = bpc.optimal;
            if (bpc.optimal) {
                solution.status = SolveStatus::kOptimal;
                solution.status_detail =
                    "branch-price-and-cut tree proved optimality";
            } else if (bpc.root_only_completed) {
                solution.status = SolveStatus::kFeasible;
                solution.status_detail =
                    "root-only mode completed a certified BPC root solve";
            } else {
                solution.status = SolveStatus::kTimeLimit;
                solution.status_detail =
                    "branch-price-and-cut stopped before proof (time or "
                    "iteration limit)";
            }
        } else if (config.exact_method == ExactMethod::kMip &&
                   !deadline.expired()) {
            solution.root_lp_lower_bound = solution.lower_bound;
            const auto mip_start = std::chrono::steady_clock::now();
            const MipResult mip = solve_compact_mip(
                solver_environment, instance, solution.assignment, solution.lower_bound,
                config, deadline);
            solution.stats.mip_seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              mip_start)
                    .count();
            solution.stats.explored_nodes = mip.explored_nodes;
            solution.mip_best_bound = mip.best_bound;
            solution.lower_bound =
                std::min(solution.upper_bound,
                         std::max(solution.lower_bound, mip.certified_lower_bound));
            if (mip.has_incumbent &&
                mip.assignment.bin_count <= solution.assignment.bin_count) {
                solution.assignment = mip.assignment;
                solution.upper_bound = mip.assignment.bin_count;
            }
            solution.lower_bound = std::min(solution.lower_bound, solution.upper_bound);
            solution.status = mip.status;
            solution.optimal = mip.optimal;
            if (mip.optimal) {
                solution.lower_bound = solution.upper_bound;
                solution.status_detail = "Gurobi compact MIP proved optimality";
            } else if (mip.status == SolveStatus::kTimeLimit) {
                solution.status_detail =
                    "time limit reached with a validated incumbent and certified bound";
            } else {
                solution.status_detail = "Gurobi compact MIP stopped before proof";
            }
        } else {
            solution.root_lp_lower_bound = solution.lower_bound;
            solution.status = SolveStatus::kTimeLimit;
            solution.status_detail =
                "global time limit reached before exact optimization";
        }
    }
#else
    }
#endif

    if (!check_assignment(instance, solution.assignment, &diagnostic)) {
        throw std::logic_error("final assignment check failed: " + diagnostic);
    }
    if (solution.optimal) {
        solution.status = SolveStatus::kOptimal;
        solution.lower_bound = solution.upper_bound;
        solution.relative_gap = 0.0;
    } else {
        solution.relative_gap = solution.upper_bound > 0
                                    ? static_cast<double>(solution.upper_bound -
                                                          solution.lower_bound) /
                                          solution.upper_bound
                                    : std::numeric_limits<double>::infinity();
    }
#if PRECPACK_HAS_GUROBI
    if (optional_root_backend_failed) {
        solution.status_detail +=
            "; optional Gurobi root strengthening was unavailable";
    }
#endif
    solution.stats.total_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start)
            .count();
    return solution;
}

}
