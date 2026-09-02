#include "precpack/solver.hpp"

#include "precpack/algorithms.hpp"
#include "precpack/bbr.hpp"
#include "precpack/build_config.hpp"
#include "precpack/initial_bounds.hpp"
#if PRECPACK_HAS_GUROBI
#include "root_column_generation.hpp"
#include <gurobi_c++.h>
#endif
#include "precpack/solver_profile.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace precpack {
namespace {

#if PRECPACK_HAS_GUROBI
void accumulate_statistics(Statistics& target, const Statistics& source) {
    target.explored_nodes += source.explored_nodes;
    target.infeasible_nodes += source.infeasible_nodes;
    target.pricing_search_nodes += source.pricing_search_nodes;
    target.rmp_count += source.rmp_count;
    target.pricing_count += source.pricing_count;
    target.cg_count += source.cg_count;
    target.cg_iterations += source.cg_iterations;
    target.generated_columns += source.generated_columns;
    target.initial_bdp_states += source.initial_bdp_states;
    target.initial_bdp_transitions += source.initial_bdp_transitions;
    target.preprocessing_seconds += source.preprocessing_seconds;
    target.lower_bound_seconds += source.lower_bound_seconds;
    target.upper_bound_seconds += source.upper_bound_seconds;
    target.rmp_seconds += source.rmp_seconds;
    target.pricing_seconds += source.pricing_seconds;
    target.cg_seconds += source.cg_seconds;
    target.initial_bdp_seconds += source.initial_bdp_seconds;
}
#endif

void accumulate_preliminary_bbr_statistics(
    BbrStatistics& target,
    const BbrStatistics& source) noexcept {
    target.attempted = target.attempted || source.attempted;
    target.exact_phase_attempted =
        target.exact_phase_attempted || source.exact_phase_attempted;
    target.item_dominance_enabled =
        target.item_dominance_enabled || source.item_dominance_enabled;
    target.generalized_item_dominance_enabled =
        target.generalized_item_dominance_enabled ||
        source.generalized_item_dominance_enabled;
    target.paper_queue_order_enabled =
        target.paper_queue_order_enabled || source.paper_queue_order_enabled;
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
    target.binlb_target_hits += source.binlb_target_hits;
    target.binlb_conflict_edges = std::max(
        target.binlb_conflict_edges, source.binlb_conflict_edges);
    target.binlb_conflict_calls += source.binlb_conflict_calls;
    target.binlb_conflict_completed += source.binlb_conflict_completed;
    target.binlb_conflict_search_nodes +=
        source.binlb_conflict_search_nodes;
    target.binlb_conflict_loads += source.binlb_conflict_loads;
    target.binlb_conflict_memo_hits += source.binlb_conflict_memo_hits;
    target.binlb_ordinary_memo_hits += source.binlb_ordinary_memo_hits;
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
    target.exact_search_seconds += source.exact_search_seconds;
    target.generalized_item_dominance_seconds +=
        source.generalized_item_dominance_seconds;
    target.binlb_seconds += source.binlb_seconds;
}

}

Solution solve(const Instance& instance, const Config& requested_config) {
    Config config = requested_config;
#if !PRECPACK_HAS_GUROBI
    if (config.require_gurobi_runtime) {
        throw std::runtime_error(
            "this run requires Gurobi, but this build excludes it");
    }
    config.bbr_enable_root_strengthening = false;
#endif
    if (config.time_limit_seconds <= 0.0) {
        throw std::invalid_argument("time limit must be positive");
    }
    const auto solve_start = std::chrono::steady_clock::now();
    Deadline deadline(config.time_limit_seconds);
    Solution solution;
    solution.gurobi_runtime_required = config.require_gurobi_runtime;

#if PRECPACK_HAS_GUROBI
    std::unique_ptr<GRBEnv> environment;
    const auto environment_provider = [&]() -> GRBEnv& {
        if (!environment) {
            environment = std::make_unique<GRBEnv>(true);
            environment->set(GRB_IntParam_OutputFlag, 0);
            environment->set(GRB_IntParam_Threads, 1);
            environment->set(GRB_IntParam_Seed, static_cast<int>(config.seed));
            environment->start();
        }
        return *environment;
    };
#endif
    Config initialization_config = config;
    initialization_config.bbr_initialization_mode = true;
    if (initialization_config.initialization_time_limit_seconds <= 0.0) {
        initialization_config.initialization_time_limit_seconds =
            std::min(3.0, 0.1 * config.time_limit_seconds);
    }
    constexpr std::uint64_t kMegabyte = 1024U * 1024U;
    solution.bbr_stats.memory_limit_bytes =
        config.bbr_memory_limit_mb * kMegabyte;
    solution.bbr_stats.initial_bdp_enabled =
        initialization_config.bbr_enable_initial_bdp;
    solution.bbr_stats.item_dominance_enabled = config.bbr_enable_jackson;
    solution.bbr_stats.generalized_item_dominance_enabled =
        config.bbr_enable_generalized_item_dominance;
    solution.bbr_stats.paper_queue_order_enabled =
        config.bbr_enable_paper_queue_order;
    solution.bbr_stats.complete_dff_enabled = config.bbr_enable_complete_dff;
    solution.bbr_stats.binlb_enabled = config.bbr_enable_binlb;
    solution.bbr_stats.conflict_binlb_enabled =
        config.bbr_enable_binlb && config.bbr_enable_conflict_binlb;
    solution.bbr_stats.binlb_call_time_limit_seconds =
        config.bbr_binlb_call_time_limit_seconds;
    solution.bbr_stats.binlb_total_time_limit_seconds =
        config.bbr_binlb_total_time_limit_seconds;
    solution.bbr_stats.binlb_node_limit = config.bbr_binlb_node_limit;
    solution.bbr_stats.binlb_load_limit = config.bbr_binlb_load_limit;
    solution.bbr_stats.binlb_memo_limit = config.bbr_binlb_memo_limit;
    solution.bbr_stats.binlb_max_items = config.bbr_binlb_max_items;
    solution.bbr_stats.conflict_binlb_call_time_limit_seconds =
        config.bbr_conflict_binlb_call_time_limit_seconds;
    solution.bbr_stats.conflict_binlb_node_limit =
        config.bbr_conflict_binlb_node_limit;
    solution.bbr_stats.seed = config.seed;
    solution.bbr_stats.time_limit_seconds = config.time_limit_seconds;
    solution.bbr_stats.initialization_time_limit_seconds =
        initialization_config.initialization_time_limit_seconds;
    const InitialBoundsResult initial = compute_initial_bounds(
        instance, initialization_config, deadline, solution.stats);
    if (initial.early_bbr_attempted) {
        accumulate_preliminary_bbr_statistics(
            solution.bbr_stats, initial.early_bbr_statistics);
        solution.stats.explored_nodes +=
            initial.early_bbr_statistics.states_expanded;
        solution.stats.infeasible_nodes +=
            initial.early_bbr_statistics.bound_prunes;
    }
    solution.lower_bound = initial.lower_bound;
    solution.assignment = initial.incumbent;
    solution.upper_bound = solution.assignment.bin_count;
    solution.bbr_stats.reverse_direction = initial.prepared.reversed;
    solution.bbr_stats.structured_preprocessing_enabled =
        initial.prepared.structured_preprocessing_enabled;
    solution.bbr_stats.structured_search_items =
        static_cast<std::uint64_t>(initial.prepared.search_instance.size());
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

    std::string diagnostic;
    if (!check_assignment(instance, solution.assignment, &diagnostic)) {
        throw std::logic_error("initial assignment check failed: " + diagnostic);
    }

#if PRECPACK_HAS_GUROBI
    const int lower_bound_before_optional_root = solution.lower_bound;
    const RootStatistics root_statistics_before_optional_root =
        solution.root_stats;
    const Statistics statistics_before_optional_root = solution.stats;
    try {
        if (config.bbr_enable_root_strengthening &&
            solution.lower_bound < solution.upper_bound &&
            initial.prepared.search_instance.size() <=
                kMaximumPositionFreeRootItems &&
            !deadline.expired()) {
            const auto run_root_model = [&](double budget,
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
                root = run_position_free_root_column_generation(
                    environment_provider(), initial.prepared.search_instance,
                    initial.prepared.search_incumbent, residual_lower_bound,
                    config, root_deadline, root_algorithm_statistics);
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

            const Instance& search_instance = initial.prepared.search_instance;
            constexpr double kMaximumRemainingTimeFraction = 0.05;
            constexpr double kRootSecondsPerItem = 0.0015;
            constexpr double kMinimumRootSeconds = 0.05;
            constexpr double kMaximumRootSeconds = 0.2;
            const double remaining_seconds = deadline.remaining_seconds();
            const double item_scaled_budget = std::max(
                kMinimumRootSeconds,
                kRootSecondsPerItem *
                    static_cast<double>(search_instance.size()));
            const double root_budget = std::min(
                {item_scaled_budget, kMaximumRootSeconds,
                 kMaximumRemainingTimeFraction * remaining_seconds});
            if (root_budget >= kMinimumRootSeconds) {
                RootStatistics root = run_root_model(
                    root_budget, solution.lower_bound);
                solution.root_stats = root;
                solution.lower_bound = std::min(
                    solution.upper_bound,
                    std::max(solution.lower_bound,
                             root.certified_lower_bound));
            }
        }
    } catch (const GRBException&) {
        if (config.require_gurobi_runtime) {
            throw;
        }
        environment.reset();
        config.bbr_enable_root_strengthening = false;
        solution.lower_bound = lower_bound_before_optional_root;
        solution.root_stats = root_statistics_before_optional_root;
        solution.stats = statistics_before_optional_root;
    }
#endif

    if (solution.lower_bound >= solution.upper_bound) {
        solution.lower_bound = solution.upper_bound;
        solution.optimal = true;
        solution.status = SolveStatus::kOptimal;
    } else {
        if (deadline.expired()) {
            solution.status = SolveStatus::kTimeLimit;
            solution.bbr_stats.timed_out = true;
        } else {
            const BbrResult bbr = run_branch_bound_remember(
                initial.prepared, solution.lower_bound, config,
                deadline);
            BbrStatistics combined_statistics = bbr.statistics;
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
            solution.optimal = bbr.optimal;
            solution.stats.explored_nodes += bbr.statistics.states_expanded;
            solution.stats.infeasible_nodes += bbr.statistics.bound_prunes;
            if (bbr.optimal) {
                solution.status = SolveStatus::kOptimal;
            } else if (bbr.memory_limited) {
                solution.status = SolveStatus::kMemoryLimit;
            } else if (bbr.timed_out || deadline.expired()) {
                solution.status = SolveStatus::kTimeLimit;
            } else {
                throw std::logic_error(
                    "exact BBR terminated without a resource limit");
            }
        }
    }

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
    solution.stats.total_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start)
            .count();
    return solution;
}

void verify_gurobi_runtime() {
#if PRECPACK_HAS_GUROBI
    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.set(GRB_IntParam_Threads, 1);
    environment.start();
#else
    throw std::runtime_error(
        "this run requires Gurobi, but this build excludes it");
#endif
}

}
