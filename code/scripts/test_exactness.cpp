#include "precpack/algorithms.hpp"
#include "precpack/bbr.hpp"
#include "precpack/bin_packing_bound.hpp"
#include "precpack/build_config.hpp"
#include "precpack/dff.hpp"
#include "precpack/initial_bounds.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#if PRECPACK_HAS_GUROBI
#include "bin_indexed_root_bound.hpp"
#include "gurobi_oracle.hpp"
#include "root_column_generation.hpp"

#include <gurobi_c++.h>
#endif

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

[[nodiscard]] precpack::Instance make_instance(
    std::vector<int> weights,
    int capacity,
    std::vector<precpack::Arc> arcs = {},
    std::string problem_type = "BPP-GP") {
    precpack::Instance instance;
    instance.problem_type = std::move(problem_type);
    instance.capacity = capacity;
    for (std::size_t item = 0; item < weights.size(); ++item) {
        instance.items.push_back({
            static_cast<int>(item), weights[item]});
    }
    instance.arcs = std::move(arcs);
    instance.initialize();
    return instance;
}

[[nodiscard]] int elementary_lower_bound(
    const precpack::Instance& instance) {
    const int capacity_bound = static_cast<int>(
        (instance.total_weight + instance.capacity - 1) /
        instance.capacity);
    int precedence_bound = 1;
    for (int item = 0; item < instance.size(); ++item) {
        precedence_bound = std::max(
            precedence_bound,
            instance.front[static_cast<std::size_t>(item)] +
                instance.back[static_cast<std::size_t>(item)] + 1);
    }
    return std::max(capacity_bound, precedence_bound);
}

[[nodiscard]] precpack::Assignment greedy_assignment(
    const precpack::Instance& instance) {
    precpack::Assignment assignment;
    assignment.bin_of_item.assign(
        static_cast<std::size_t>(instance.size()), -1);
    std::vector<int> loads;
    for (const int item : instance.topological_order) {
        int earliest = 0;
        for (const auto& [predecessor, separation] :
             instance.predecessor_arcs[static_cast<std::size_t>(item)]) {
            earliest = std::max(
                earliest,
                assignment.bin_of_item[
                    static_cast<std::size_t>(predecessor)] + separation);
        }
        int bin = earliest;
        while (true) {
            if (bin >= static_cast<int>(loads.size())) {
                loads.resize(static_cast<std::size_t>(bin + 1), 0);
            }
            if (loads[static_cast<std::size_t>(bin)] +
                    instance.items[static_cast<std::size_t>(item)].weight <=
                instance.capacity) {
                break;
            }
            ++bin;
        }
        assignment.bin_of_item[static_cast<std::size_t>(item)] = bin;
        loads[static_cast<std::size_t>(bin)] +=
            instance.items[static_cast<std::size_t>(item)].weight;
    }
    assignment.bin_count = 1 + *std::max_element(
        assignment.bin_of_item.begin(), assignment.bin_of_item.end());
    return assignment;
}

[[nodiscard]] bool find_assignment_with_horizon(
    const precpack::Instance& instance,
    int horizon,
    precpack::Assignment* solution) {
    std::vector<int> bin_of_item(
        static_cast<std::size_t>(instance.size()), -1);
    std::vector<int> loads(static_cast<std::size_t>(horizon), 0);
    std::vector<std::int64_t> remaining_weight(
        static_cast<std::size_t>(instance.size() + 1), 0);
    for (int position = instance.size() - 1; position >= 0; --position) {
        const int item = instance.topological_order[
            static_cast<std::size_t>(position)];
        remaining_weight[static_cast<std::size_t>(position)] =
            remaining_weight[static_cast<std::size_t>(position + 1)] +
            instance.items[static_cast<std::size_t>(item)].weight;
    }

    const std::function<bool(int)> search = [&](int position) {
        if (position == instance.size()) {
            precpack::Assignment candidate;
            candidate.bin_of_item = bin_of_item;
            candidate.bin_count = 1 + *std::max_element(
                bin_of_item.begin(), bin_of_item.end());
            std::string diagnostic;
            if (!precpack::check_assignment(
                    instance, candidate, &diagnostic)) {
                throw std::logic_error(
                    "brute-force oracle produced an invalid assignment: " +
                    diagnostic);
            }
            *solution = std::move(candidate);
            return true;
        }

        std::int64_t free_capacity = 0;
        for (const int load : loads) {
            free_capacity += instance.capacity - load;
        }
        if (remaining_weight[static_cast<std::size_t>(position)] >
            free_capacity) {
            return false;
        }

        const int item = instance.topological_order[
            static_cast<std::size_t>(position)];
        int earliest = 0;
        for (const auto& [predecessor, separation] :
             instance.predecessor_arcs[static_cast<std::size_t>(item)]) {
            earliest = std::max(
                earliest,
                bin_of_item[static_cast<std::size_t>(predecessor)] +
                    separation);
        }
        for (int bin = earliest; bin < horizon; ++bin) {
            const int weight =
                instance.items[static_cast<std::size_t>(item)].weight;
            if (loads[static_cast<std::size_t>(bin)] + weight >
                instance.capacity) {
                continue;
            }
            bin_of_item[static_cast<std::size_t>(item)] = bin;
            loads[static_cast<std::size_t>(bin)] += weight;
            if (search(position + 1)) {
                return true;
            }
            loads[static_cast<std::size_t>(bin)] -= weight;
            bin_of_item[static_cast<std::size_t>(item)] = -1;
        }
        return false;
    };
    return search(0);
}

[[nodiscard]] precpack::Assignment brute_force_optimum(
    const precpack::Instance& instance) {
    const precpack::Assignment upper = greedy_assignment(instance);
    for (int horizon = elementary_lower_bound(instance);
         horizon <= upper.bin_count; ++horizon) {
        precpack::Assignment assignment;
        if (find_assignment_with_horizon(instance, horizon, &assignment)) {
            return assignment;
        }
    }
    throw std::logic_error("brute-force oracle failed to recover a feasible upper bound");
}

[[nodiscard]] precpack::Config exact_config(
    precpack::ProblemKind problem = precpack::ProblemKind::kBppGp) {
    precpack::Config config =
        precpack::make_solver_config(problem, 10.0, 256, 1);
    config.bbr_enable_root_strengthening = false;
    return config;
}

[[nodiscard]] precpack::PreparedInstance make_identity_prepared(
    const precpack::Instance& instance,
    precpack::Assignment incumbent) {
    precpack::PreparedInstance prepared;
    prepared.search_instance = instance;
    prepared.search_incumbent = std::move(incumbent);
    prepared.search_to_original.resize(
        static_cast<std::size_t>(instance.size()));
    std::iota(prepared.search_to_original.begin(),
              prepared.search_to_original.end(), 0);
    prepared.original_to_search = prepared.search_to_original;
    prepared.original_item_count = instance.size();
    return prepared;
}

[[nodiscard]] precpack::Assignment spaced_assignment(
    const precpack::Instance& instance,
    int spacing = 4) {
    precpack::Assignment assignment;
    assignment.bin_of_item.assign(
        static_cast<std::size_t>(instance.size()), -1);
    int position = 0;
    for (const int item : instance.topological_order) {
        assignment.bin_of_item[static_cast<std::size_t>(item)] =
            spacing * position++;
    }
    assignment.bin_count = spacing * (instance.size() - 1) + 1;
    std::string diagnostic;
    if (!precpack::check_assignment(instance, assignment, &diagnostic)) {
        throw std::logic_error(
            "spaced test incumbent is invalid: " + diagnostic);
    }
    return assignment;
}

void require_solver_matches_oracle(
    const precpack::Instance& instance,
    precpack::ProblemKind problem) {
    const precpack::Assignment oracle = brute_force_optimum(instance);
    const precpack::Solution solution =
        precpack::solve(instance, exact_config(problem));
    std::string diagnostic;
    require(solution.status == precpack::SolveStatus::kOptimal &&
                solution.optimal &&
                solution.lower_bound == oracle.bin_count &&
                solution.upper_bound == oracle.bin_count &&
                precpack::check_assignment(
                    instance, solution.assignment, &diagnostic),
            "PrecPack disagreed with the brute-force oracle: expected " +
                std::to_string(oracle.bin_count) + ", received [" +
                std::to_string(solution.lower_bound) + "," +
                std::to_string(solution.upper_bound) + "]: " + diagnostic);
}

void test_assignment_checker() {
    const precpack::Instance instance =
        make_instance({1, 1}, 10, {{0, 1, 2}});
    const precpack::Assignment valid{{0, 2}, 3};
    std::string diagnostic;
    require(precpack::check_assignment(instance, valid, &diagnostic),
            "checker rejected a valid assignment: " + diagnostic);
    precpack::Assignment invalid = valid;
    invalid.bin_of_item[1] = 1;
    invalid.bin_count = 2;
    require(!precpack::check_assignment(instance, invalid, &diagnostic),
            "checker accepted a violated generalized precedence constraint");
}

void test_precedence_distance_range() {
    const int maximum_supported_separation =
        std::numeric_limits<int>::max() - 1;
    const precpack::Instance boundary = make_instance(
        {1, 1}, 1, {{0, 1, maximum_supported_separation}});
    require(boundary.front[1] == maximum_supported_separation &&
                boundary.back[0] == maximum_supported_separation,
            "maximum supported precedence distance changed");

    const int half_range = std::numeric_limits<int>::max() / 2 + 1;
    bool rejected = false;
    try {
        static_cast<void>(make_instance(
            {1, 1, 1}, 1,
            {{0, 1, half_range}, {1, 2, half_range}}));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected,
            "an out-of-range precedence path was not rejected");
}

void test_complete_dff_dual_feasibility() {
    const std::vector<int> weights{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    constexpr int capacity = 10;
    const precpack::DffTransformSet transforms =
        precpack::build_complete_dff_transforms(weights, capacity, true);
    require(transforms.item_count == static_cast<int>(weights.size()) &&
                transforms.size() > 100U,
            "complete DFF family was not generated");

    bool identity_found = false;
    for (std::size_t transform = 0;
         transform < transforms.size(); ++transform) {
        const std::int64_t* row = transforms.row(transform);
        identity_found = identity_found ||
            (transforms.capacities[transform] == capacity &&
             std::equal(weights.begin(), weights.end(), row));
        for (std::uint32_t mask = 0;
             mask < (1U << static_cast<unsigned>(weights.size())); ++mask) {
            int original_sum = 0;
            std::int64_t transformed_sum = 0;
            for (std::size_t item = 0; item < weights.size(); ++item) {
                if ((mask & (1U << static_cast<unsigned>(item))) != 0U) {
                    original_sum += weights[item];
                    transformed_sum += row[item];
                }
            }
            if (original_sum <= capacity) {
                require(
                    transformed_sum <= transforms.capacities[transform],
                    "generated transformation is not dual feasible");
            }
        }
    }
    require(identity_found, "complete DFF family lost the identity transform");
}

void test_initial_dff_exact_arithmetic() {
    require(precpack::compute_initial_dff_lower_bound(
                {670, 330}, 1000) == 1,
            "initial DFF rounded an exact one-bin boundary above one");
    require(precpack::compute_initial_dff_lower_bound(
                {335, 335, 330}, 1000) == 1,
            "initial DFF rounded an exact three-item fill above one");
    require(precpack::compute_initial_dff_lower_bound(
                {600, 600}, 1000) == 2,
            "initial DFF missed two items above half capacity");
    for (int capacity = 2; capacity <= 40; ++capacity) {
        for (int first = 1; first < capacity; ++first) {
            require(precpack::compute_initial_dff_lower_bound(
                        {first, capacity - first}, capacity) == 1,
                    "initial DFF violated dual feasibility");
        }
    }
}

void test_initialization_bounds() {
    const precpack::Instance instance = make_instance(
        {4, 3, 6, 2, 5}, 10,
        {{0, 2, 1}, {1, 3, 0}, {2, 4, 1}});
    const precpack::Assignment oracle = brute_force_optimum(instance);
    precpack::Config config = exact_config();
    config.bbr_enable_early_exact_probe = false;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    const precpack::InitialBoundsResult initial =
        precpack::compute_initial_bounds(
            instance, config, deadline, statistics);
    std::string diagnostic;
    require(initial.lower_bound <= oracle.bin_count &&
                initial.incumbent.bin_count >= oracle.bin_count &&
                precpack::check_assignment(
                    instance, initial.incumbent, &diagnostic),
            "initialization returned invalid bounds or incumbent: " +
                diagnostic);
}

void test_generalized_lower_bound_regression() {
    const precpack::Instance instance = make_instance(
        {4, 8, 8, 1, 11, 11, 1}, 11,
        {{0, 3, 3}, {1, 5, 3}, {2, 3, 2}, {3, 4, 0},
         {3, 5, 0}, {4, 6, 2}, {5, 6, 2}});
    const precpack::Assignment oracle = brute_force_optimum(instance);
    require(oracle.bin_count == 8,
            "generalized lower-bound regression oracle changed");

    precpack::Config config = exact_config();
    config.bbr_enable_early_exact_probe = false;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    const precpack::InitialBoundsResult initial =
        precpack::compute_initial_bounds(
            instance, config, deadline, statistics);
    require(initial.lower_bound <= oracle.bin_count,
            "generalized initialization produced an invalid lower bound");

    require_solver_matches_oracle(
        instance, precpack::ProblemKind::kBppGp);
}

void test_exact_state_memory() {
    const precpack::Instance instance = make_instance(
        {5, 3, 4, 1, 4, 3, 4, 4}, 10,
        {{0, 1, 3}, {1, 2, 1}, {1, 4, 3}, {3, 6, 3}, {5, 7, 3}});
    const int optimum = brute_force_optimum(instance).bin_count;
    const precpack::Assignment incumbent = spaced_assignment(instance);

    precpack::Config config = exact_config();
    config.bbr_enable_generalized_item_dominance = false;
    config.bbr_enable_profile_dominance = false;
    config.bbr_enable_paper_queue_order = false;
    config.bbr_enable_complete_dff = false;
    precpack::Deadline deadline(5.0);
    const precpack::BbrResult result =
        precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent),
            elementary_lower_bound(instance), config, deadline);
    require(result.optimal && result.incumbent.bin_count == optimum &&
                result.statistics.states_reopened > 0,
            "exact-state memory did not safely exercise reopen: optimal=" +
                std::to_string(result.optimal) +
                " ub=" + std::to_string(result.incumbent.bin_count) +
                " expected=" + std::to_string(optimum) +
                " reopened=" +
                std::to_string(result.statistics.states_reopened));

    const precpack::Instance collision_instance = make_instance(
        {6, 4, 6, 4, 6, 4, 6, 4, 6, 4, 6, 4}, 10);
    precpack::Config collision_config = exact_config();
    collision_config.bbr_enable_jackson = false;
    collision_config.bbr_enable_no_successor = false;
    collision_config.bbr_enable_superset_memory = false;
    collision_config.bbr_enable_profile_dominance = false;
    collision_config.bbr_enable_generalized_item_dominance = false;
    collision_config.bbr_enable_paper_queue_order = false;
    collision_config.bbr_enable_complete_dff = false;
    collision_config.bbr_enable_closure_bound = false;
    precpack::Deadline collision_deadline(5.0);
    const precpack::BbrResult collision_result =
        precpack::run_branch_bound_remember(
            make_identity_prepared(
                collision_instance, spaced_assignment(collision_instance)),
            6, collision_config, collision_deadline);
    require(collision_result.optimal &&
                collision_result.incumbent.bin_count == 6 &&
                collision_result.statistics.hash_probes >
                    collision_result.statistics.hash_lookups,
            "exact-state table collision probing was not exercised: "
            "lookups=" +
                std::to_string(collision_result.statistics.hash_lookups) +
                " probes=" +
                std::to_string(collision_result.statistics.hash_probes));
}

void test_profile_and_item_dominance() {
    const precpack::Instance profile_instance = make_instance(
        {2, 6, 3, 4, 6, 1, 3, 1}, 10,
        {{0, 5, 3}, {0, 7, 2}, {1, 2, 3}, {1, 5, 3}, {4, 5, 1}});
    const int profile_optimum =
        brute_force_optimum(profile_instance).bin_count;
    const precpack::Assignment profile_incumbent =
        spaced_assignment(profile_instance);
    precpack::Config profile_config = exact_config();
    profile_config.bbr_enable_jackson = false;
    profile_config.bbr_enable_generalized_item_dominance = false;
    profile_config.bbr_enable_paper_queue_order = false;
    profile_config.bbr_enable_complete_dff = false;
    precpack::Deadline profile_deadline(5.0);
    const precpack::BbrResult profile_result =
        precpack::run_branch_bound_remember(
            make_identity_prepared(profile_instance, profile_incumbent),
            elementary_lower_bound(profile_instance), profile_config,
            profile_deadline);
    require(profile_result.optimal &&
                profile_result.incumbent.bin_count == profile_optimum &&
                profile_result.statistics.profile_dominance_prunes > 0,
            "profile dominance was not exercised safely");

    profile_config.bbr_enable_profile_dominance = false;
    precpack::Deadline profile_baseline_deadline(5.0);
    const precpack::BbrResult profile_baseline =
        precpack::run_branch_bound_remember(
            make_identity_prepared(profile_instance, profile_incumbent),
            elementary_lower_bound(profile_instance), profile_config,
            profile_baseline_deadline);
    require(profile_baseline.optimal &&
                profile_baseline.incumbent.bin_count == profile_optimum,
            "profile-dominance baseline returned a different optimum");

    const precpack::Instance dominance_instance = make_instance(
        {3, 2, 4, 2}, 5, {{0, 1, 1}, {2, 3, 1}});
    const precpack::Assignment dominance_incumbent{{0, 2, 4, 6}, 7};
    precpack::Config baseline_config = exact_config();
    baseline_config.bbr_enable_jackson = false;
    baseline_config.bbr_enable_no_successor = false;
    baseline_config.bbr_enable_superset_memory = false;
    baseline_config.bbr_enable_paper_queue_order = false;
    baseline_config.bbr_enable_complete_dff = false;
    baseline_config.bbr_enable_generalized_item_dominance = false;
    precpack::Deadline baseline_deadline(5.0);
    const precpack::BbrResult baseline =
        precpack::run_branch_bound_remember(
            make_identity_prepared(
                dominance_instance, dominance_incumbent),
            2, baseline_config, baseline_deadline);

    baseline_config.bbr_enable_generalized_item_dominance = true;
    precpack::Deadline dominance_deadline(5.0);
    const precpack::BbrResult dominance =
        precpack::run_branch_bound_remember(
            make_identity_prepared(
                dominance_instance, dominance_incumbent),
            2, baseline_config, dominance_deadline);
    require(baseline.optimal && dominance.optimal &&
                baseline.incumbent.bin_count == 3 &&
                dominance.incumbent.bin_count == 3 &&
                dominance.statistics.generalized_item_dominance_pairs > 0 &&
                dominance.statistics.generalized_item_dominance_checks > 0 &&
                dominance.statistics.generalized_item_dominance_prunes > 0,
            "generalized item dominance was not exercised safely");

    std::mt19937 random(20260806U);
    std::uint64_t item_dominance_prunes = 0U;
    for (int trial = 0; trial < 20; ++trial) {
        constexpr int item_count = 8;
        std::vector<int> weights;
        weights.reserve(item_count);
        for (int item = 0; item < item_count; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < item_count; ++from) {
            for (int to = from + 1; to < item_count; ++to) {
                if (random() % 100U < 24U) {
                    arcs.push_back({from, to, 1});
                }
            }
        }
        const precpack::Instance instance = make_instance(
            std::move(weights), 10, std::move(arcs), "BPP-P");
        const int optimum = brute_force_optimum(instance).bin_count;
        const precpack::Assignment incumbent = spaced_assignment(instance, 2);
        precpack::Config config = exact_config(
            precpack::ProblemKind::kBppP);
        config.bbr_enable_generalized_item_dominance = false;
        config.bbr_enable_complete_dff = false;
        precpack::Deadline deadline(5.0);
        const precpack::BbrResult result =
            precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent),
                elementary_lower_bound(instance), config, deadline);
        require(result.optimal && result.incumbent.bin_count == optimum,
                "item dominance disagreed with brute force in trial " +
                    std::to_string(trial));
        item_dominance_prunes += result.statistics.jackson_prunes;
    }
    require(item_dominance_prunes > 0U,
            "labeled item dominance was not exercised");
}

void test_structured_preprocessing() {
    const precpack::Instance instance = make_instance(
        {2, 2, 10, 6, 6, 3, 6}, 10,
        {{0, 2, 1}, {2, 4, 1}, {1, 3, 1}, {3, 5, 1},
         {1, 6, 1}, {6, 5, 1}},
        "BPP-P");
    const int optimum = brute_force_optimum(instance).bin_count;

    precpack::Config config = exact_config(
        precpack::ProblemKind::kBppP);
    config.bbr_initialization_mode = true;
    config.bbr_enable_early_exact_probe = false;
    config.bbr_enable_initial_bdp = false;
    config.initialization_time_limit_seconds = 1.0;
    precpack::Deadline initialization_deadline(5.0);
    precpack::Statistics initialization_statistics;
    const precpack::InitialBoundsResult reduced =
        precpack::compute_initial_bounds(
            instance, config, initialization_deadline,
            initialization_statistics);
    require(reduced.prepared.structured_preprocessing_enabled &&
                !reduced.prepared.full_capacity_removals.empty() &&
                !reduced.prepared.fixed_prefix_bins.empty() &&
                !reduced.prepared.fixed_suffix_bins.empty() &&
                reduced.prepared.fixed_bin_offset >= 3 &&
                reduced.prepared.search_instance.size() < instance.size(),
            "structured BPP-P preprocessing did not exercise every reduction");

    const precpack::Assignment expanded_initial =
        precpack::map_prepared_assignment_to_original(
            reduced.prepared, reduced.prepared.search_incumbent);
    std::string diagnostic;
    require(precpack::check_assignment(
                instance, expanded_initial, &diagnostic),
            "structured initial assignment is invalid: " + diagnostic);
    precpack::Deadline reduced_deadline(5.0);
    const precpack::BbrResult reduced_result =
        precpack::run_branch_bound_remember(
            reduced.prepared, reduced.lower_bound, config,
            reduced_deadline);
    const precpack::Assignment reduced_assignment =
        precpack::map_prepared_assignment_to_original(
            reduced.prepared, reduced_result.incumbent);
    require(reduced_result.optimal &&
                reduced_result.certified_lower_bound == optimum &&
                reduced_assignment.bin_count == optimum &&
                precpack::check_assignment(
                    instance, reduced_assignment, &diagnostic),
            "structured BPP-P BBR disagreed with brute force: " + diagnostic);

    config.bbr_enable_structured_preprocessing = false;
    precpack::Deadline baseline_initialization_deadline(5.0);
    precpack::Statistics baseline_statistics;
    const precpack::InitialBoundsResult baseline =
        precpack::compute_initial_bounds(
            instance, config, baseline_initialization_deadline,
            baseline_statistics);
    require(!baseline.prepared.structured_preprocessing_enabled &&
                baseline.prepared.fixed_bin_offset == 0 &&
                baseline.prepared.search_instance.size() == instance.size(),
            "structured-preprocessing disable switch was ignored");
    precpack::Deadline baseline_deadline(5.0);
    const precpack::BbrResult baseline_result =
        precpack::run_branch_bound_remember(
            baseline.prepared, baseline.lower_bound, config,
            baseline_deadline);
    require(baseline_result.optimal &&
                baseline_result.certified_lower_bound == optimum &&
                baseline_result.incumbent.bin_count == optimum,
            "unreduced BPP-P baseline disagreed with brute force");
}

void test_bin_packing_bound() {
    const precpack::Instance instance =
        make_instance({6, 4, 6, 4, 5, 5}, 10);
    precpack::BinPackingBoundLimits limits;
    limits.call_time_limit_seconds = 2.0;
    limits.search_node_limit = 1'000'000U;
    limits.nondominated_load_limit_per_state = 100'000U;
    limits.memo_entry_limit = 100'000U;
    limits.maximum_item_count = 20;
    precpack::BinPackingBound bound(instance, limits);
    std::vector<std::uint64_t> remaining(
        bound.bit_block_count(), std::numeric_limits<std::uint64_t>::max());
    if (instance.size() % 64 != 0) {
        remaining.back() =
            (std::uint64_t{1} << static_cast<unsigned>(instance.size())) - 1U;
    }
    precpack::Deadline deadline(5.0);
    const precpack::BinPackingBoundResult result =
        bound.solve(remaining.data(), deadline);
    require(result.attempted && result.completed && result.optimum == 3,
            "ordinary BINLB returned the wrong optimum");
    int memoized = 0;
    require(bound.lookup_exact(remaining.data(), &memoized) &&
                memoized == 3,
            "ordinary BINLB did not retain its exact memo entry");
}

void test_small_exact_cases() {
    require_solver_matches_oracle(
        make_instance({4, 4}, 10, {{0, 1, 0}}, "SALBP-I"),
        precpack::ProblemKind::kSalbpI);
    require_solver_matches_oracle(
        make_instance({4, 4}, 10, {{0, 1, 1}}, "BPP-P"),
        precpack::ProblemKind::kBppP);
    for (int separation = 0; separation <= 3; ++separation) {
        require_solver_matches_oracle(
            make_instance({4, 4}, 10, {{0, 1, separation}}),
            precpack::ProblemKind::kBppGp);
    }
    require_solver_matches_oracle(
        make_instance({10, 4, 6}, 10,
                      {{0, 1, 1}, {1, 2, 1}}, "BPP-P"),
        precpack::ProblemKind::kBppP);
}

void test_random_bruteforce_oracle() {
    std::mt19937 random(20260821U);
    for (int case_index = 0; case_index < 36; ++case_index) {
        const int n = 4 + static_cast<int>(random() % 3U);
        const int capacity = 7 + static_cast<int>(random() % 6U);
        std::vector<int> weights(static_cast<std::size_t>(n));
        for (int& weight : weights) {
            weight = 1 + static_cast<int>(random() %
                static_cast<std::uint32_t>(capacity));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 24U) {
                    arcs.push_back({
                        from, to, static_cast<int>(random() % 3U)});
                }
            }
        }
        require_solver_matches_oracle(
            make_instance(std::move(weights), capacity, std::move(arcs)),
            precpack::ProblemKind::kBppGp);
    }
}

void test_resource_statuses() {
    const precpack::Instance instance =
        make_instance(std::vector<int>(12U, 1), 3);
    precpack::Assignment incumbent;
    incumbent.bin_of_item.resize(12U);
    std::iota(incumbent.bin_of_item.begin(),
              incumbent.bin_of_item.end(), 0);
    incumbent.bin_count = 12;
    precpack::PreparedInstance prepared;
    prepared.search_instance = instance;
    prepared.search_incumbent = incumbent;
    prepared.search_to_original.resize(12U);
    std::iota(prepared.search_to_original.begin(),
              prepared.search_to_original.end(), 0);
    prepared.original_to_search = prepared.search_to_original;
    prepared.original_item_count = instance.size();

    precpack::Config time_config = exact_config();
    precpack::Deadline time_deadline(1e-12);
    const precpack::BbrResult timed_out =
        precpack::run_branch_bound_remember(
            prepared, 4, time_config, time_deadline);
    require(!timed_out.optimal &&
                timed_out.timed_out,
            "time-limit termination was misclassified");
    require(std::string(precpack::to_string(precpack::SolveStatus::kOptimal)) ==
                    "OPTIMAL" &&
                std::string(precpack::to_string(
                    precpack::SolveStatus::kTimeLimit)) == "TIME_LIMIT" &&
                std::string(precpack::to_string(
                    precpack::SolveStatus::kMemoryLimit)) == "MEMORY_LIMIT",
            "public result statuses changed");
}

#if !PRECPACK_HAS_GUROBI
void test_gurobi_free_fallback() {
    const precpack::Instance instance =
        make_instance({6, 4, 6, 4}, 10);
    precpack::Config config = exact_config();
    config.bbr_enable_root_strengthening = true;
    const precpack::Solution solution = precpack::solve(instance, config);
    require(solution.optimal && !solution.root_stats.attempted,
            "Gurobi-free build did not fall back to exact BBR");

    config.require_gurobi_runtime = true;
    bool rejected = false;
    try {
        static_cast<void>(precpack::solve(instance, config));
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    require(rejected,
            "strict Gurobi runtime request was accepted by a Gurobi-free build");
}
#else
[[nodiscard]] GRBEnv make_environment() {
    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.set(GRB_IntParam_Threads, 1);
    environment.set(GRB_IntParam_Seed, 1);
    environment.start();
    return environment;
}

void test_gurobi_oracle_and_root_bounds() {
    GRBEnv environment = make_environment();
    const precpack::Instance packing =
        make_instance({6, 4, 6, 4}, 10);
    const precpack::Assignment packing_incumbent{{0, 1, 2, 3}, 4};
    precpack::Config config = exact_config();
    precpack::Deadline mip_deadline(10.0);
    const precpack::test::MipResult mip =
        precpack::test::solve_compact_mip(
            environment, packing, packing_incumbent, 2,
            config, mip_deadline);
    require(mip.optimal && mip.assignment.bin_count == 2 &&
                mip.certified_lower_bound == 2,
            "test-only compact MIP oracle returned the wrong optimum");

    const precpack::Instance generalized =
        make_instance({1, 1}, 10, {{0, 1, 2}});
    const precpack::Assignment incumbent{{0, 2}, 3};
    precpack::Deadline position_free_deadline(10.0);
    precpack::Statistics position_free_statistics;
    const precpack::RootStatistics position_free =
        precpack::run_position_free_root_column_generation(
            environment, generalized, incumbent, 1, config,
            position_free_deadline, position_free_statistics);
    require(position_free.completed &&
                position_free.certified_lower_bound == 2,
            "position-free root bound returned the wrong certificate");

    precpack::Deadline direct_deadline(10.0);
    precpack::Statistics direct_statistics;
    const precpack::BinIndexedRootBoundResult direct =
        precpack::run_bin_indexed_root_bound(
            environment, generalized, incumbent, 1, config,
            direct_deadline, direct_statistics);
    require(direct.completed && direct.certified_lower_bound == 3,
            "direct-precedence root bound returned the wrong certificate");
}
#endif

}

int main() {
    try {
        test_assignment_checker();
        test_precedence_distance_range();
        test_complete_dff_dual_feasibility();
        test_initial_dff_exact_arithmetic();
        test_initialization_bounds();
        test_generalized_lower_bound_regression();
        test_exact_state_memory();
        test_profile_and_item_dominance();
        test_structured_preprocessing();
        test_bin_packing_bound();
        test_small_exact_cases();
        test_random_bruteforce_oracle();
        test_resource_statuses();
#if PRECPACK_HAS_GUROBI
        test_gurobi_oracle_and_root_bounds();
#else
        test_gurobi_free_fallback();
#endif
        std::cout << "PrecPack exactness tests passed.\n";
        return 0;
#if PRECPACK_HAS_GUROBI
    } catch (const GRBException& error) {
        std::cerr << "Gurobi error " << error.getErrorCode() << ": "
                  << error.getMessage() << '\n';
        return 1;
#endif
    } catch (const std::exception& error) {
        std::cerr << "Exactness test failed: " << error.what() << '\n';
        return 1;
    }
}
