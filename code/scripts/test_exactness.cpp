#include "precpack/algorithms.hpp"
#include "precpack/bbr.hpp"
#include "precpack/bin_packing_bound.hpp"
#include "precpack/build_config.hpp"
#include "precpack/conflict_bin_packing.hpp"
#include "precpack/dff.hpp"
#include "precpack/initial_bounds.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#if PRECPACK_HAS_GUROBI
#include "precpack/branch_price.hpp"
#include "precpack/column_generation.hpp"
#include "precpack/gurobi_solver.hpp"
#include "precpack/root_column_generation.hpp"

#include <gurobi_c++.h>
#endif

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

#if defined(__has_feature)
#if __has_feature(address_sanitizer)
constexpr bool kSanitizerBuild = true;
#else
constexpr bool kSanitizerBuild = false;
#endif
#else
constexpr bool kSanitizerBuild = false;
#endif

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

[[nodiscard]] precpack::Instance make_instance(std::vector<int> weights,
                                            int capacity,
                                            std::vector<precpack::Arc> arcs = {}) {
    precpack::Instance instance;
    instance.capacity = capacity;
    instance.problem_type = "BPP-GP";
    for (std::size_t i = 0; i < weights.size(); ++i) {
        instance.items.push_back(
            precpack::Item{static_cast<int>(i), weights[static_cast<std::size_t>(i)]});
    }
    instance.arcs = std::move(arcs);
    instance.initialize();
    return instance;
}

#if !PRECPACK_HAS_GUROBI
[[nodiscard]] std::filesystem::path repository_root() {
    return std::filesystem::path(__FILE__)
        .parent_path()
        .parent_path()
        .parent_path();
}

void require_exact_solution(const precpack::Instance& instance,
                            precpack::ProblemKind problem,
                            int expected_optimum) {
    const precpack::Config config =
        precpack::make_solver_config(problem, 10.0, 512, 1);
    const precpack::Solution solution = precpack::solve(instance, config);
    std::string diagnostic;
    require(solution.optimal &&
                solution.status == precpack::SolveStatus::kOptimal &&
                solution.lower_bound == solution.upper_bound &&
                solution.upper_bound == expected_optimum &&
                precpack::check_assignment(instance, solution.assignment,
                                           &diagnostic),
            "Gurobi-free BBR failed to return a validated optimum: " +
                diagnostic);
    require(!solution.root_stats.attempted &&
                solution.bbr_stats.root_cg_mode ==
                    precpack::BbrRootCgMode::kNone,
            "the Gurobi-free build attempted commercial root strengthening");
}

void test_no_gurobi_public_profiles() {
    require(!precpack::kHasGurobiSupport,
            "the Gurobi-free test was built with Gurobi enabled");
    const std::filesystem::path root = repository_root();
    const std::filesystem::path bpp_instance =
        root / "data/instances/otto/n_0020/instance_n=20_1.txt";

    const precpack::Config requested_bpp_gp = precpack::make_solver_config(
        precpack::ProblemKind::kBppGp, 10.0, 512, 1);
    require(requested_bpp_gp.bbr_root_cg_mode ==
                precpack::BbrRootCgMode::kPriceAndSwitch,
            "the public BPP-GP profile no longer requests optional root "
            "strengthening");

    require_exact_solution(
        precpack::read_instance(root / "data/instances/scholl/"
                                       "Jackson/Jackson_c7.txt",
                                std::nullopt, "SALBP-I"),
        precpack::ProblemKind::kSalbpI, 8);
    require_exact_solution(
        precpack::read_instance(bpp_instance, std::nullopt, "BPP-P"),
        precpack::ProblemKind::kBppP, 5);
    require_exact_solution(
        precpack::read_instance(
            bpp_instance,
            root / "data/bpp-gp-graphs/separation-01/n_0020/"
                   "instance_n=20_1.graph",
            "BPP-GP"),
        precpack::ProblemKind::kBppGp, 3);
}

void test_unavailable_internal_method_is_rejected() {
    precpack::Instance instance;
    instance.capacity = 2;
    instance.items = {{0, 1}, {1, 1}};
    instance.initialize();
    precpack::Config config;
    config.exact_method = precpack::ExactMethod::kMip;
    bool rejected = false;
    try {
        static_cast<void>(precpack::solve(instance, config));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected,
            "a Gurobi-dependent internal method was accepted without Gurobi");
}
#endif

void test_bounds_heuristic_and_checker() {
    const precpack::Instance instance =
        make_instance({1, 1}, 10, {{0, 1, 2}});
    const precpack::LowerBounds bounds = precpack::compute_simple_lower_bounds(instance);
    require(bounds.capacity == 1, "capacity lower bound mismatch");
    require(bounds.precedence == 3, "precedence lower bound mismatch");
    std::mt19937 random(1);
    const precpack::Assignment assignment =
        precpack::construct_initial_assignment(instance, bounds.combined, random, 2);
    require(assignment.bin_count == 3, "heuristic failed to attain obvious optimum");
    std::string diagnostic;
    require(precpack::check_assignment(instance, assignment, &diagnostic),
            "heuristic assignment invalid: " + diagnostic);

    precpack::Assignment bad = assignment;
    bad.bin_of_item[1] = 1;
    bad.bin_count = 2;
    require(!precpack::check_assignment(instance, bad, &diagnostic),
            "checker accepted a violated generalized precedence constraint");
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
    for (std::size_t transform = 0; transform < transforms.size(); ++transform) {
        const std::int64_t* row = transforms.row(transform);
        if (transforms.capacities[transform] == capacity &&
            std::equal(weights.begin(), weights.end(), row)) {
            identity_found = true;
        }
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
                require(transformed_sum <= transforms.capacities[transform],
                        "generated transformation is not dual feasible");
            }
        }
    }
    require(identity_found, "complete DFF family lost the identity transform");
}

void test_initial_lb4_dff_exact_arithmetic() {
    require(precpack::compute_initial_dff_lower_bound({670, 330}, 1000) == 1,
            "LB4 DFF rounded an exact one-bin boundary above one");
    require(precpack::compute_initial_dff_lower_bound({335, 335, 330}, 1000) == 1,
            "LB4 DFF rounded a three-item exact fill above one");
    require(precpack::compute_initial_dff_lower_bound({600, 600}, 1000) == 2,
            "LB4 DFF failed to detect two items above half capacity");

    for (int capacity = 2; capacity <= 100; ++capacity) {
        for (int first = 1; first < capacity; ++first) {
            const std::vector<int> weights{first, capacity - first};
            require(precpack::compute_initial_dff_lower_bound(weights, capacity) == 1,
                    "LB4 DFF violated dual feasibility for capacity=" +
                        std::to_string(capacity) + " split=" +
                        std::to_string(first));
        }
    }

    for (int capacity = 2; capacity <= 24; ++capacity) {
        std::vector<int> pattern;
        const std::function<void(int, int)> enumerate =
            [&](int remaining, int minimum_weight) {
                if (remaining == 0) {
                    require(
                        precpack::compute_initial_dff_lower_bound(pattern, capacity) == 1,
                        "LB4 DFF violated dual feasibility for an integer partition "
                        "of capacity=" + std::to_string(capacity));
                    return;
                }
                for (int weight = minimum_weight; weight <= remaining; ++weight) {
                    pattern.push_back(weight);
                    enumerate(remaining - weight, weight);
                    pattern.pop_back();
                }
            };
        enumerate(capacity, 1);
    }
}

#if PRECPACK_HAS_GUROBI
GRBEnv make_environment() {
    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.set(GRB_IntParam_Threads, 1);
    environment.set(GRB_IntParam_Seed, 1);
    environment.start();
    return environment;
}

void test_compact_mip() {
    const precpack::Instance instance = make_instance({6, 4, 6, 4}, 10);
    precpack::Assignment incumbent{{0, 1, 2, 3}, 4};
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.gurobi_log = false;
    precpack::Deadline deadline(10.0);
    GRBEnv environment = make_environment();
    const precpack::MipResult result = precpack::solve_compact_mip(
        environment, instance, incumbent, 2, config, deadline);
    require(result.optimal, "compact MIP did not prove a four-item test case");
    require(result.assignment.bin_count == 2 && result.certified_lower_bound == 2,
            "compact MIP returned the wrong optimum");
}

void test_initial_bpp_column_generation() {
    const precpack::Instance instance = make_instance({6, 4, 6, 4}, 10);
    const precpack::Assignment incumbent{{0, 1, 2, 3}, 4};
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.enable_sr_cuts = true;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const precpack::ColumnGenerationResult result =
        precpack::run_initial_bpp_column_generation(
            environment, instance, incumbent, 1, config, deadline, statistics);
    require(result.converged && result.pricing_proven,
            "initial BPP column generation did not converge");
    require(result.integer_lower_bound == 2,
            "initial BPP column generation returned the wrong bound");
    require(statistics.generated_columns > 0 &&
                statistics.pricing_search_nodes > 0,
            "initial pricing path was not exercised");
}

void test_initial_pipeline_reference_oracle() {
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const std::filesystem::path data = repository_root / "data/instances";
    const precpack::Instance instance = precpack::read_instance(
        data / "otto/n_0020/instance_n=20_6.txt",
        repository_root /
            "data/bpp-gp-graphs/separation-01/n_0020/instance_n=20_6.graph",
        "BPP-GP-01", 6);
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const auto provider = [&environment]() -> GRBEnv& { return environment; };
    const precpack::InitialBoundsResult result =
        precpack::compute_initial_bounds(
            instance, config, deadline, statistics, provider);
    require(result.lb1 == 3 && result.lb2 == 3 && result.lb3 == 3 &&
                result.lb4 == 3 && result.lower_bound == 3,
            "LB1-LB4 initialization disagreed with the fixed oracle");
    require(result.incumbent.bin_count == 4,
            "initial heuristics disagreed with the fixed oracle");
    require(result.initial_column_generation_attempted && statistics.cg_count == 1 &&
                statistics.generated_columns > 0,
            "fixed oracle did not exercise initial column/row generation");
    std::string diagnostic;
    require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
            "mapped initial incumbent is invalid: " + diagnostic);

    const precpack::Instance alns_instance = precpack::read_instance(
        data / "otto/n_0050/instance_n=50_12.txt",
        repository_root /
            "data/bpp-gp-graphs/separation-01/n_0050/instance_n=50_12.graph",
        "BPP-GP-01", 12);
    precpack::Deadline alns_deadline(10.0);
    precpack::Statistics alns_statistics;
    const precpack::InitialBoundsResult alns_result =
        precpack::compute_initial_bounds(
            alns_instance, config, alns_deadline, alns_statistics, provider);
    const std::vector<int> expected_bins{
        0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 4, 5, 3, 0, 1, 3,
        2, 2, 1, 2, 5, 1, 4, 1, 5, 4, 3, 1, 2, 4, 5, 2, 3, 5,
        5, 1, 3, 3, 4, 2, 2, 2, 4, 5, 3, 2, 4, 3, 3};
    require(alns_result.preprocessing_reversed &&
                alns_result.lower_bound == 6 &&
                alns_result.incumbent.bin_count == 6 &&
                alns_result.incumbent.bin_of_item == expected_bins,
            "reversal/ALNS initialization disagreed with the fixed oracle");
    require(precpack::check_assignment(alns_instance, alns_result.incumbent,
                                     &diagnostic),
            "mapped ALNS incumbent is invalid: " + diagnostic);

    precpack::Config no_alns_config = config;
    no_alns_config.enable_initial_alns = false;
    precpack::Deadline no_alns_deadline(10.0);
    precpack::Statistics no_alns_statistics;
    const precpack::InitialBoundsResult no_alns_result =
        precpack::compute_initial_bounds(
            alns_instance, no_alns_config, no_alns_deadline,
            no_alns_statistics, provider);
    require(no_alns_result.lower_bound == 6 &&
                no_alns_result.incumbent.bin_count == 7,
            "--no-initial-alns did not preserve the pre-ALNS incumbent");
    require(precpack::check_assignment(alns_instance, no_alns_result.incumbent,
                                     &diagnostic),
            "mapped pre-ALNS incumbent is invalid: " + diagnostic);
}
#endif

void test_bbr12_mhh_bounded_portfolio() {
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const precpack::Instance instance = precpack::read_instance(
        repository_root /
            "data/instances/scholl/Barthol2/Barthol2_c146.txt",
        std::nullopt, "SALBP-I", 56);
    const auto forbidden_environment = []() -> GRBEnv& {
        throw std::logic_error(
            "bounded BBR12 MHH regression unexpectedly requested Gurobi");
    };
    const auto run_initialization = [&](bool portfolio) {
        precpack::Config config;
        config.time_limit_seconds = 10.0;
        config.initialization_time_limit_seconds = 10.0;
        config.enable_initial_column_generation = false;
        config.enable_initial_alns = false;
        config.bbr_initialization_mode = true;
        config.bbr_enable_early_exact_probe = false;
        config.bbr_enable_initial_bdp = false;
        config.bbr_enable_initial_alns = false;
        config.bbr_enable_bbr12_mhh = true;
        config.bbr_enable_bbr12_mhh_portfolio = portfolio;
        config.bbr12_mhh_portfolio_max_items = 200;
        config.bbr12_mhh_full_load_limit = portfolio ? 1'000 : 50;
        precpack::Deadline deadline(10.0);
        precpack::Statistics statistics;
        return precpack::compute_initial_bounds(
            instance, config, deadline, statistics, forbidden_environment);
    };

    const precpack::InitialBoundsResult capped = run_initialization(false);
    const precpack::InitialBoundsResult portfolio = run_initialization(true);
    std::string diagnostic;
    require(capped.lower_bound == 29 && capped.incumbent.bin_count == 30,
            "50-load BBR12 MHH regression oracle changed");
    require(portfolio.lower_bound == 29 &&
                portfolio.incumbent.bin_count == 29,
            "bounded BBR12 MHH portfolio no longer repairs Barthol2/c146");
    require(precpack::check_assignment(instance, capped.incumbent, &diagnostic),
            "50-load BBR12 MHH incumbent is invalid: " + diagnostic);
    require(precpack::check_assignment(instance, portfolio.incumbent, &diagnostic),
            "bounded BBR12 MHH portfolio incumbent is invalid: " + diagnostic);
}

#if PRECPACK_HAS_GUROBI
void test_branch_price_tree() {
    const precpack::Instance instance = make_instance({6, 4, 6, 4}, 10);
    const precpack::Assignment incumbent{{0, 1, 2, 3}, 4};
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.enable_sr_cuts = false;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const precpack::BranchPriceResult result =
        precpack::run_branch_price_and_cut(environment, instance, incumbent, 1,
                                        config, deadline, statistics);
    require(result.optimal && result.incumbent.bin_count == 2 &&
                result.certified_lower_bound == 2,
            "branch-price-and-cut failed on the four-item test case");
    require(statistics.explored_nodes > 0 && statistics.generated_columns > 0 &&
                statistics.phase_one_count > 0,
            "branch-price test did not exercise the node RMP and pricing path");
}

void test_exact_bppc_bound() {
    const precpack::Instance instance = make_instance(
        {6, 6, 4, 4}, 10, {{0, 1, 0}, {1, 2, 0}, {1, 3, 0}});
    const precpack::Assignment original_incumbent{{0, 1, 1, 2}, 3};
    std::string diagnostic;
    require(precpack::check_assignment(instance, original_incumbent, &diagnostic),
            "BPPC test incumbent is invalid: " + diagnostic);

    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.cg_time_limit_seconds = 10.0;
    config.enable_sr_cuts = true;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const precpack::BppcBoundResult result =
        precpack::run_bppc_branch_price_bound(
            environment, instance, original_incumbent, 2, config, deadline,
            statistics);
    require(result.optimal && result.certified_lower_bound == 3 &&
                result.incumbent_value == 3,
            "exact BPPC tree failed to prove the closure-conflict bound: "
            "optimal=" + std::to_string(result.optimal) +
                " timed_out=" + std::to_string(result.timed_out) +
                " LB=" + std::to_string(result.certified_lower_bound) +
                " UB=" + std::to_string(result.incumbent_value) +
                " root_LB=" +
                std::to_string(result.root_integer_lower_bound) +
                " root_LP=" + std::to_string(result.root_lp_value));
    require(statistics.precedence_check_count == 0 &&
                statistics.generated_precedence_rows == 0 &&
                statistics.position_branches == 0,
            "BPPC tree used a forbidden precedence check/row/position branch");

    precpack::Config solver_config = config;
    solver_config.exact_method = precpack::ExactMethod::kBppc;
    const precpack::Solution solution = precpack::solve(instance, solver_config);
    require(solution.optimal && solution.lower_bound == 3 &&
                solution.upper_bound == 3 &&
                !solution.bppc_stats.attempted,
            "the public solve() BPPC mode did not preserve the exact "
            "initial-bound early exit");
    require(precpack::check_assignment(instance, solution.assignment, &diagnostic),
            "the public solve() BPPC path replaced the original feasible "
            "assignment: " + diagnostic);
}

void test_public_exact_bppc_mode() {
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const std::filesystem::path data = repository_root / "data/instances";
    const precpack::Instance instance = precpack::read_instance(
        data / "otto/n_0020/instance_n=20_6.txt",
        repository_root /
            "data/bpp-gp-graphs/separation-01/n_0020/instance_n=20_6.graph",
        "BPP-GP-01", 6);
    precpack::Config config;
    config.exact_method = precpack::ExactMethod::kBppc;
    config.time_limit_seconds = 10.0;
    config.cg_time_limit_seconds = 10.0;
    config.enable_initial_alns = false;
    const precpack::Solution solution = precpack::solve(instance, config);
    std::string diagnostic;
    require(!solution.optimal && solution.status == precpack::SolveStatus::kFeasible &&
                solution.lower_bound == 3 && solution.upper_bound == 4 &&
                solution.initial_upper_bound == 4 &&
                solution.bppc_stats.attempted &&
                solution.bppc_stats.completed &&
                solution.bppc_stats.certified_lower_bound == 3 &&
                solution.bppc_stats.incumbent_value == 3,
            "the public solve() BPPC path returned the wrong relaxation "
            "bound or promoted its incumbent to an original-problem upper "
            "bound");
    require(solution.stats.precedence_check_count == 0 &&
                solution.stats.generated_precedence_rows == 0 &&
                solution.stats.position_branches == 0,
            "the public solve() BPPC path used a forbidden precedence "
            "operation");
    require(precpack::check_assignment(instance, solution.assignment, &diagnostic),
            "the public solve() BPPC path returned an invalid original "
            "incumbent: " + diagnostic);
}

void test_bpc_timeout_is_not_a_proof() {
    const precpack::Instance instance = make_instance({6, 4, 6, 4}, 10);
    const precpack::Assignment incumbent{{0, 1, 2, 3}, 4};
    precpack::Config config;
    config.time_limit_seconds = 1.0;
    config.cg_time_limit_seconds = 1e-9;
    config.enable_sr_cuts = false;
    precpack::Deadline deadline(1.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const precpack::BranchPriceResult result =
        precpack::run_branch_price_and_cut(environment, instance, incumbent, 1,
                                        config, deadline, statistics);
    require(!result.optimal && result.timed_out &&
                result.certified_lower_bound == 1,
            "an interrupted BPC tree was incorrectly reported as a proof");
}

void test_root_only_stops_before_branching() {
    const precpack::Instance instance = make_instance({6, 4, 6, 4}, 10);
    const precpack::Assignment incumbent{{0, 1, 2, 3}, 4};
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.cg_time_limit_seconds = 10.0;
    config.root_node_only = true;
    config.enable_sr_cuts = false;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const precpack::BranchPriceResult result =
        precpack::run_branch_price_and_cut(environment, instance, incumbent, 1,
                                        config, deadline, statistics);
    require(result.root_only_completed && !result.timed_out &&
                result.root_integer_lower_bound == result.certified_lower_bound &&
                statistics.explored_nodes == 1 && statistics.rf_branches == 0 &&
                statistics.position_branches == 0,
            "root-only mode did not stop after the certified root node");
}

void test_three_set_covering_root_models() {
    const precpack::Instance instance =
        make_instance({1, 1}, 10, {{0, 1, 2}});
    const precpack::Assignment incumbent{{0, 2}, 3};
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.cg_time_limit_seconds = 10.0;
    config.enable_sr_cuts = false;
    config.max_columns_per_pricing = 4;
    GRBEnv environment = make_environment();

    precpack::Deadline fixed_deadline(10.0);
    precpack::Statistics fixed_statistics;
    const precpack::RootStatistics fixed =
        precpack::run_position_free_root_column_generation(
            environment, instance, incumbent, 1,
            precpack::RootModelKind::kFixedK, config, fixed_deadline,
            fixed_statistics);
    require(fixed.completed && !fixed.timed_out &&
                fixed.certified_lower_bound == 2 && fixed.iterations >= 2 &&
                fixed.pricing_count >= 2 &&
                fixed.fixed_point_scale_min > 0 &&
                fixed.fixed_point_scale_min <= fixed.fixed_point_scale_max &&
                fixed.fixed_point_scale_max <= (std::uint64_t{1} << 62U) &&
                fixed.certificate_scale_min >=
                    fixed.fixed_point_scale_max &&
                fixed.rmp_objective_multiplier >= 1.0 &&
                static_cast<long double>(fixed.fixed_point_scale_max) *
                        1.5e-9L <=
                    static_cast<long double>(
                        fixed.rmp_objective_multiplier),
            "fixed-K set-covering root returned the wrong certified bound");

    precpack::Deadline m_deadline(10.0);
    precpack::Statistics m_statistics;
    const precpack::RootStatistics m =
        precpack::run_position_free_root_column_generation(
            environment, instance, incumbent, 1,
            precpack::RootModelKind::kM, config, m_deadline, m_statistics);
    require(m.completed && !m.timed_out &&
                m.certified_lower_bound == 2 && m.pricing_count >= 2 &&
                m.fixed_point_scale_min > 0 &&
                m.fixed_point_scale_min <= m.fixed_point_scale_max &&
                m.fixed_point_scale_max <= (std::uint64_t{1} << 62U) &&
                m.certificate_scale_min >= m.fixed_point_scale_max &&
                m.rmp_objective_multiplier >= 1.0 &&
                static_cast<long double>(m.fixed_point_scale_max) *
                        1.5e-9L <=
                    static_cast<long double>(m.rmp_objective_multiplier),
            "M set-covering root returned the wrong certified bound");

    precpack::Config direct_config = config;
    direct_config.root_node_only = true;
    direct_config.set_covering_master = true;
    precpack::Deadline direct_deadline(10.0);
    precpack::Statistics direct_statistics;
    const precpack::BranchPriceResult direct =
        precpack::run_branch_price_and_cut(
            environment, instance, incumbent, 1, direct_config,
            direct_deadline, direct_statistics);
    require(direct.root_only_completed && !direct.timed_out &&
                direct.root_integer_lower_bound == 3 &&
                direct.root_column_count >= 2,
            "direct-precedence set-covering root returned the wrong bound");
}

void test_end_to_end_bounds_close() {
    const precpack::Instance instance =
        make_instance({1, 1}, 10, {{0, 1, 2}});
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.gurobi_log = false;
    const precpack::Solution solution = precpack::solve(instance, config);
    require(solution.optimal && solution.lower_bound == 3 && solution.upper_bound == 3,
            "end-to-end solver failed on a bound-closing test case");
}

void test_initial_only_stops_before_bpc() {
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const std::filesystem::path data = repository_root / "data/instances";
    const precpack::Instance instance = precpack::read_instance(
        data / "otto/n_0020/instance_n=20_6.txt",
        repository_root /
            "data/bpp-gp-graphs/separation-01/n_0020/instance_n=20_6.graph",
        "BPP-GP-01", 6);
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.initialization_only = true;
    const precpack::Solution solution = precpack::solve(instance, config);
    require(!solution.optimal && solution.status == precpack::SolveStatus::kFeasible &&
                solution.lower_bound == 3 && solution.upper_bound == 4,
            "initial-only mode returned the wrong fixed-oracle bounds/status");
    require(solution.root_lp_lower_bound == solution.initial_lower_bound &&
                solution.stats.explored_nodes == 0 &&
                solution.stats.rf_branches == 0 &&
                solution.stats.position_branches == 0,
            "initial-only mode entered the exact BPC/compact-MIP stage");
}

void test_m_root_sr_cut_and_tree() {
    const precpack::Instance instance = make_instance({4, 4, 4}, 8);
    const precpack::Assignment incumbent{{0, 1, 2}, 3};
    precpack::Config config;
    config.root_model = precpack::RootModelKind::kM;
    config.time_limit_seconds = 10.0;
    config.cg_time_limit_seconds = 10.0;
    config.enable_sr_cuts = true;
    config.sr_enumeration_item_limit = 100;
    precpack::Deadline deadline(10.0);
    precpack::Statistics statistics;
    GRBEnv environment = make_environment();
    const precpack::BranchPriceResult result =
        precpack::run_m_branch_price_and_cut(
            environment, instance, incumbent, 1, config, deadline,
            statistics);
    std::string diagnostic;
    require(result.optimal && result.certified_lower_bound == 2 &&
                result.incumbent.bin_count == 2,
            "M branch-price-and-cut failed the triplet SR oracle");
    require(statistics.generated_sr_rows > 0,
            "M root did not separate the violated triplet SR cut");
    require(statistics.generated_precedence_rows == 0 &&
                statistics.position_branches == 0,
            "M tree generated a forbidden precedence row/position branch");
    require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
            "M tree returned an invalid SR-oracle assignment: " + diagnostic);
}
#endif

[[nodiscard]] int brute_force_optimum(const precpack::Instance& instance);
[[nodiscard]] int brute_force_bppc_optimum(const precpack::Instance& instance);

[[nodiscard]] precpack::PreparedInstance make_identity_prepared(
    const precpack::Instance& instance,
    precpack::Assignment incumbent) {
    precpack::PreparedInstance prepared;
    prepared.search_instance = instance;
    prepared.search_incumbent = std::move(incumbent);
    prepared.search_to_original.resize(static_cast<std::size_t>(instance.size()));
    std::iota(prepared.search_to_original.begin(),
              prepared.search_to_original.end(), 0);
    prepared.original_to_search = prepared.search_to_original;
    return prepared;
}

void test_bbr_separation_transitions() {
    for (int separation = 0; separation <= 3; ++separation) {
        const precpack::Instance instance =
            make_instance({4, 4}, 10, {{0, 1, separation}});
        precpack::Assignment incumbent;
        incumbent.bin_of_item = {0, std::max(1, 2 * separation)};
        incumbent.bin_count = incumbent.bin_of_item[1] + 1;
        std::string diagnostic;
        require(precpack::check_assignment(instance, incumbent, &diagnostic),
                "BBR separation test incumbent is invalid: " + diagnostic);
        const int expected = separation == 0 ? 1 : separation + 1;
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        precpack::Config config;
        config.time_limit_seconds = 2.0;
        config.bbr_memory_limit_mb = 64;
        config.bbr_state_limit = 100'000;
        precpack::Deadline deadline(2.0);
        const precpack::BbrResult result = precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent), bounds.combined,
            config, deadline);
        require(result.optimal && result.certified_lower_bound == expected &&
                    result.incumbent.bin_count == expected,
                "BBR returned the wrong optimum for separation " +
                    std::to_string(separation));
        require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
                "BBR separation incumbent is invalid: " + diagnostic);
        if (separation >= 2) {
            require(result.statistics.forced_empty_transitions > 0,
                    "BBR did not exercise a forced cooldown empty bin");
        }
    }
}

void test_bbr_requested_configuration_metadata() {
    const auto run_case = [](std::string problem_type,
                             std::vector<precpack::Arc> arcs,
                             precpack::Assignment incumbent) {
        precpack::Instance instance = make_instance({6, 4}, 10, std::move(arcs));
        instance.problem_type = std::move(problem_type);
        precpack::Config config;
        config.time_limit_seconds = 2.0;
        config.bbr_memory_limit_mb = 64;
        config.bbr_state_limit = 100'000;
        config.bbr_enable_complete_dff = true;
        config.bbr_enable_generalized_item_dominance = true;
        config.bbr_enable_binlb = true;
        config.bbr_enable_conflict_binlb = true;
        for (const std::uint64_t heuristic_load_limit : {0U, 1U}) {
            config.bbr_heuristic_load_limit = heuristic_load_limit;
            precpack::Deadline deadline(2.0);
            const precpack::BbrResult result = precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent),
                incumbent.bin_count, config, deadline);
            require(
                result.optimal &&
                    result.certified_lower_bound == incumbent.bin_count &&
                    result.statistics.complete_dff_enabled &&
                    result.statistics.generalized_item_dominance_enabled &&
                    result.statistics.binlb_enabled &&
                    result.statistics.conflict_binlb_enabled &&
                    result.statistics.configured_state_limit ==
                        config.bbr_state_limit,
                "BBR early return lost requested configuration metadata for " +
                    instance.problem_type + " with heuristic load limit " +
                    std::to_string(heuristic_load_limit));
        }
    };

    run_case("SALBP-I", {{0, 1, 0}}, {{0, 0}, 1});
    run_case("BPP-P", {{0, 1, 1}}, {{0, 1}, 2});
    run_case("BPP-GP", {{0, 1, 1}}, {{0, 1}, 2});
}

void test_bbr_cooldown_profile_dominance() {
    const precpack::Instance instance = make_instance(
        {2, 6, 3, 4, 6, 1, 3, 1}, 10,
        {{0, 5, 3}, {0, 7, 2}, {1, 2, 3}, {1, 5, 3}, {4, 5, 1}});
    const int optimum = brute_force_optimum(instance);
    const precpack::LowerBounds bounds =
        precpack::compute_simple_lower_bounds(instance);
    std::mt19937 heuristic_random(29);
    precpack::Assignment incumbent = precpack::construct_initial_assignment(
        instance, bounds.combined, heuristic_random, 2);
    for (int& bin : incumbent.bin_of_item) {
        bin *= 2;
    }
    incumbent.bin_count =
        1 + *std::max_element(incumbent.bin_of_item.begin(),
                              incumbent.bin_of_item.end());
    std::string diagnostic;
    require(precpack::check_assignment(instance, incumbent, &diagnostic),
            "cooldown-profile test incumbent is invalid: " + diagnostic);

    precpack::Config config;
    config.time_limit_seconds = 2.0;
    config.bbr_memory_limit_mb = 64;
    config.bbr_state_limit = 100'000;
    config.bbr_enable_jackson = false;
    config.bbr_enable_paper_queue_order = false;
    config.bbr_enable_complete_dff = false;
    precpack::Deadline deadline(2.0);
    const precpack::BbrResult result = precpack::run_branch_bound_remember(
        make_identity_prepared(instance, incumbent), bounds.combined, config,
        deadline);
    require(result.optimal && result.incumbent.bin_count == optimum &&
                result.statistics.profile_dominance_prunes > 0,
            "BBR cooldown-profile dominance was not exercised safely: ub=" +
                std::to_string(result.incumbent.bin_count) +
                " states=" +
                std::to_string(result.statistics.states_created) +
                " profile=" +
                std::to_string(result.statistics.profile_dominance_prunes) +
                " bound=" + std::to_string(result.statistics.bound_prunes));

    config.bbr_enable_profile_dominance = false;
    precpack::Deadline baseline_deadline(2.0);
    const precpack::BbrResult baseline = precpack::run_branch_bound_remember(
        make_identity_prepared(instance, incumbent), bounds.combined, config,
        baseline_deadline);
    require(baseline.optimal && baseline.incumbent.bin_count == optimum,
            "BBR profile-disabled baseline returned a different optimum");
}

void test_bbr_exact_memory_reopen() {
    const precpack::Instance instance = make_instance(
        {5, 3, 4, 1, 4, 3, 4, 4}, 10,
        {{0, 1, 3}, {1, 2, 1}, {1, 4, 3}, {3, 6, 3}, {5, 7, 3}});
    const int optimum = brute_force_optimum(instance);
    const precpack::LowerBounds bounds =
        precpack::compute_simple_lower_bounds(instance);
    std::mt19937 heuristic_random(29);
    precpack::Assignment incumbent = precpack::construct_initial_assignment(
        instance, bounds.combined, heuristic_random, 2);
    for (int& bin : incumbent.bin_of_item) {
        bin *= 2;
    }
    incumbent.bin_count =
        1 + *std::max_element(incumbent.bin_of_item.begin(),
                              incumbent.bin_of_item.end());

    precpack::Config config;
    config.time_limit_seconds = 2.0;
    config.bbr_memory_limit_mb = 64;
    config.bbr_state_limit = 100'000;
    config.bbr_enable_profile_dominance = false;
    config.bbr_enable_paper_queue_order = false;
    config.bbr_enable_complete_dff = false;
    precpack::Deadline deadline(2.0);
    const precpack::BbrResult result = precpack::run_branch_bound_remember(
        make_identity_prepared(instance, incumbent), bounds.combined, config,
        deadline);
    require(result.optimal && result.incumbent.bin_count == optimum &&
                result.statistics.states_reopened > 0,
            "BBR did not safely reopen a better exact-state label");
}

void test_bbr_reverse_mapping_and_limits() {
    const precpack::Instance original =
        make_instance({4, 4, 6}, 10, {{0, 1, 2}});
    precpack::PreparedInstance reversed;
    reversed.search_instance =
        make_instance({4, 4, 6}, 10, {{1, 0, 2}});
    reversed.search_incumbent = {{4, 0, 4}, 5};
    reversed.search_to_original = {0, 1, 2};
    reversed.reversed = true;
    precpack::Config config;
    config.time_limit_seconds = 2.0;
    config.bbr_memory_limit_mb = 64;
    config.bbr_state_limit = 100'000;
    precpack::Deadline reverse_deadline(2.0);
    const precpack::BbrResult reverse_result =
        precpack::run_branch_bound_remember(reversed, 3, config,
                                         reverse_deadline);
    const precpack::Assignment mapped =
        precpack::map_prepared_assignment_to_original(reversed,
                                                    reverse_result.incumbent);
    std::string diagnostic;
    require(reverse_result.optimal && mapped.bin_count == 3 &&
                precpack::check_assignment(original, mapped, &diagnostic),
            "BBR reverse search/mapping failed: " + diagnostic);

    const precpack::Instance restart_instance = make_instance(
        {4, 4, 4, 1, 1, 1}, 10,
        {{0, 3, 1}, {1, 4, 1}, {2, 5, 1}});
    const precpack::PreparedInstance restart_prepared =
        make_identity_prepared(restart_instance, {{0, 2, 4, 6, 8, 10}, 11});
    precpack::Config two_phase_config = config;
    two_phase_config.bbr_heuristic_load_limit = 1;
    precpack::Deadline two_phase_deadline(2.0);
    const precpack::BbrResult two_phase_result =
        precpack::run_branch_bound_remember(restart_prepared, 1,
                                         two_phase_config,
                                         two_phase_deadline);
    require(two_phase_result.optimal &&
                two_phase_result.incumbent.bin_count == 3 &&
                two_phase_result.statistics.heuristic_phase_attempted &&
                two_phase_result.statistics.exact_phase_attempted &&
                two_phase_result.statistics.load_generation_truncated &&
                two_phase_result.statistics.truncated_load_states > 0 &&
                two_phase_result.statistics.exact_search_seconds > 0.0,
            "BBR did not restart a clean unrestricted proof phase");

    const precpack::Instance limit_instance = make_instance({6, 4, 6, 4}, 10);
    const precpack::PreparedInstance limit_prepared =
        make_identity_prepared(limit_instance, {{0, 2, 4, 6}, 7});
    precpack::Config memory_config = config;
    memory_config.bbr_state_limit = 1;
    precpack::Deadline memory_deadline(2.0);
    const precpack::BbrResult memory_result =
        precpack::run_branch_bound_remember(limit_prepared, 2, memory_config,
                                         memory_deadline);
    require(!memory_result.optimal && memory_result.state_limited &&
                !memory_result.memory_limited &&
                !memory_result.timed_out &&
                memory_result.certified_lower_bound <= 2,
            "BBR global state limit was incorrectly classified or certified");

    precpack::Deadline timeout_deadline(1e-12);
    const precpack::BbrResult timeout_result =
        precpack::run_branch_bound_remember(limit_prepared, 2, config,
                                         timeout_deadline);
    require(!timeout_result.optimal && timeout_result.timed_out &&
                !timeout_result.state_limited &&
                !timeout_result.memory_limited &&
                timeout_result.certified_lower_bound <= 2,
            "BBR timeout was incorrectly treated as an optimal proof");

    const precpack::Instance memory_instance = make_instance(
        std::vector<int>(16U, 1), 8);
    precpack::Assignment loose_incumbent;
    loose_incumbent.bin_of_item.resize(16U);
    for (int item = 0; item < 16; ++item) {
        loose_incumbent.bin_of_item[static_cast<std::size_t>(item)] = 2 * item;
    }
    loose_incumbent.bin_count = 31;
    precpack::Config byte_limit_config = config;
    byte_limit_config.bbr_memory_limit_mb = 1;
    byte_limit_config.bbr_state_limit = 60'000'000;
    byte_limit_config.bbr_heuristic_load_limit = 0;
    byte_limit_config.bbr_enable_jackson = false;
    byte_limit_config.bbr_enable_no_successor = false;
    byte_limit_config.bbr_enable_superset_memory = false;
    byte_limit_config.bbr_enable_profile_dominance = false;
    precpack::Deadline byte_limit_deadline(5.0);
    const precpack::BbrResult byte_limit_result =
        precpack::run_branch_bound_remember(
            make_identity_prepared(memory_instance, loose_incumbent), 2,
            byte_limit_config, byte_limit_deadline);
    require(!byte_limit_result.optimal && byte_limit_result.memory_limited &&
                !byte_limit_result.timed_out &&
                byte_limit_result.certified_lower_bound <= 2 &&
                byte_limit_result.statistics.peak_memory_bytes <=
                    1024U * 1024U,
            "BBR byte budget was exceeded or incorrectly certified: optimal=" +
                std::to_string(byte_limit_result.optimal) +
                " memory=" + std::to_string(byte_limit_result.memory_limited) +
                " timeout=" + std::to_string(byte_limit_result.timed_out) +
                " lb=" +
                std::to_string(byte_limit_result.certified_lower_bound) +
                " peak=" + std::to_string(
                    byte_limit_result.statistics.peak_memory_bytes) +
                " states=" + std::to_string(
                    byte_limit_result.statistics.states_created));
}

void test_bbr_bppp_structured_preprocessing() {
    precpack::Instance instance = make_instance(
        {2, 2, 10, 6, 6, 3, 6}, 10,
        {{0, 2, 1}, {2, 4, 1}, {1, 3, 1}, {3, 5, 1},
         {1, 6, 1}, {6, 5, 1}});
    instance.problem_type = "BPP-P";
    const int optimum = brute_force_optimum(instance);

    precpack::Config config;
    config.time_limit_seconds = 5.0;
    config.initialization_time_limit_seconds = 1.0;
    config.enable_initial_alns = false;
    config.enable_initial_column_generation = false;
    config.bbr_initialization_mode = true;
    config.bbr_enable_structured_preprocessing = true;
    config.bbr_memory_limit_mb = 64;
    config.bbr_state_limit = 100'000;

    const auto forbidden_environment = []() -> GRBEnv& {
        throw std::logic_error(
            "structured preprocessing unexpectedly requested Gurobi");
    };
    precpack::Deadline initialization_deadline(5.0);
    precpack::Statistics initialization_statistics;
    const precpack::InitialBoundsResult reduced =
        precpack::compute_initial_bounds(
            instance, config, initialization_deadline,
            initialization_statistics, forbidden_environment);
    require(reduced.prepared.structured_preprocessing_enabled &&
                !reduced.prepared.full_capacity_removals.empty() &&
                !reduced.prepared.fixed_prefix_bins.empty() &&
                !reduced.prepared.fixed_suffix_bins.empty() &&
                reduced.prepared.fixed_bin_offset >= 3 &&
                reduced.prepared.search_instance.size() < instance.size(),
            "2016 structured BPP-P preprocessing did not exercise every reduction");

    const precpack::Assignment expanded_initial =
        precpack::map_prepared_assignment_to_original(
            reduced.prepared, reduced.prepared.search_incumbent);
    std::string diagnostic;
    require(precpack::check_assignment(instance, expanded_initial, &diagnostic) &&
                expanded_initial.bin_count ==
                    reduced.prepared.search_incumbent.bin_count +
                        reduced.prepared.fixed_bin_offset,
            "structured BPP-P initial expansion is invalid: " + diagnostic);

    precpack::Deadline reduced_deadline(5.0);
    const precpack::BbrResult reduced_result =
        precpack::run_branch_bound_remember(
            reduced.prepared, reduced.lower_bound, config, reduced_deadline);
    const precpack::Assignment reduced_assignment =
        precpack::map_prepared_assignment_to_original(
            reduced.prepared, reduced_result.incumbent);
    require(reduced_result.optimal &&
                reduced_result.certified_lower_bound == optimum &&
                reduced_assignment.bin_count == optimum &&
                precpack::check_assignment(instance, reduced_assignment,
                                         &diagnostic),
            "structured BPP-P BBR disagreed with brute force: " + diagnostic);

    config.bbr_enable_structured_preprocessing = false;
    precpack::Deadline baseline_initialization_deadline(5.0);
    precpack::Statistics baseline_initialization_statistics;
    const precpack::InitialBoundsResult baseline =
        precpack::compute_initial_bounds(
            instance, config, baseline_initialization_deadline,
            baseline_initialization_statistics,
            forbidden_environment);
    require(!baseline.prepared.structured_preprocessing_enabled &&
                baseline.prepared.fixed_bin_offset == 0 &&
                baseline.prepared.search_instance.size() == instance.size(),
            "structured BPP-P preprocessing disable switch was ignored");
    precpack::Deadline baseline_deadline(5.0);
    const precpack::BbrResult baseline_result =
        precpack::run_branch_bound_remember(
            baseline.prepared, baseline.lower_bound, config,
            baseline_deadline);
    require(baseline_result.optimal &&
                baseline_result.certified_lower_bound == optimum &&
                baseline_result.incumbent.bin_count == optimum,
            "unreduced BPP-P baseline disagreed with the structured oracle");
}

void test_bbr_2016_generalized_item_dominance() {
    {
        const precpack::Instance instance = make_instance(
            {3, 2, 4, 2}, 5, {{0, 1, 1}, {2, 3, 1}});
        precpack::Assignment incumbent{{0, 2, 4, 6}, 7};
        std::string diagnostic;
        require(precpack::check_assignment(instance, incumbent, &diagnostic),
                "generalized-dominance incumbent is invalid: " + diagnostic);

        precpack::Config baseline;
        baseline.time_limit_seconds = 3.0;
        baseline.bbr_memory_limit_mb = 64;
        baseline.bbr_state_limit = 500'000;
        baseline.bbr_heuristic_load_limit = 0;
        baseline.bbr_enable_jackson = false;
        baseline.bbr_enable_no_successor = false;
        baseline.bbr_enable_superset_memory = false;
        baseline.bbr_enable_paper_queue_order = false;
        baseline.bbr_enable_complete_dff = false;
        precpack::Deadline baseline_deadline(3.0);
        const precpack::BbrResult baseline_result =
            precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent), 2, baseline,
                baseline_deadline);

        precpack::Config variant = baseline;
        variant.bbr_enable_generalized_item_dominance = true;
        precpack::Deadline variant_deadline(3.0);
        const precpack::BbrResult result = precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent), 2, variant,
            variant_deadline);
        require(baseline_result.optimal && result.optimal &&
                    baseline_result.incumbent.bin_count == 3 &&
                    result.incumbent.bin_count == 3 &&
                    result.statistics.generalized_item_dominance_pairs > 0 &&
                    result.statistics.generalized_item_dominance_checks > 0 &&
                    result.statistics.generalized_item_dominance_prunes > 0,
                "2016 generalized item dominance was not exercised safely");
    }

    std::mt19937 random(20260806);
    std::uint64_t basic_item_dominance_prunes = 0;
    for (int trial = 0; trial < 30; ++trial) {
        constexpr int n = 8;
        std::vector<int> weights;
        weights.reserve(n);
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 24U) {
                    arcs.push_back({from, to, 1});
                }
            }
        }
        const precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        const int optimum = brute_force_optimum(instance);
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        std::mt19937 heuristic_random(41);
        precpack::Assignment incumbent = precpack::construct_initial_assignment(
            instance, bounds.combined, heuristic_random, 2);
        for (int& bin : incumbent.bin_of_item) {
            bin *= 2;
        }
        incumbent.bin_count = 1 + *std::max_element(
            incumbent.bin_of_item.begin(), incumbent.bin_of_item.end());

        precpack::Config config;
        config.time_limit_seconds = 2.0;
        config.bbr_memory_limit_mb = 64;
        config.bbr_state_limit = 500'000;
        config.bbr_heuristic_load_limit = 0;
        config.bbr_enable_jackson = false;
        config.bbr_enable_no_successor = false;
        config.bbr_enable_superset_memory = false;
        config.bbr_enable_generalized_item_dominance = true;
        precpack::Deadline deadline(2.0);
        const precpack::BbrResult result = precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent), bounds.combined,
            config, deadline);
        require(result.optimal && result.incumbent.bin_count == optimum,
                "generalized BPP-P dominance disagreed with brute force in trial " +
                    std::to_string(trial));

        precpack::Config basic_config = config;
        basic_config.bbr_enable_jackson = true;
        basic_config.bbr_enable_generalized_item_dominance = false;
        precpack::Deadline basic_deadline(2.0);
        const precpack::BbrResult basic_result =
            precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent), bounds.combined,
                basic_config, basic_deadline);
        require(basic_result.optimal &&
                    basic_result.incumbent.bin_count == optimum,
                "basic Pereira BPP-P item dominance disagreed with brute force in trial " +
                    std::to_string(trial));
        basic_item_dominance_prunes +=
            basic_result.statistics.jackson_prunes;
    }
    require(basic_item_dominance_prunes > 0U,
            "basic Pereira BPP-P item dominance was not exercised");
}

void test_bbr_generalized_item_dominance_deadline() {
    const precpack::Instance instance = make_instance(
        {3, 2, 4, 2}, 5, {{0, 1, 1}, {2, 3, 1}});
    const precpack::Assignment incumbent{{0, 2, 4, 6}, 7};
    std::string diagnostic;
    require(precpack::check_assignment(instance, incumbent, &diagnostic),
            "deadline-test incumbent is invalid: " + diagnostic);

    precpack::Config config;
    config.time_limit_seconds = 0.0;
    config.bbr_memory_limit_mb = 64;
    config.bbr_state_limit = 500'000;
    config.bbr_heuristic_load_limit = 0;
    config.bbr_enable_complete_dff = false;
    config.bbr_enable_generalized_item_dominance = true;
    precpack::Deadline deadline(0.0);
    const precpack::BbrResult result = precpack::run_branch_bound_remember(
        make_identity_prepared(instance, incumbent), 2, config, deadline);

    require(result.timed_out && !result.optimal &&
                result.incumbent.bin_count == incumbent.bin_count &&
                result.statistics.generalized_item_dominance_search_nodes == 0,
            "expired deadline did not stop generalized-dominance preprocessing");
    require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
            "deadline-test result lost its valid incumbent: " + diagnostic);
}

#if PRECPACK_HAS_GUROBI
void test_bbr_root_dw_bound_integration() {
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const precpack::Instance instance = precpack::read_instance(
        repository_root /
            "data/instances/otto/n_0100/instance_n=100_54.txt",
        std::nullopt, "BPP-P", 54);
    precpack::Config config;
    config.time_limit_seconds = 5.0;
    config.initialization_time_limit_seconds = 2.0;
    config.bbr_enable_paper_queue_order = false;
    config.bbr_enable_complete_dff = false;
    config.bbr_root_cg_mode = precpack::BbrRootCgMode::kM;
    config.bbr_root_cg_time_limit_seconds = 1.0;
    const precpack::Solution solution = precpack::solve(instance, config);
    require(solution.optimal && solution.lower_bound == 51 &&
                solution.upper_bound == 51 &&
                solution.root_stats.attempted &&
                solution.root_stats.completed &&
                solution.root_stats.certified_lower_bound == 51 &&
                solution.bbr_stats.stop_reason == "ROOT_CG_OPTIMAL",
            "the integrated 2016 DW root bound did not safely close BPP-P 54");
}
#endif

void test_bbr_random_bruteforce_oracle() {
    std::mt19937 random(20260803);
    int searched_cases = 0;
    std::uint64_t profile_prunes = 0;
    std::uint64_t separation_item_prunes = 0;
    for (int trial = 0; trial < 40; ++trial) {
        constexpr int n = 8;
        std::vector<int> weights;
        weights.reserve(n);
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 24U) {
                    arcs.push_back(
                        {from, to, static_cast<int>(random() % 4U)});
                }
            }
        }
        precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        instance.id = trial + 1;
        const int optimum = brute_force_optimum(instance);
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        std::mt19937 heuristic_random(29);
        precpack::Assignment incumbent = precpack::construct_initial_assignment(
            instance, bounds.combined, heuristic_random, 2);
        for (int& bin : incumbent.bin_of_item) {
            bin *= 2;
        }
        incumbent.bin_count =
            1 + *std::max_element(incumbent.bin_of_item.begin(),
                                  incumbent.bin_of_item.end());

        precpack::Config config;
        config.time_limit_seconds = 3.0;
        config.bbr_memory_limit_mb = 64;
        config.bbr_state_limit = 500'000;
        precpack::Deadline deadline(3.0);
        const precpack::BbrResult result = precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent), bounds.combined,
            config, deadline);
        require(result.optimal && result.certified_lower_bound == optimum &&
                    result.incumbent.bin_count == optimum,
                "BBR disagreed with brute force in randomized trial " +
                    std::to_string(trial));
        std::string diagnostic;
        require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
                "randomized BBR incumbent is invalid: " + diagnostic);
        searched_cases += result.statistics.states_expanded > 0 ? 1 : 0;
        profile_prunes += result.statistics.profile_dominance_prunes;
        separation_item_prunes += result.statistics.jackson_prunes;

        if (trial < 12) {
            precpack::Config baseline = config;
            baseline.bbr_enable_jackson = false;
            baseline.bbr_enable_no_successor = false;
            baseline.bbr_enable_superset_memory = false;
            precpack::Deadline baseline_deadline(3.0);
            const precpack::BbrResult baseline_result =
                precpack::run_branch_bound_remember(
                    make_identity_prepared(instance, incumbent),
                    bounds.combined, baseline, baseline_deadline);
            require(baseline_result.optimal &&
                        baseline_result.incumbent.bin_count == optimum,
                    "BBR dominance-free baseline disagreed with brute force");
        }
    }
    require(searched_cases > 0,
            "randomized BBR oracle did not exercise the search engine");
    require(separation_item_prunes > 0U,
            "randomized BBR oracle did not exercise labeled-separation item dominance");
    (void)profile_prunes;
}

void test_bbr_salbp_dominance_switches() {
    std::mt19937 random(20260804);
    std::array<std::uint64_t, 3> exercised{};
    for (int trial = 0; trial < 100; ++trial) {
        constexpr int n = 10;
        std::vector<int> weights;
        weights.reserve(n);
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 27U) {
                    arcs.push_back({from, to, 0});
                }
            }
        }
        const precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        std::mt19937 heuristic_random(31);
        precpack::Assignment incumbent = precpack::construct_initial_assignment(
            instance, bounds.combined, heuristic_random, 2);
        for (int& bin : incumbent.bin_of_item) {
            bin *= 2;
        }
        incumbent.bin_count = 1 + *std::max_element(
            incumbent.bin_of_item.begin(), incumbent.bin_of_item.end());

        precpack::Config baseline;
        baseline.time_limit_seconds = 3.0;
        baseline.bbr_memory_limit_mb = 64;
        baseline.bbr_state_limit = 500'000;
        baseline.bbr_heuristic_load_limit = 0;
        baseline.bbr_enable_jackson = false;
        baseline.bbr_enable_no_successor = false;
        baseline.bbr_enable_superset_memory = false;
        baseline.bbr_enable_paper_queue_order = false;
        baseline.bbr_enable_complete_dff = false;
        precpack::Deadline baseline_deadline(3.0);
        const precpack::BbrResult baseline_result =
            precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent), bounds.combined,
                baseline, baseline_deadline);
        require(baseline_result.optimal,
                "dominance-free SALBP tree did not complete");
        const int optimum = baseline_result.incumbent.bin_count;

        for (int rule = 0; rule < 3; ++rule) {
            precpack::Config variant = baseline;
            variant.bbr_enable_jackson = rule == 0;
            variant.bbr_enable_no_successor = rule == 1;
            variant.bbr_enable_superset_memory = rule == 2;
            precpack::Deadline variant_deadline(3.0);
            const precpack::BbrResult result = precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent), bounds.combined,
                variant, variant_deadline);
            require(result.optimal && result.incumbent.bin_count == optimum,
                    "an isolated SALBP dominance switch deleted the optimum");
            exercised[static_cast<std::size_t>(rule)] +=
                rule == 0 ? result.statistics.jackson_prunes
                : rule == 1 ? result.statistics.no_successor_prunes
                            : result.statistics.superset_memory_prunes;
        }
    }
    require(std::all_of(exercised.begin(), exercised.end(),
                        [](std::uint64_t count) { return count > 0U; }),
            "isolated SALBP dominance tests did not exercise every rule: " +
                std::to_string(exercised[0]) + "," +
                std::to_string(exercised[1]) + "," +
                std::to_string(exercised[2]));
}

void test_bbr_bppp_generalized_maximum_load_dominance() {
    std::mt19937 random(20260814);
    std::uint64_t exercised = 0;
    for (int trial = 0; trial < 100; ++trial) {
        constexpr int n = 10;
        std::vector<int> weights;
        weights.reserve(n);
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 27U) {
                    arcs.push_back({from, to, 1});
                }
            }
        }
        const precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        std::mt19937 heuristic_random(43);
        precpack::Assignment incumbent = precpack::construct_initial_assignment(
            instance, bounds.combined, heuristic_random, 2);
        for (int& bin : incumbent.bin_of_item) {
            bin *= 2;
        }
        incumbent.bin_count = 1 + *std::max_element(
            incumbent.bin_of_item.begin(), incumbent.bin_of_item.end());

        precpack::Config baseline;
        baseline.time_limit_seconds = 3.0;
        baseline.bbr_memory_limit_mb = 64;
        baseline.bbr_state_limit = 500'000;
        baseline.bbr_heuristic_load_limit = 0;
        baseline.bbr_enable_jackson = false;
        baseline.bbr_enable_no_successor = false;
        baseline.bbr_enable_superset_memory = false;
        baseline.bbr_enable_paper_queue_order = false;
        baseline.bbr_enable_complete_dff = false;
        precpack::Deadline baseline_deadline(3.0);
        const precpack::BbrResult baseline_result =
            precpack::run_branch_bound_remember(
                make_identity_prepared(instance, incumbent), bounds.combined,
                baseline, baseline_deadline);
        require(baseline_result.optimal,
                "dominance-free BPP-P tree did not complete");

        precpack::Config variant = baseline;
        variant.bbr_enable_superset_memory = true;
        precpack::Deadline variant_deadline(3.0);
        const precpack::BbrResult result = precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent), bounds.combined,
            variant, variant_deadline);
        require(result.optimal &&
                    result.incumbent.bin_count ==
                        baseline_result.incumbent.bin_count,
                "2016 generalized maximum-load dominance deleted the BPP-P optimum");
        exercised += result.statistics.superset_memory_prunes;

    }
    require(exercised > 0U,
            "BPP-P generalized maximum-load dominance was not exercised");
}

#if PRECPACK_HAS_GUROBI
void test_bbr_random_mip_oracle() {
    std::mt19937 random(20260805);
    GRBEnv environment = make_environment();
    for (int trial = 0; trial < 9; ++trial) {
        const int n = 4 + trial;
        std::vector<int> weights;
        weights.reserve(static_cast<std::size_t>(n));
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 18U) {
                    arcs.push_back(
                        {from, to, static_cast<int>(random() % 4U)});
                }
            }
        }
        const precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        std::mt19937 heuristic_random(37);
        precpack::Assignment incumbent = precpack::construct_initial_assignment(
            instance, bounds.combined, heuristic_random, 3);
        for (int& bin : incumbent.bin_of_item) {
            bin *= 2;
        }
        incumbent.bin_count = 1 + *std::max_element(
            incumbent.bin_of_item.begin(), incumbent.bin_of_item.end());

        precpack::Config config;
        config.time_limit_seconds = 5.0;
        config.bbr_memory_limit_mb = 64;
        config.bbr_state_limit = 500'000;
        config.bbr_heuristic_load_limit = 0;
        precpack::Deadline bbr_deadline(5.0);
        const precpack::BbrResult bbr = precpack::run_branch_bound_remember(
            make_identity_prepared(instance, incumbent), bounds.combined,
            config, bbr_deadline);
        precpack::Deadline mip_deadline(5.0);
        const precpack::MipResult mip = precpack::solve_compact_mip(
            environment, instance, incumbent, bounds.combined, config,
            mip_deadline);
        std::string diagnostic;
        require(mip.optimal && bbr.optimal &&
                    bbr.incumbent.bin_count == mip.assignment.bin_count &&
                    bbr.certified_lower_bound == mip.certified_lower_bound,
                "BBR disagreed with compact MIP at random n=" +
                    std::to_string(n));
        require(precpack::check_assignment(instance, bbr.incumbent, &diagnostic),
                "BBR/MIP oracle produced an invalid BBR incumbent: " +
                    diagnostic);
    }
}
#endif

void test_bbr_salbp_464_without_dominance() {
    if (kSanitizerBuild) {
        return;
    }
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const precpack::Instance instance = precpack::read_instance(
        repository_root /
            "data/instances/otto/n_0100/instance_n=100_464.txt",
        std::nullopt, "SALBP-I", 464);
    precpack::Config initialization_config;
    initialization_config.time_limit_seconds = 10.0;
    initialization_config.initialization_time_limit_seconds = 3.0;
    initialization_config.enable_initial_column_generation = false;
    initialization_config.enable_initial_alns = false;
    initialization_config.bbr_initialization_mode = true;
    initialization_config.bbr_enable_early_exact_probe = false;
    precpack::Deadline initialization_deadline(10.0);
    precpack::Statistics statistics;
    const auto forbidden_environment = []() -> GRBEnv& {
        throw std::logic_error("BBR initialization unexpectedly requested Gurobi");
    };
    const precpack::InitialBoundsResult initial =
        precpack::compute_initial_bounds(
            instance, initialization_config, initialization_deadline,
            statistics, forbidden_environment);
    require(initial.lower_bound == 24 &&
                initial.prepared.search_incumbent.bin_count == 25,
            "SALBP-I 464 initialization oracle changed");

    precpack::Config probe_config = initialization_config;
    probe_config.bbr_enable_early_exact_probe = true;
    precpack::Deadline probe_deadline(10.0);
    precpack::Statistics probe_statistics;
    const precpack::InitialBoundsResult probed =
        precpack::compute_initial_bounds(
            instance, probe_config, probe_deadline, probe_statistics,
            forbidden_environment);
    require(probed.early_bbr_attempted && probed.early_bbr_optimal &&
                probed.lower_bound == 25 &&
                probed.incumbent.bin_count == 25 &&
                probed.early_bbr_statistics.exact_phase_attempted &&
                !probed.early_bbr_statistics.load_generation_truncated,
            "SALBP-I 464 early exact probe failed to close the known gap");
    std::string probe_diagnostic;
    require(precpack::check_assignment(instance, probed.incumbent,
                                     &probe_diagnostic),
            "SALBP-I 464 early probe incumbent is invalid: " +
                probe_diagnostic);

    precpack::Config exact_config = initialization_config;
    exact_config.bbr_memory_limit_mb = 256;
    exact_config.bbr_state_limit = 1'000'000;
    exact_config.bbr_heuristic_load_limit = 0;
    exact_config.bbr_enable_jackson = false;
    exact_config.bbr_enable_no_successor = false;
    exact_config.bbr_enable_superset_memory = false;
    exact_config.bbr_enable_profile_dominance = false;
    precpack::Deadline exact_deadline(10.0);
    const precpack::BbrResult exact = precpack::run_branch_bound_remember(
        initial.prepared, initial.lower_bound, exact_config, exact_deadline);
    require(exact.optimal && exact.certified_lower_bound == 25 &&
                exact.incumbent.bin_count == 25 &&
                exact.statistics.jackson_prunes == 0 &&
                exact.statistics.no_successor_prunes == 0 &&
                exact.statistics.superset_memory_prunes == 0,
            "dominance-free BBR did not independently prove SALBP-I 464 = 25");
}

void test_bbr_bppp_early_probe() {
    if (kSanitizerBuild) {
        return;
    }
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const precpack::Instance instance = precpack::read_instance(
        repository_root /
            "data/instances/otto/n_0100/instance_n=100_464.txt",
        std::nullopt, "BPP-P", 464);
    precpack::Config config;
    config.time_limit_seconds = 5.0;
    config.initialization_time_limit_seconds = 1.5;
    config.enable_initial_column_generation = false;
    config.enable_initial_alns = false;
    config.bbr_initialization_mode = true;
    config.bbr_enable_early_exact_probe = true;
    config.bbr_enable_generalized_item_dominance = true;
    precpack::Deadline deadline(5.0);
    precpack::Statistics statistics;
    const auto forbidden_environment = []() -> GRBEnv& {
        throw std::logic_error(
            "BPP-P early-probe test unexpectedly requested Gurobi");
    };
    const precpack::InitialBoundsResult result =
        precpack::compute_initial_bounds(
            instance, config, deadline, statistics, forbidden_environment);
    require(result.early_bbr_attempted && result.early_bbr_optimal &&
                result.lower_bound == 34 &&
                result.incumbent.bin_count == 34 &&
                result.early_bbr_statistics.exact_phase_attempted &&
                !result.early_bbr_statistics
                     .generalized_item_dominance_enabled &&
                !result.early_bbr_statistics.load_generation_truncated,
            "BPP-P 464 early exact probe failed to prove the known optimum");
    std::string diagnostic;
    require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
            "BPP-P 464 early-probe incumbent is invalid: " + diagnostic);

    config.bbr_root_cg_mode = precpack::BbrRootCgMode::kNone;
    const precpack::Solution end_to_end = precpack::solve(instance, config);
    require(end_to_end.optimal &&
                end_to_end.bbr_stats.stop_reason == "EARLY_BBR_OPTIMAL" &&
                end_to_end.bbr_stats.paper_queue_order_enabled &&
                end_to_end.bbr_stats.complete_dff_enabled &&
                end_to_end.bbr_stats.dff_transform_count >= 7U,
            "early BBR statistics lost the retained queue/DFF configuration");
}

void test_bbr_precpack_unified_startup() {
    const std::filesystem::path repository_root =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    const precpack::Instance instance = precpack::read_instance(
        repository_root /
            "data/instances/otto/n_0020/instance_n=20_6.txt",
        repository_root /
            "data/bpp-gp-graphs/separation-01/n_0020/instance_n=20_6.graph",
        "BPP-GP-01", 6);
    precpack::Config config;
    require(!config.bbr_enable_initial_alns &&
                config.bbr_enable_paper_queue_order &&
                config.bbr_enable_complete_dff &&
                !config.bbr_enable_binlb &&
                !config.bbr_enable_conflict_binlb &&
                config.bbr_root_cg_mode ==
                    precpack::BbrRootCgMode::kPriceAndSwitch &&
                config.bbr_root_cg_time_limit_seconds == 5.0,
            "the common BBR startup defaults changed");
    config.time_limit_seconds = 2.0;
    config.initialization_time_limit_seconds = 0.2;
    config.enable_initial_column_generation = false;
    config.enable_initial_alns = false;
    config.bbr_initialization_mode = true;
    precpack::Deadline deadline(2.0);
    precpack::Statistics statistics;
    const auto forbidden_environment = []() -> GRBEnv& {
        throw std::logic_error(
            "BPP-GP early-probe test unexpectedly requested Gurobi");
    };
    const precpack::InitialBoundsResult result =
        precpack::compute_initial_bounds(
            instance, config, deadline, statistics, forbidden_environment);
    require(result.early_bbr_attempted && result.early_bbr_optimal &&
                result.lower_bound == 3 && result.incumbent.bin_count == 3 &&
                result.early_bbr_statistics.exact_phase_attempted &&
                !result.early_bbr_statistics.load_generation_truncated,
            "BPP-GP-01 did not use the common exact startup probe");
    std::string diagnostic;
    require(precpack::check_assignment(instance, result.incumbent, &diagnostic),
            "BPP-GP-01 early-probe incumbent is invalid: " + diagnostic);
}

[[nodiscard]] int brute_force_optimum(const precpack::Instance& instance) {
    const precpack::LowerBounds bounds = precpack::compute_simple_lower_bounds(instance);
    std::mt19937 random(17);
    const precpack::Assignment initial =
        precpack::construct_initial_assignment(instance, bounds.combined, random, 2);
    int best = initial.bin_count;
    std::vector<int> bin_of_item(static_cast<std::size_t>(instance.size()), -1);
    std::vector<int> load(static_cast<std::size_t>(initial.bin_count), 0);

    std::function<void(int, int)> search = [&](int position, int last_bin) {
        if (position == instance.size()) {
            best = std::min(best, last_bin + 1);
            return;
        }
        if (last_bin + 1 >= best) {
            return;
        }
        const int item =
            instance.topological_order[static_cast<std::size_t>(position)];
        int earliest = 0;
        for (const auto& [predecessor, separation] :
             instance.predecessor_arcs[static_cast<std::size_t>(item)]) {
            earliest = std::max(
                earliest,
                bin_of_item[static_cast<std::size_t>(predecessor)] + separation);
        }
        for (int bin = earliest; bin + 1 < best; ++bin) {
            const int weight =
                instance.items[static_cast<std::size_t>(item)].weight;
            if (load[static_cast<std::size_t>(bin)] + weight > instance.capacity) {
                continue;
            }
            bin_of_item[static_cast<std::size_t>(item)] = bin;
            load[static_cast<std::size_t>(bin)] += weight;
            search(position + 1, std::max(last_bin, bin));
            load[static_cast<std::size_t>(bin)] -= weight;
            bin_of_item[static_cast<std::size_t>(item)] = -1;
        }
    };
    search(0, -1);
    return best;
}

[[nodiscard]] int brute_force_bppc_optimum(
    const precpack::Instance& instance) {
    const int n = instance.size();
    std::vector<unsigned char> reachable(
        static_cast<std::size_t>(n) * n, 0U);
    for (int item = 0; item < n; ++item) {
        reachable[static_cast<std::size_t>(item) * n + item] = 1U;
    }
    for (const precpack::Arc& arc : instance.arcs) {
        reachable[static_cast<std::size_t>(arc.from) * n + arc.to] = 1U;
    }
    for (int via = 0; via < n; ++via) {
        for (int from = 0; from < n; ++from) {
            if (reachable[static_cast<std::size_t>(from) * n + via] == 0U) {
                continue;
            }
            for (int to = 0; to < n; ++to) {
                if (reachable[static_cast<std::size_t>(via) * n + to] != 0U) {
                    reachable[static_cast<std::size_t>(from) * n + to] = 1U;
                }
            }
        }
    }
    std::vector<unsigned char> conflicts(
        static_cast<std::size_t>(n) * n, 0U);
    for (int lhs = 0; lhs < n; ++lhs) {
        for (int rhs = lhs + 1; rhs < n; ++rhs) {
            bool conflict = instance.separation(lhs, rhs) > 0 ||
                            instance.separation(rhs, lhs) > 0;
            if (!conflict) {
                int from = -1;
                int to = -1;
                if (reachable[static_cast<std::size_t>(lhs) * n + rhs] != 0U) {
                    from = lhs;
                    to = rhs;
                } else if (reachable[static_cast<std::size_t>(rhs) * n + lhs] !=
                           0U) {
                    from = rhs;
                    to = lhs;
                }
                if (from >= 0) {
                    int interval_load = 0;
                    for (int item = 0; item < n; ++item) {
                        if (reachable[static_cast<std::size_t>(from) * n + item] !=
                                0U &&
                            reachable[static_cast<std::size_t>(item) * n + to] !=
                                0U) {
                            interval_load += instance.items[
                                static_cast<std::size_t>(item)].weight;
                        }
                    }
                    conflict = interval_load > instance.capacity;
                }
            }
            if (conflict) {
                conflicts[static_cast<std::size_t>(lhs) * n + rhs] = 1U;
                conflicts[static_cast<std::size_t>(rhs) * n + lhs] = 1U;
            }
        }
    }

    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        const int lhs_degree = static_cast<int>(std::count(
            conflicts.begin() + static_cast<std::size_t>(lhs) * n,
            conflicts.begin() + static_cast<std::size_t>(lhs + 1) * n,
            static_cast<unsigned char>(1U)));
        const int rhs_degree = static_cast<int>(std::count(
            conflicts.begin() + static_cast<std::size_t>(rhs) * n,
            conflicts.begin() + static_cast<std::size_t>(rhs + 1) * n,
            static_cast<unsigned char>(1U)));
        if (lhs_degree != rhs_degree) {
            return lhs_degree > rhs_degree;
        }
        return instance.items[static_cast<std::size_t>(lhs)].weight >
               instance.items[static_cast<std::size_t>(rhs)].weight;
    });

    int best = n;
    std::vector<int> bin_of_item(static_cast<std::size_t>(n), -1);
    std::vector<int> load(static_cast<std::size_t>(n), 0);
    const std::function<void(int, int)> search =
        [&](int position, int used_bins) {
            if (used_bins >= best) {
                return;
            }
            if (position == n) {
                best = used_bins;
                return;
            }
            const int item = order[static_cast<std::size_t>(position)];
            for (int bin = 0; bin <= used_bins; ++bin) {
                if (bin == used_bins && used_bins + 1 >= best) {
                    break;
                }
                const int item_weight =
                    instance.items[static_cast<std::size_t>(item)].weight;
                if (load[static_cast<std::size_t>(bin)] + item_weight >
                    instance.capacity) {
                    continue;
                }
                bool compatible = true;
                for (int other = 0; other < n; ++other) {
                    if (bin_of_item[static_cast<std::size_t>(other)] == bin &&
                        conflicts[static_cast<std::size_t>(item) * n + other] !=
                            0U) {
                        compatible = false;
                        break;
                    }
                }
                if (!compatible) {
                    continue;
                }
                bin_of_item[static_cast<std::size_t>(item)] = bin;
                load[static_cast<std::size_t>(bin)] += item_weight;
                search(position + 1, std::max(used_bins, bin + 1));
                load[static_cast<std::size_t>(bin)] -= item_weight;
                bin_of_item[static_cast<std::size_t>(item)] = -1;
            }
        };
    search(0, 0);
    return best;
}

void test_ordinary_binlb() {
    std::mt19937 random(20260819);
    for (int trial = 0; trial < 32; ++trial) {
        constexpr int n = 9;
        std::vector<int> weights;
        weights.reserve(n);
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 9U));
        }
        const precpack::Instance instance =
            make_instance(std::move(weights), 10, {});
        const int expected = brute_force_bppc_optimum(instance);

        precpack::BinPackingBoundLimits limits;
        limits.call_time_limit_seconds = 2.0;
        limits.search_node_limit = 5'000'000U;
        limits.nondominated_load_limit_per_state = 1'000'000U;
        limits.memo_entry_limit = 100'000U;
        limits.maximum_item_count = n;
        precpack::BinPackingBound bound(instance, limits);
        std::vector<std::uint64_t> remaining(
            1U, (std::uint64_t{1} << n) - 1U);
        precpack::Deadline deadline(2.0);
        const precpack::BinPackingBoundResult result =
            bound.solve(remaining.data(), deadline);
        require(result.attempted && result.completed &&
                    result.optimum == expected,
                "ordinary BINLB disagreed with the brute-force BPP oracle in "
                "trial " + std::to_string(trial));
    }

    const precpack::Instance duplicate_instance =
        make_instance({6, 6, 4, 4, 3}, 10, {});
    precpack::BinPackingBoundLimits duplicate_limits;
    duplicate_limits.call_time_limit_seconds = 2.0;
    duplicate_limits.search_node_limit = 1'000'000U;
    duplicate_limits.nondominated_load_limit_per_state = 100'000U;
    duplicate_limits.memo_entry_limit = 10'000U;
    duplicate_limits.maximum_item_count = 10;
    precpack::BinPackingBound duplicate_bound(
        duplicate_instance, duplicate_limits);
    std::vector<std::uint64_t> first(
        1U, (std::uint64_t{1} << 0U) | (std::uint64_t{1} << 2U) |
                (std::uint64_t{1} << 4U));
    std::vector<std::uint64_t> equivalent(
        1U, (std::uint64_t{1} << 1U) | (std::uint64_t{1} << 3U) |
                (std::uint64_t{1} << 4U));
    precpack::Deadline duplicate_deadline(2.0);
    const precpack::BinPackingBoundResult first_result =
        duplicate_bound.solve(first.data(), duplicate_deadline);
    const precpack::BinPackingBoundResult equivalent_result =
        duplicate_bound.solve(equivalent.data(), duplicate_deadline);
    require(first_result.completed && equivalent_result.completed &&
                first_result.optimum == 2 &&
                equivalent_result.optimum == 2 &&
                equivalent_result.memo_hits > 0U,
            "ordinary BINLB did not reuse an equal-size multiset state");

    const precpack::Instance abort_instance =
        make_instance({7, 7, 7, 7, 4, 4, 4}, 10, {});
    precpack::BinPackingBoundLimits abort_limits;
    abort_limits.call_time_limit_seconds = 2.0;
    abort_limits.search_node_limit = 1U;
    abort_limits.nondominated_load_limit_per_state = 100U;
    abort_limits.memo_entry_limit = 100U;
    abort_limits.maximum_item_count = 20;
    precpack::BinPackingBound abort_bound(abort_instance, abort_limits);
    std::vector<std::uint64_t> all_items(
        1U, (std::uint64_t{1} << abort_instance.size()) - 1U);
    precpack::Deadline abort_deadline(0.0);
    const precpack::BinPackingBoundResult aborted =
        abort_bound.solve(all_items.data(), abort_deadline);
    require(aborted.attempted && !aborted.completed &&
                aborted.timed_out && aborted.optimum == 0,
            "an unfinished ordinary BINLB call exposed an integer optimum");

    const precpack::Instance barthol2_c85 = make_instance(
        {16, 30, 7, 47, 29, 8, 39, 37, 32, 29, 17, 11, 32, 15, 53,
         53, 8, 24, 24, 8, 7, 8, 14, 13, 10, 25, 11, 25, 11, 29, 25,
         10, 14, 41, 42, 47, 7, 80, 7, 41, 47, 16, 32, 66, 80, 7,
         41, 13, 47, 33, 34, 11, 18, 25, 7, 28, 12, 52, 14, 3, 3, 8,
         16, 33, 8, 18, 10, 14, 28, 11, 18, 25, 40, 40, 1, 5, 28, 8,
         81, 7, 26, 10, 21, 26, 20, 21, 47, 23, 13, 19, 15, 35, 26,
         46, 20, 31, 19, 34, 51, 39, 30, 26, 13, 45, 58, 28, 8, 83,
         40, 34, 23, 62, 11, 19, 14, 31, 32, 26, 55, 31, 32, 26, 19,
         14, 19, 48, 55, 8, 11, 27, 18, 36, 23, 20, 46, 64, 22, 15,
         34, 22, 51, 48, 64, 70, 37, 64, 78, 78},
        85, {});
    precpack::BinPackingBoundLimits c85_limits;
    c85_limits.call_time_limit_seconds = 0.02;
    c85_limits.search_node_limit = 5'000'000U;
    c85_limits.nondominated_load_limit_per_state = 1'000'000U;
    c85_limits.memo_entry_limit = 100'000U;
    c85_limits.maximum_item_count = barthol2_c85.size();
    precpack::BinPackingBound c85_bound(barthol2_c85, c85_limits);
    std::vector<std::uint64_t> c85_remaining(
        (static_cast<std::size_t>(barthol2_c85.size()) + 63U) / 64U,
        ~std::uint64_t{0});
    const unsigned c85_tail =
        static_cast<unsigned>(barthol2_c85.size()) & 63U;
    if (c85_tail != 0U) {
        c85_remaining.back() =
            (std::uint64_t{1} << c85_tail) - 1U;
    }
    precpack::Deadline c85_deadline(0.02);
    const precpack::BinPackingBoundResult c85_result =
        c85_bound.solve(c85_remaining.data(), c85_deadline);
    require(c85_result.lower_bound <= 50 &&
                (!c85_result.completed || c85_result.optimum == 50),
            "ordinary BINLB produced a false Barthol2/c85 bound");
}

void test_conflict_aware_binlb() {
    std::mt19937 random(20260818);
    for (int trial = 0; trial < 24; ++trial) {
        constexpr int n = 8;
        std::vector<int> weights;
        weights.reserve(n);
        for (int item = 0; item < n; ++item) {
            weights.push_back(1 + static_cast<int>(random() % 7U));
        }
        std::vector<precpack::Arc> arcs;
        for (int from = 0; from < n; ++from) {
            for (int to = from + 1; to < n; ++to) {
                if (random() % 100U < 28U) {
                    arcs.push_back(
                        {from, to, static_cast<int>(random() % 3U)});
                }
            }
        }
        const precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        const int expected = brute_force_bppc_optimum(instance);

        precpack::ConflictBinPackingLimits limits;
        limits.call_time_limit_seconds = 2.0;
        limits.total_time_limit_seconds = 2.0;
        limits.search_node_limit = 5'000'000U;
        limits.maximal_load_limit_per_state = 1'000'000U;
        limits.memo_entry_limit = 100'000U;
        limits.maximum_item_count = n;
        std::vector<std::uint64_t> remaining(1U,
            (std::uint64_t{1} << n) - 1U);
        precpack::Deadline deadline(2.0);

        precpack::BinPackingBoundLimits ordinary_limits;
        ordinary_limits.call_time_limit_seconds = 2.0;
        ordinary_limits.search_node_limit = 5'000'000U;
        ordinary_limits.nondominated_load_limit_per_state = 1'000'000U;
        ordinary_limits.memo_entry_limit = 100'000U;
        ordinary_limits.maximum_item_count = n;
        precpack::BinPackingBound ordinary(instance, ordinary_limits);
        const precpack::BinPackingBoundResult ordinary_result =
            ordinary.solve(remaining.data(), deadline);
        require(ordinary_result.completed,
                "ordinary BINLB did not complete before BPPC integration");
        precpack::ConflictBinPackingBound bound(instance, limits, &ordinary);
        const precpack::ConflictBinPackingResult result =
            bound.solve(remaining.data(), deadline);
        require(result.attempted && result.completed &&
                    result.optimum == expected &&
                    result.ordinary_memo_hits > 0U,
                "conflict-aware BINLB disagreed with the brute-force BPPC "
                "oracle in trial " + std::to_string(trial));

        const precpack::ConflictBinPackingResult cached =
            bound.solve(remaining.data(), deadline);
        require(cached.completed && cached.optimum == expected &&
                    cached.memo_hits > 0U,
                "conflict-aware BINLB did not reuse its collision-safe memo");

    }

    const precpack::Instance instance = make_instance(
        {7, 7, 7, 7, 4, 4, 4}, 10, {});
    precpack::ConflictBinPackingLimits abort_limits;
    abort_limits.call_time_limit_seconds = 2.0;
    abort_limits.total_time_limit_seconds = 2.0;
    abort_limits.search_node_limit = 1U;
    abort_limits.maximal_load_limit_per_state = 100U;
    abort_limits.memo_entry_limit = 100U;
    abort_limits.maximum_item_count = 20;
    precpack::ConflictBinPackingBound abort_bound(instance, abort_limits);
    std::vector<std::uint64_t> all_items(
        1U, (std::uint64_t{1} << instance.size()) - 1U);
    precpack::Deadline abort_deadline(2.0);
    const precpack::ConflictBinPackingResult aborted =
        abort_bound.solve(all_items.data(), abort_deadline);
    require(aborted.attempted && !aborted.completed &&
                aborted.node_limited && aborted.optimum == 0,
            "an unfinished conflict-aware BINLB call exposed an integer "
            "optimum");

    precpack::Assignment incumbent{{0, 2, 4, 6, 8, 8, 10}, 11};
    std::string diagnostic;
    require(precpack::check_assignment(instance, incumbent, &diagnostic),
            "conflict-aware BINLB BBR incumbent is invalid: " + diagnostic);
    precpack::Config config;
    config.time_limit_seconds = 2.0;
    config.bbr_memory_limit_mb = 64U;
    config.bbr_state_limit = 100'000U;
    config.bbr_enable_complete_dff = false;
    config.bbr_enable_conflict_binlb = true;
    config.bbr_conflict_binlb_call_time_limit_seconds = 1.0;
    config.bbr_conflict_binlb_total_time_limit_seconds = 1.0;
    config.bbr_conflict_binlb_node_limit = 1'000'000U;
    config.bbr_conflict_binlb_load_limit = 100'000U;
    config.bbr_conflict_binlb_max_items = 20;
    precpack::Deadline bbr_deadline(2.0);
    const precpack::BbrResult bbr = precpack::run_branch_bound_remember(
        make_identity_prepared(instance, incumbent), 4, config,
        bbr_deadline);
    require(bbr.optimal && bbr.incumbent.bin_count == 6 &&
                bbr.statistics.binlb_enabled &&
                bbr.statistics.binlb_calls > 0U &&
                bbr.statistics.binlb_completed > 0U &&
                bbr.statistics.conflict_binlb_enabled &&
                bbr.statistics.conflict_binlb_calls > 0U &&
                bbr.statistics.conflict_binlb_completed > 0U,
            "BBR did not safely exercise the conflict-aware BINLB: opt=" +
                std::to_string(bbr.optimal) + " ub=" +
                std::to_string(bbr.incumbent.bin_count) + " enabled=" +
                std::to_string(bbr.statistics.conflict_binlb_enabled) +
                " calls=" +
                std::to_string(bbr.statistics.conflict_binlb_calls) +
                " completed=" +
                std::to_string(bbr.statistics.conflict_binlb_completed) +
                " aborted=" +
                std::to_string(bbr.statistics.conflict_binlb_aborted) +
                " improvements=" + std::to_string(
                    bbr.statistics.conflict_binlb_bound_improvements));

    const precpack::Instance prune_instance =
        make_instance({4, 8, 10, 4, 3, 4, 8, 5}, 10, {});
    precpack::Assignment optimal_incumbent{{3, 1, 0, 4, 5, 4, 2, 3}, 6};
    require(precpack::check_assignment(prune_instance, optimal_incumbent,
                                    &diagnostic),
            "ordinary BINLB prune incumbent is invalid: " + diagnostic);
    precpack::Deadline prune_deadline(2.0);
    const precpack::BbrResult pruned = precpack::run_branch_bound_remember(
        make_identity_prepared(prune_instance, optimal_incumbent), 5, config,
        prune_deadline);
    require(pruned.optimal && pruned.statistics.binlb_calls == 1U &&
                pruned.statistics.binlb_completed == 1U &&
                pruned.statistics.binlb_prunes == 1U &&
                pruned.statistics.conflict_binlb_calls == 0U,
            "BBR did not skip BPPC after the ordinary BINLB pruned the root: "
            "bin_calls=" + std::to_string(pruned.statistics.binlb_calls) +
            " bin_prunes=" + std::to_string(pruned.statistics.binlb_prunes) +
            " conflict_calls=" +
            std::to_string(pruned.statistics.conflict_binlb_calls));

}

#if PRECPACK_HAS_GUROBI
void test_random_small_oracle() {
    std::mt19937 random(20260731);
    GRBEnv environment = make_environment();
    const auto environment_provider = [&environment]() -> GRBEnv& {
        return environment;
    };
    int solved_bpc_cases = 0;
    int solved_m_bpc_cases = 0;
    int solved_bppc_cases = 0;
    std::uint64_t bpc_rf_branches = 0;
    std::uint64_t bpc_position_branches = 0;
    std::uint64_t m_precedence_checks = 0;
    for (int trial = 0; trial < 12; ++trial) {
        constexpr int n = 6;
        std::vector<int> weights;
        weights.reserve(n);
        for (int i = 0; i < n; ++i) {
            weights.push_back(1 + static_cast<int>(random() % 6U));
        }
        std::vector<precpack::Arc> arcs;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (random() % 100U < 22U) {
                    arcs.push_back(
                        {i, j, static_cast<int>(random() % 3U)});
                }
            }
        }
        precpack::Instance instance =
            make_instance(std::move(weights), 10, std::move(arcs));
        instance.id = trial + 1;
        const int optimum = brute_force_optimum(instance);
        const precpack::LowerBounds bounds =
            precpack::compute_simple_lower_bounds(instance);
        require(bounds.combined <= optimum,
                "simple lower bound exceeded brute-force optimum");

        std::mt19937 heuristic_random(11);
        const precpack::Assignment incumbent = precpack::construct_initial_assignment(
            instance, bounds.combined, heuristic_random, 3);
        precpack::Config config;
        config.time_limit_seconds = 5.0;
        config.cg_time_limit_seconds = 2.0;
        config.enable_sr_cuts = false;
        std::string diagnostic;

        precpack::Deadline initial_deadline(2.0);
        precpack::Statistics initial_statistics;
        const precpack::InitialBoundsResult initial_bounds =
            precpack::compute_initial_bounds(
                instance, config, initial_deadline, initial_statistics,
                environment_provider);
        require(initial_bounds.lower_bound <= optimum,
                "initial bound exceeded brute-force optimum in trial " +
                    std::to_string(trial));
        require(precpack::check_assignment(instance, initial_bounds.incumbent,
                                         &diagnostic),
                "random initial incumbent is invalid: " + diagnostic);

        precpack::Deadline mip_deadline(5.0);
        const precpack::MipResult mip = precpack::solve_compact_mip(
            environment, instance, incumbent, bounds.combined, config,
            mip_deadline);
        require(mip.optimal && mip.assignment.bin_count == optimum,
                "compact MIP disagreed with the brute-force oracle in trial " +
                    std::to_string(trial));

        precpack::Assignment stretched = incumbent;
        int last_bin = -1;
        for (int& bin : stretched.bin_of_item) {
            bin *= 2;
            last_bin = std::max(last_bin, bin);
        }
        stretched.bin_count = last_bin + 1;
        require(precpack::check_assignment(instance, stretched, &diagnostic),
                "stretched BPC incumbent is invalid: " + diagnostic);
        if (trial < 6) {
            const int bppc_optimum = brute_force_bppc_optimum(instance);
            require(bounds.capacity <= bppc_optimum &&
                        bppc_optimum <= optimum,
                    "BPPC oracle violated the relaxation ordering");
            precpack::Config bppc_config = config;
            bppc_config.enable_sr_cuts = true;
            precpack::Deadline bppc_deadline(5.0);
            precpack::Statistics bppc_statistics;
            const precpack::BppcBoundResult bppc =
                precpack::run_bppc_branch_price_bound(
                    environment, instance, stretched, bounds.capacity,
                    bppc_config, bppc_deadline, bppc_statistics);
            require(bppc.optimal &&
                        bppc.certified_lower_bound == bppc_optimum &&
                        bppc.incumbent_value == bppc_optimum,
                    "BPPC tree disagreed with the brute-force oracle in trial " +
                        std::to_string(trial));
            require(bppc_statistics.precedence_check_count == 0 &&
                        bppc_statistics.generated_precedence_rows == 0 &&
                        bppc_statistics.position_branches == 0,
                    "BPPC random oracle used a precedence-specific operation");
            ++solved_bppc_cases;

            std::array<int, 2> position_free_bounds{};
            const std::array<precpack::RootModelKind, 2> root_kinds{
                precpack::RootModelKind::kFixedK,
                precpack::RootModelKind::kM};
            for (std::size_t model = 0; model < root_kinds.size(); ++model) {
                precpack::Deadline root_deadline(2.0);
                precpack::Statistics root_statistics;
                const precpack::RootStatistics root =
                    precpack::run_position_free_root_column_generation(
                        environment, instance, stretched, bounds.combined,
                        root_kinds[model], config, root_deadline,
                        root_statistics);
                require(root.completed && root.certified_lower_bound <= optimum,
                        "position-free root bound exceeded the brute-force oracle");
                position_free_bounds[model] = root.certified_lower_bound;
            }
            require(position_free_bounds[0] == position_free_bounds[1],
                    "completed fixed-K and M roots disagreed on integer LB");

            precpack::Config direct_root_config = config;
            direct_root_config.root_node_only = true;
            direct_root_config.set_covering_master = true;
            precpack::Deadline direct_root_deadline(2.0);
            precpack::Statistics direct_root_statistics;
            const precpack::BranchPriceResult direct_root =
                precpack::run_branch_price_and_cut(
                    environment, instance, stretched, bounds.combined,
                    direct_root_config, direct_root_deadline,
                    direct_root_statistics);
            require(direct_root.root_only_completed &&
                        direct_root.root_integer_lower_bound <= optimum,
                    "direct covering root bound exceeded the brute-force oracle");
        }
        precpack::Deadline bpc_deadline(5.0);
        precpack::Statistics bpc_statistics;
        const precpack::BranchPriceResult bpc =
            precpack::run_branch_price_and_cut(
                environment, instance, stretched, bounds.combined, config,
                bpc_deadline, bpc_statistics);
        require(bpc.optimal && bpc.incumbent.bin_count == optimum &&
                    bpc.certified_lower_bound == optimum,
                "BPC disagreed with the brute-force oracle in trial " +
                    std::to_string(trial));
        ++solved_bpc_cases;
        bpc_rf_branches += bpc_statistics.rf_branches;
        bpc_position_branches += bpc_statistics.position_branches;

        if (trial < 6) {
            precpack::Config m_config = config;
            m_config.root_model = precpack::RootModelKind::kM;
            m_config.enable_sr_cuts = true;
            m_config.sr_enumeration_item_limit = 100;
            precpack::Deadline m_deadline(5.0);
            precpack::Statistics m_statistics;
            const precpack::BranchPriceResult m_bpc =
                precpack::run_m_branch_price_and_cut(
                    environment, instance, stretched, bounds.combined,
                    m_config, m_deadline, m_statistics);
            require(m_bpc.optimal &&
                        m_bpc.incumbent.bin_count == optimum &&
                        m_bpc.certified_lower_bound == optimum,
                    "M BPC disagreed with the brute-force oracle in trial " +
                        std::to_string(trial));
            require(m_statistics.generated_precedence_rows == 0 &&
                        m_statistics.position_branches == 0,
                    "M BPC used a forbidden row or position branch");
            require(precpack::check_assignment(instance, m_bpc.incumbent,
                                             &diagnostic),
                    "M BPC random incumbent is invalid: " + diagnostic);
            ++solved_m_bpc_cases;
            m_precedence_checks += m_statistics.precedence_check_count;
        }
    }
    require(solved_bpc_cases > 0,
            "random oracle suite did not exercise the BPC solver");
    require(solved_m_bpc_cases == 6 && m_precedence_checks > 0,
            "random oracle suite did not exercise the M BPC check path");
    require(solved_bppc_cases == 6,
            "random oracle suite did not exercise the exact BPPC path");
    require(bpc_rf_branches > 0 && bpc_position_branches > 0,
            "random oracle suite did not exercise both BPC branch types");
    std::cout << "Random BPC oracle cases=" << solved_bpc_cases
              << " RF=" << bpc_rf_branches
              << " position=" << bpc_position_branches
              << " M-cases=" << solved_m_bpc_cases
              << " M-checks=" << m_precedence_checks
              << " BPPC-cases=" << solved_bppc_cases
              << '\n';
}
#endif

}

int main() {
    try {
#if !PRECPACK_HAS_GUROBI
        test_no_gurobi_public_profiles();
        test_unavailable_internal_method_is_rejected();
#endif
        test_bounds_heuristic_and_checker();
        test_complete_dff_dual_feasibility();
        test_initial_lb4_dff_exact_arithmetic();
        test_bbr12_mhh_bounded_portfolio();
        test_bbr_separation_transitions();
        test_bbr_requested_configuration_metadata();
        test_bbr_cooldown_profile_dominance();
        test_bbr_exact_memory_reopen();
        test_bbr_reverse_mapping_and_limits();
        test_bbr_bppp_structured_preprocessing();
        test_bbr_2016_generalized_item_dominance();
        test_bbr_generalized_item_dominance_deadline();
        test_bbr_random_bruteforce_oracle();
        test_bbr_salbp_dominance_switches();
        test_bbr_bppp_generalized_maximum_load_dominance();
        test_bbr_salbp_464_without_dominance();
        test_bbr_bppp_early_probe();
        test_bbr_precpack_unified_startup();
        test_ordinary_binlb();
        test_conflict_aware_binlb();

#if PRECPACK_HAS_GUROBI
        test_compact_mip();
        test_initial_bpp_column_generation();
        test_initial_pipeline_reference_oracle();
        test_branch_price_tree();
        test_exact_bppc_bound();
        test_public_exact_bppc_mode();
        test_bpc_timeout_is_not_a_proof();
        test_root_only_stops_before_branching();
        test_three_set_covering_root_models();
        test_end_to_end_bounds_close();
        test_initial_only_stops_before_bpc();
        test_m_root_sr_cut_and_tree();
        test_bbr_root_dw_bound_integration();
        test_bbr_random_mip_oracle();
        test_random_small_oracle();
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
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
}
