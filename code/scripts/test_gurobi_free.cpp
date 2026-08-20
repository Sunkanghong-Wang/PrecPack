#include "precpack/algorithms.hpp"
#include "precpack/build_config.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

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
    require(
        solution.optimal &&
            solution.status == precpack::SolveStatus::kOptimal &&
            solution.lower_bound == solution.upper_bound &&
            solution.upper_bound == expected_optimum &&
            precpack::check_assignment(instance, solution.assignment,
                                       &diagnostic),
        "Gurobi-free BBR failed to return a validated optimum: " + diagnostic);
    require(
        !solution.root_stats.attempted &&
            solution.bbr_stats.root_cg_mode == precpack::BbrRootCgMode::kNone,
        "the Gurobi-free build attempted commercial root strengthening");
}

void test_gurobi_free_exact_solver() {
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
        precpack::read_instance(root / "data/instances/scholl269/"
                                       "Jackson/Jackson_c7.txt",
                                std::nullopt, "SALBP-I"),
        precpack::ProblemKind::kSalbpI, 8);
    require_exact_solution(
        precpack::read_instance(bpp_instance, std::nullopt, "BPP-P"),
        precpack::ProblemKind::kBppP, 5);
    require_exact_solution(
        precpack::read_instance(
            bpp_instance,
            root / "data/bpp-gp-graphs/separation-01/n_0020/instance_n=20_1.graph",
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

}

int main() {
    try {
        test_gurobi_free_exact_solver();
        test_unavailable_internal_method_is_rejected();
        std::cout << "Gurobi-free exact-solver tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Gurobi-free test failed: " << error.what() << '\n';
        return 1;
    }
}
