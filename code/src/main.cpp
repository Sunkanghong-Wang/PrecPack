#include "precpack/build_config.hpp"
#include "precpack/cli.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#if PRECPACK_HAS_GUROBI
#include <gurobi_c++.h>
#endif

#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

[[nodiscard]] std::string sanitize_filename(std::string value) {
    for (char& character : value) {
        const unsigned char code = static_cast<unsigned char>(character);
        if (std::isalnum(code) == 0 && character != '-' && character != '_' &&
            character != '.') {
            character = '_';
        }
    }
    return value.empty() ? "instance" : value;
}

[[nodiscard]] std::string instance_key(
    const precpack::CommandLineOptions& options) {
    std::string key = options.instance_path.stem().string();
    if (options.graph_path.has_value()) {
        std::filesystem::path collection =
            options.graph_path->parent_path().parent_path().filename();
        if (collection.empty()) {
            collection = options.graph_path->stem();
        }
        key += "__" + collection.string();
    }
    return sanitize_filename(std::move(key));
}

[[nodiscard]] std::string problem_slug(precpack::ProblemKind problem) {
    switch (problem) {
        case precpack::ProblemKind::kSalbpI:
            return "salbp-i";
        case precpack::ProblemKind::kBppP:
            return "bpp-p";
        case precpack::ProblemKind::kBppGp:
            return "bpp-gp";
    }
    throw std::logic_error("unknown problem kind");
}

}

int main(int argc, char** argv) {
    try {
        const precpack::CommandLineOptions options =
            precpack::parse_command_line(argc, argv);
        if (options.show_help) {
            precpack::print_help(std::cout, argv[0]);
            return 0;
        }

        const precpack::Config config = precpack::make_solver_config(
            options.problem, options.time_limit_seconds,
            options.memory_limit_mb, options.threads);
        const precpack::Instance instance =
            precpack::read_instance(options.instance_path, options.graph_path,
                                    precpack::to_string(options.problem));
        const precpack::Solution solution = precpack::solve(instance, config);

        const std::string key = instance_key(options);
        const std::filesystem::path assignment_path =
            options.output_directory / "solutions" /
            (problem_slug(options.problem) + "__" + key + ".sol");
        const std::filesystem::path result_path =
            options.output_directory /
            (std::string(precpack::to_string(options.problem)) +
             "_Results.csv");
        const std::filesystem::path assignment_reference =
            assignment_path.lexically_relative(options.output_directory);
        precpack::write_assignment(assignment_path, instance, solution);
        precpack::append_result_csv(
            result_path, key, options.instance_path, options.graph_path,
            instance, solution, assignment_reference);

        std::cout << "status=" << precpack::to_string(solution.status)
                  << " problem=" << precpack::to_string(options.problem)
                  << " LB=" << solution.lower_bound
                  << " UB=" << solution.upper_bound
                  << " gap=" << std::setprecision(8) << solution.relative_gap
                  << " time=" << solution.stats.total_seconds << "s"
                  << " threads=" << solution.threads << '\n';
        return solution.assignment.complete() ? 0 : 2;
#if PRECPACK_HAS_GUROBI
    } catch (const GRBException& error) {
        std::cerr << "Gurobi error " << error.getErrorCode() << ": "
                  << error.getMessage() << '\n';
        return 1;
#endif
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
