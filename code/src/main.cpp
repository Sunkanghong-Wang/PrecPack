#include "precpack/build_config.hpp"
#include "precpack/batch.hpp"
#include "precpack/cli.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/output_lock.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#if PRECPACK_HAS_GUROBI
#include <gurobi_c++.h>
#endif

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        const precpack::CommandLineOptions options =
            precpack::parse_command_line(argc, argv);
        if (options.show_help) {
            precpack::print_help(std::cout, argv[0]);
            return 0;
        }
        if (options.batch_mode) {
            return precpack::run_batch(options);
        }

        const precpack::Config config = precpack::make_solver_config(
            options.problem, options.time_limit_seconds,
            options.memory_limit_mb, options.threads);
        const precpack::OutputLock output_lock(options.output_directory);
        const std::string key = precpack::make_instance_key(
            options.instance_path, options.graph_path);
        const std::filesystem::path assignment_path =
            options.output_directory / "solutions" /
            (std::string(precpack::to_slug(options.problem)) + "__" + key +
             ".sol");
        const std::filesystem::path result_path =
            options.output_directory /
            (std::string(precpack::to_string(options.problem)) +
             "_Results.csv");
        const std::filesystem::path assignment_reference =
            assignment_path.lexically_relative(options.output_directory);
        precpack::require_unused_instance_key(result_path, key);

        const precpack::Instance instance =
            precpack::read_instance(options.instance_path, options.graph_path,
                                    precpack::to_string(options.problem));
        const precpack::Solution solution = precpack::solve(instance, config);
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
