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

        const std::string instance_set = precpack::make_instance_set(
            options.instance_path, options.graph_path);
        const std::string instance_name =
            options.instance_path.stem().string();
        const std::filesystem::path assignment_path =
            options.output_directory / "solutions" /
            std::filesystem::path(instance_set) / (instance_name + ".sol");
        const std::filesystem::path result_path =
            options.output_directory /
            (std::string(precpack::to_string(options.problem)) +
             "_Results.csv");
        if (precpack::has_nonempty_solution(assignment_path)) {
            std::cout << "skip instance=" << instance_name
                      << " solution=" << assignment_path.string() << '\n';
            return 0;
        }
        const precpack::OutputLock output_lock(result_path);
        if (precpack::has_nonempty_solution(assignment_path)) {
            std::cout << "skip instance=" << instance_name
                      << " solution=" << assignment_path.string() << '\n';
            return 0;
        }
        precpack::validate_result_csv(result_path);

        const precpack::Config config =
            precpack::make_command_line_solver_config(
                options, options.time_limit_seconds);
        const precpack::Instance instance =
            precpack::read_instance(options.instance_path, options.graph_path,
                                    precpack::to_string(options.problem));
        const precpack::Solution solution = precpack::solve(instance, config);
        precpack::write_assignment(assignment_path, instance, solution);
        precpack::append_result_csv(
            result_path, instance_set, instance_name, instance, solution);

        std::cout << "status=" << precpack::to_string(solution.status)
                  << " problem=" << precpack::to_string(options.problem)
                  << " LB=" << solution.lower_bound
                  << " UB=" << solution.upper_bound
                  << " gap=" << std::setprecision(8) << solution.relative_gap
                  << " time=" << solution.stats.total_seconds << "s\n";
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
