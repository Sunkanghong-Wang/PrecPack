#include "precpack/batch.hpp"

#include "precpack/build_config.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/output_lock.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#include "batch_schedule.hpp"
#include "environment.hpp"

#if PRECPACK_HAS_GUROBI
#include <gurobi_c++.h>
#endif

#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace precpack {
namespace {

constexpr std::array<const char*, 7> kOttoBaseDirectories = {
    "n_0020", "n_0050", "n_0100", "n_0250",
    "n_0500", "n_0750", "n_1000",
};

[[nodiscard]] bool require_gurobi_runtime() {
    const std::optional<std::string> value =
        internal::environment_value("PRECPACK_REQUIRE_GUROBI_RUNTIME");
    if (!value.has_value() || *value == "0") {
        return false;
    }
    if (*value != "1") {
        throw std::invalid_argument(
            "PRECPACK_REQUIRE_GUROBI_RUNTIME must be 0 or 1");
    }
    return true;
}

[[nodiscard]] std::filesystem::path repository_root() {
    const std::optional<std::string> value =
        internal::environment_value("PRECPACK_REPOSITORY_ROOT");
    if (!value.has_value()) {
        return std::filesystem::current_path();
    }
    return std::filesystem::absolute(*value).lexically_normal();
}

[[nodiscard]] std::filesystem::path caller_directory() {
    const std::optional<std::string> value =
        internal::environment_value("PRECPACK_CALLER_DIRECTORY");
    if (!value.has_value()) {
        return std::filesystem::current_path();
    }
    return std::filesystem::absolute(*value).lexically_normal();
}

[[nodiscard]] std::filesystem::path resolve_from(
    const std::filesystem::path& path,
    const std::filesystem::path& base) {
    if (path.is_absolute()) {
        return path.lexically_normal();
    }
    return (base / path).lexically_normal();
}

[[nodiscard]] std::optional<std::filesystem::path> resolve_from(
    const std::optional<std::filesystem::path>& path,
    const std::filesystem::path& base) {
    if (!path.has_value()) {
        return std::nullopt;
    }
    return resolve_from(*path, base);
}

[[nodiscard]] std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](char character) {
        return static_cast<char>(
            std::tolower(static_cast<unsigned char>(character)));
    });
    return value;
}

[[nodiscard]] std::filesystem::path normalized_existing_path(
    const std::filesystem::path& path) {
    return std::filesystem::weakly_canonical(path);
}

[[nodiscard]] std::vector<std::filesystem::path> files_below(
    const std::filesystem::path& path,
    std::string_view required_extension) {
    if (std::filesystem::is_regular_file(path)) {
        if (lowercase(path.extension().string()) !=
            lowercase(std::string(required_extension))) {
            throw std::invalid_argument(
                "expected a " + std::string(required_extension) +
                " file: " + path.string());
        }
        return {normalized_existing_path(path)};
    }
    if (!std::filesystem::is_directory(path)) {
        throw std::invalid_argument("input does not exist: " + path.string());
    }

    std::vector<std::filesystem::path> files;
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::recursive_directory_iterator(path)) {
        if (entry.is_regular_file() &&
            lowercase(entry.path().extension().string()) ==
                lowercase(std::string(required_extension))) {
            files.push_back(normalized_existing_path(entry.path()));
        }
    }
    std::sort(files.begin(), files.end(), [](const auto& left,
                                             const auto& right) {
        return left.generic_string() < right.generic_string();
    });
    return files;
}

[[nodiscard]] std::vector<BatchCase> collect_bpp_gp_cases(
    const std::filesystem::path& item_path,
    const std::filesystem::path& graph_path) {
    const std::vector<std::filesystem::path> instances =
        files_below(item_path, ".txt");
    const std::vector<std::filesystem::path> graphs =
        files_below(graph_path, ".graph");

    if (std::filesystem::is_regular_file(item_path)) {
        const std::filesystem::path& instance = instances.front();
        std::vector<BatchCase> cases;
        for (const std::filesystem::path& graph : graphs) {
            if (graph.stem() == instance.stem()) {
                cases.push_back(BatchCase{instance, graph});
            }
        }
        if (cases.empty()) {
            throw std::invalid_argument(
                "no graph with stem '" + instance.stem().string() +
                "' below " + graph_path.string());
        }
        return cases;
    }

    std::map<std::string, std::vector<std::filesystem::path>> by_stem;
    for (const std::filesystem::path& instance : instances) {
        by_stem[instance.stem().string()].push_back(instance);
    }

    const std::filesystem::path item_root =
        normalized_existing_path(item_path);
    std::vector<BatchCase> cases;
    cases.reserve(graphs.size());
    for (const std::filesystem::path& graph : graphs) {
        const auto found = by_stem.find(graph.stem().string());
        if (found == by_stem.end()) {
            throw std::invalid_argument(
                "no matching instance file for " + graph.string() +
                " below " + item_path.string());
        }
        const std::vector<std::filesystem::path>& candidates = found->second;
        const std::array<std::filesystem::path, 2> direct_candidates = {
            item_root / graph.parent_path().filename() /
                (graph.stem().string() + ".txt"),
            item_root / (graph.stem().string() + ".txt"),
        };

        std::optional<std::filesystem::path> instance;
        for (const std::filesystem::path& candidate : direct_candidates) {
            if (!std::filesystem::is_regular_file(candidate)) {
                continue;
            }
            const std::filesystem::path normalized =
                normalized_existing_path(candidate);
            if (std::find(candidates.begin(), candidates.end(), normalized) !=
                candidates.end()) {
                instance = normalized;
                break;
            }
        }
        if (!instance.has_value()) {
            std::vector<std::filesystem::path> same_size;
            for (const std::filesystem::path& candidate : candidates) {
                if (candidate.parent_path().filename() ==
                    graph.parent_path().filename()) {
                    same_size.push_back(candidate);
                }
            }
            if (same_size.size() == 1U) {
                instance = same_size.front();
            } else if (candidates.size() == 1U) {
                instance = candidates.front();
            } else {
                std::string matches;
                for (const std::filesystem::path& candidate : candidates) {
                    if (!matches.empty()) {
                        matches += ", ";
                    }
                    matches += candidate.string();
                }
                throw std::invalid_argument(
                    "ambiguous instance match for " + graph.string() +
                    ": " + matches);
            }
        }
        cases.push_back(BatchCase{*instance, graph});
    }
    return cases;
}

[[nodiscard]] std::filesystem::path solution_reference(
    const BatchCase& batch_case) {
    const std::string instance_set = make_instance_set(
        batch_case.instance_path, batch_case.graph_path);
    return std::filesystem::path("solutions") /
           std::filesystem::path(instance_set) /
           (batch_case.instance_path.stem().string() + ".sol");
}

[[nodiscard]] std::string single_line(std::string_view text) {
    std::string result(text);
    for (char& character : result) {
        if (character == '\n' || character == '\r') {
            character = ' ';
        }
    }
    return result;
}

void record_batch_error(const std::filesystem::path& output_directory,
                        ProblemKind problem,
                        const BatchCase& batch_case,
                        std::size_t index,
                        std::size_t total,
                        double time_limit_seconds,
                        std::string_view category,
                        std::string_view diagnostic) noexcept {
    const std::filesystem::path log_path = output_directory / "errors.log";
    try {
        std::filesystem::create_directories(output_directory);
        std::ofstream output(log_path, std::ios::app | std::ios::binary);
        if (!output) {
            std::cerr << "  warning: cannot open batch error log: "
                      << log_path.string() << '\n';
            return;
        }
        output << "index=" << index + 1U << '/' << total
               << " problem=" << to_string(problem)
               << " instance=" << batch_case.instance_path.generic_string();
        if (batch_case.graph_path.has_value()) {
            output << " graph=" << batch_case.graph_path->generic_string();
        }
        output << " time_limit_seconds=" << std::setprecision(12)
               << time_limit_seconds << " category=" << category
               << " diagnostic=" << single_line(diagnostic) << '\n';
        output.flush();
        if (!output) {
            std::cerr << "  warning: cannot write batch error log: "
                      << log_path.string() << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "  warning: cannot record batch error in "
                  << log_path.string() << ": " << error.what() << '\n';
    } catch (...) {
        std::cerr << "  warning: cannot record batch error in "
                  << log_path.string() << '\n';
    }
}

}

std::vector<BatchCase> collect_batch_cases(
    ProblemKind problem,
    const std::optional<std::filesystem::path>& input_path,
    const std::optional<std::filesystem::path>& graph_path,
    const std::filesystem::path& instance_root,
    const std::filesystem::path& graph_root) {
    if (problem != ProblemKind::kBppGp && graph_path.has_value()) {
        throw std::invalid_argument("--graph-dir is only valid for bpp-gp");
    }
    if (problem == ProblemKind::kSalbpI) {
        std::vector<std::filesystem::path> roots;
        if (input_path.has_value()) {
            roots.push_back(*input_path);
        } else {
            roots = {instance_root / "otto", instance_root / "scholl"};
        }
        std::vector<BatchCase> cases;
        for (const std::filesystem::path& root : roots) {
            for (const std::filesystem::path& instance :
                 files_below(root, ".txt")) {
                cases.push_back(BatchCase{instance, std::nullopt});
            }
        }
        return cases;
    }
    if (problem == ProblemKind::kBppP) {
        std::vector<std::filesystem::path> roots;
        if (input_path.has_value()) {
            roots.push_back(*input_path);
        } else {
            for (const char* directory : kOttoBaseDirectories) {
                roots.push_back(instance_root / "otto" / directory);
            }
            roots.push_back(instance_root / "scholl");
        }
        std::vector<BatchCase> cases;
        for (const std::filesystem::path& root : roots) {
            for (const std::filesystem::path& instance :
                 files_below(root, ".txt")) {
                cases.push_back(BatchCase{instance, std::nullopt});
            }
        }
        return cases;
    }
    return collect_bpp_gp_cases(input_path.value_or(instance_root / "otto"),
                                graph_path.value_or(graph_root));
}

int run_batch(const CommandLineOptions& options) {
    const std::filesystem::path root = repository_root();
    const std::filesystem::path caller = caller_directory();
    const std::optional<std::filesystem::path> input_path =
        resolve_from(options.input_path, caller);
    const std::optional<std::filesystem::path> graph_directory =
        resolve_from(options.graph_directory, caller);
    const bool strict_gurobi = require_gurobi_runtime();
    const std::filesystem::path output_directory =
        resolve_from(options.output_directory, caller);
    const std::filesystem::path result_path =
        output_directory /
        (std::string(to_string(options.problem)) + "_Results.csv");
    std::vector<BatchCase> cases = collect_batch_cases(
        options.problem, input_path, graph_directory,
        root / "data/instances", root / "data/bpp-gp-graphs");
    if (cases.empty()) {
        throw std::invalid_argument("no instances selected");
    }
    std::vector<double> time_limits;
    time_limits.reserve(cases.size());
    std::set<double> distinct_time_limits;
    for (const BatchCase& batch_case : cases) {
        const double limit = internal::effective_batch_time_limit(
            options, batch_case, root / "data/instances");
        time_limits.push_back(limit);
        distinct_time_limits.insert(limit);
    }

    std::set<std::string> output_names;
    std::vector<std::filesystem::path> solution_references;
    solution_references.reserve(cases.size());
    std::vector<bool> completed(cases.size(), false);
    std::size_t complete_count = 0U;
    for (std::size_t index = 0U; index < cases.size(); ++index) {
        const std::filesystem::path reference =
            solution_reference(cases[index]);
        if (!output_names.insert(reference.generic_string()).second) {
            throw std::invalid_argument(
                "selected files contain duplicate solution names");
        }
        solution_references.push_back(reference);
        if (has_nonempty_solution(output_directory / reference)) {
            completed[index] = true;
            ++complete_count;
        }
    }

    std::unique_ptr<OutputLock> output_lock;
    if (complete_count < cases.size()) {
        output_lock = std::make_unique<OutputLock>(result_path);
        complete_count = 0U;
        for (std::size_t index = 0U; index < cases.size(); ++index) {
            completed[index] = has_nonempty_solution(
                output_directory / solution_references[index]);
            complete_count += completed[index] ? 1U : 0U;
        }
        if (complete_count < cases.size()) {
            validate_result_csv(result_path);
            if (strict_gurobi) {
                verify_gurobi_runtime();
            }
        }
    }

    std::cout << "PrecPack batch\n"
              << "  problem   : " << to_string(options.problem) << '\n'
              << "  selected  : " << cases.size() << '\n'
              << "  completed : " << complete_count << '\n'
              << "  pending   : " << cases.size() - complete_count << '\n'
              << "  output    : " << output_directory.string() << '\n'
              << "  time limit: ";
    bool first_limit = true;
    for (const double limit : distinct_time_limits) {
        if (!first_limit) {
            std::cout << ',';
        }
        std::cout << std::setprecision(12) << limit;
        first_limit = false;
    }
    std::cout << "s\n";

    std::size_t failures = 0U;
    for (std::size_t index = 0U; index < cases.size(); ++index) {
        const BatchCase& batch_case = cases[index];
        const double time_limit_seconds = time_limits[index];
        const std::filesystem::path& reference = solution_references[index];
        if (completed[index]) {
            std::cout << '[' << index + 1U << '/' << cases.size()
                      << "] skip "
                      << batch_case.instance_path.filename().string() << '\n';
            continue;
        }
        std::cout << '[' << index + 1U << '/' << cases.size() << "] "
                  << batch_case.instance_path.filename().string()
                  << " (limit=" << std::setprecision(12)
                  << time_limit_seconds << "s)" << std::endl;
        try {
            Config config = make_command_line_solver_config(
                options, time_limit_seconds);
            config.require_gurobi_runtime = strict_gurobi;
            const Instance instance = read_instance(
                batch_case.instance_path, batch_case.graph_path,
                to_string(options.problem));
            const Solution solution = solve(instance, config);
            if (!solution.assignment.complete()) {
                throw std::runtime_error(
                    "solver returned no complete validated assignment");
            }
            write_assignment(output_directory / reference, instance, solution);
            append_result_csv(
                result_path,
                make_instance_set(batch_case.instance_path,
                                  batch_case.graph_path),
                batch_case.instance_path.stem().string(), instance, solution);
            completed[index] = true;
            std::cout << "  status=" << to_string(solution.status)
                      << " LB=" << solution.lower_bound
                      << " UB=" << solution.upper_bound
                      << " gap=" << std::setprecision(8)
                      << solution.relative_gap
                      << " time=" << solution.stats.total_seconds << "s\n";
#if PRECPACK_HAS_GUROBI
        } catch (const GRBException& error) {
            std::cerr << "  Gurobi error " << error.getErrorCode() << ": "
                      << error.getMessage() << '\n';
            record_batch_error(
                output_directory, options.problem, batch_case, index,
                cases.size(), time_limit_seconds, "gurobi_error",
                std::to_string(error.getErrorCode()) + ": " +
                    error.getMessage());
            if (strict_gurobi) {
                throw;
            }
            ++failures;
#endif
        } catch (const std::exception& error) {
            std::cerr << "  error: " << error.what() << '\n';
            record_batch_error(output_directory, options.problem, batch_case,
                               index, cases.size(), time_limit_seconds,
                               "error", error.what());
            ++failures;
        } catch (...) {
            std::cerr << "  error: non-standard exception\n";
            record_batch_error(output_directory, options.problem, batch_case,
                               index, cases.size(), time_limit_seconds,
                               "error", "non-standard exception");
            ++failures;
        }
    }
    return failures == 0U ? 0 : 1;
}

}
