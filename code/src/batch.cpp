#include "precpack/batch.hpp"

#include "precpack/build_config.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/output_lock.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver.hpp"
#include "precpack/solver_profile.hpp"

#if PRECPACK_HAS_GUROBI
#include <gurobi_c++.h>
#endif

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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
    const char* value = std::getenv("PRECPACK_REQUIRE_GUROBI_RUNTIME");
    if (value == nullptr || *value == '\0' || std::string_view(value) == "0") {
        return false;
    }
    if (std::string_view(value) != "1") {
        throw std::invalid_argument(
            "PRECPACK_REQUIRE_GUROBI_RUNTIME must be 0 or 1");
    }
    return true;
}

[[nodiscard]] std::filesystem::path repository_root() {
    const char* value = std::getenv("PRECPACK_REPOSITORY_ROOT");
    if (value == nullptr || *value == '\0') {
        return std::filesystem::current_path();
    }
    return std::filesystem::absolute(value).lexically_normal();
}

[[nodiscard]] std::filesystem::path caller_directory() {
    const char* value = std::getenv("PRECPACK_CALLER_DIRECTORY");
    if (value == nullptr || *value == '\0') {
        return std::filesystem::current_path();
    }
    return std::filesystem::absolute(value).lexically_normal();
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
    std::transform(value.begin(), value.end(), value.begin(), [](char value) {
        return static_cast<char>(
            std::tolower(static_cast<unsigned char>(value)));
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
    ProblemKind problem,
    const BatchCase& batch_case) {
    const std::string key =
        make_instance_key(batch_case.instance_path, batch_case.graph_path);
    return std::filesystem::path("solutions") /
           (std::string(to_slug(problem)) + "__" + key + ".sol");
}

[[nodiscard]] bool nonempty_regular_file(
    const std::filesystem::path& path) {
    std::error_code error;
    const bool regular = std::filesystem::is_regular_file(path, error);
    if (error || !regular) {
        return false;
    }
    const std::uintmax_t size = std::filesystem::file_size(path, error);
    return !error && size > 0U;
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
    std::vector<BatchCase> cases = collect_batch_cases(
        options.problem, input_path, graph_directory,
        root / "data/instances", root / "data/bpp-gp-graphs");
    if (cases.empty()) {
        throw std::invalid_argument("no instances selected");
    }

    const bool strict_gurobi = require_gurobi_runtime();
    if (strict_gurobi && !options.check_only) {
        verify_gurobi_runtime();
    }
    const std::filesystem::path output_directory =
        resolve_from(options.output_directory, caller);
    std::unique_ptr<OutputLock> output_lock;
    if (!options.check_only || std::filesystem::exists(output_directory)) {
        output_lock = std::make_unique<OutputLock>(output_directory);
    }
    const std::filesystem::path result_path =
        output_directory /
        (std::string(to_string(options.problem)) + "_Results.csv");
    Config config = make_solver_config(
        options.problem, options.time_limit_seconds, options.memory_limit_mb,
        options.threads);
    config.require_gurobi_runtime = strict_gurobi;
    const std::vector<ResultReference> references =
        read_result_references(result_path);
    std::map<std::string, ResultReference> recorded_results;
    for (const ResultReference& reference : references) {
        if (!recorded_results.emplace(reference.instance_key, reference).second) {
            throw std::runtime_error(
                "duplicate instance_key in result CSV: " +
                reference.instance_key);
        }
    }

    std::set<std::string> output_names;
    std::set<std::string> completed;
    std::size_t complete_count = 0U;
    for (const BatchCase& batch_case : cases) {
        const std::string key =
            make_instance_key(batch_case.instance_path, batch_case.graph_path);
        const std::filesystem::path reference =
            solution_reference(options.problem, batch_case);
        if (!output_names.insert(reference.generic_string()).second) {
            throw std::invalid_argument(
                "selected files contain duplicate output names");
        }
        const auto recorded = recorded_results.find(key);
        if (recorded != recorded_results.end()) {
            const ResultReference& result = recorded->second;
            if (result.problem != to_string(options.problem) ||
                std::abs(result.time_limit_seconds -
                         options.time_limit_seconds) >
                    1e-9 * std::max(1.0, options.time_limit_seconds) ||
                result.threads != resolve_thread_count(options.threads) ||
                result.state_limit != config.bbr_state_limit ||
                result.memory_limit_mb != options.memory_limit_mb ||
                result.gurobi_enabled != kHasGurobiSupport ||
                result.gurobi_required != strict_gurobi ||
                result.solution_file.generic_string() !=
                    reference.generic_string()) {
                throw std::runtime_error(
                    "existing result profile does not match the requested "
                    "batch for instance_key=" + key +
                    "; use a different output directory");
            }
            if (!nonempty_regular_file(output_directory / reference)) {
                throw std::runtime_error(
                    "result row has no matching nonempty solution for "
                    "instance_key=" + key);
            }
            completed.insert(key);
            ++complete_count;
        }
    }

    std::cout << "PrecPack batch\n"
              << "  problem   : " << to_string(options.problem) << '\n'
              << "  selected  : " << cases.size() << '\n'
              << "  completed : " << complete_count << '\n'
              << "  pending   : " << cases.size() - complete_count << '\n'
              << "  output    : " << output_directory.string() << '\n';
    if (options.check_only) {
        return 0;
    }

    std::size_t failures = 0U;
    for (std::size_t index = 0U; index < cases.size(); ++index) {
        const BatchCase& batch_case = cases[index];
        const std::string key =
            make_instance_key(batch_case.instance_path, batch_case.graph_path);
        const std::filesystem::path reference =
            solution_reference(options.problem, batch_case);
        if (completed.contains(key)) {
            std::cout << '[' << index + 1U << '/' << cases.size() << "] skip "
                      << batch_case.instance_path.filename().string() << '\n';
            continue;
        }
        std::cout << '[' << index + 1U << '/' << cases.size() << "] "
                  << batch_case.instance_path.filename().string() << std::endl;
        try {
            const Instance instance = read_instance(
                batch_case.instance_path, batch_case.graph_path,
                to_string(options.problem));
            const Solution solution = solve(instance, config);
            write_assignment(output_directory / reference, instance, solution);
            append_result_csv(result_path, key, batch_case.instance_path,
                              batch_case.graph_path, instance, solution,
                              reference);
            completed.insert(key);
            std::cout << "  status=" << to_string(solution.status)
                      << " LB=" << solution.lower_bound
                      << " UB=" << solution.upper_bound
                      << " gap=" << std::setprecision(8)
                      << solution.relative_gap
                      << " time=" << solution.stats.total_seconds << "s\n";
            if (!solution.assignment.complete()) {
                ++failures;
            }
#if PRECPACK_HAS_GUROBI
        } catch (const GRBException& error) {
            std::cerr << "  Gurobi error " << error.getErrorCode() << ": "
                      << error.getMessage() << '\n';
            if (strict_gurobi) {
                throw;
            }
            ++failures;
#endif
        } catch (const std::exception& error) {
            std::cerr << "  error: " << error.what() << '\n';
            ++failures;
        }
    }
    return failures == 0U ? 0 : 1;
}

}
