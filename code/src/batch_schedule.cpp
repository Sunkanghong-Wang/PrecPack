#include "batch_schedule.hpp"

#include <optional>
#include <stdexcept>
#include <string>

namespace precpack::internal {
namespace {

[[nodiscard]] std::optional<std::filesystem::path> relative_path_below(
    const std::filesystem::path& path,
    const std::filesystem::path& root) {
    const std::filesystem::path normalized_path =
        std::filesystem::weakly_canonical(path);
    const std::filesystem::path normalized_root =
        std::filesystem::weakly_canonical(root);
    const std::filesystem::path relative =
        normalized_path.lexically_relative(normalized_root);
    if (relative.empty() || relative == "." || relative.is_absolute()) {
        return std::nullopt;
    }
    const auto first = relative.begin();
    if (first == relative.end() || *first == "..") {
        return std::nullopt;
    }
    return relative;
}

[[nodiscard]] std::string first_directory(
    const std::filesystem::path& relative) {
    const auto first = relative.begin();
    return first == relative.end() ? std::string{} : first->string();
}

}

double effective_batch_time_limit(
    const CommandLineOptions& options,
    const BatchCase& batch_case,
    const std::filesystem::path& instance_root) {
    if (options.time_limit_was_set) {
        return options.time_limit_seconds;
    }

    const std::optional<std::filesystem::path> family_path =
        relative_path_below(batch_case.instance_path, instance_root);
    if (!family_path.has_value()) {
        return options.time_limit_seconds;
    }
    const std::string family = first_directory(*family_path);
    if (family == "scholl") {
        if (options.problem == ProblemKind::kSalbpI) {
            return 350.0;
        }
        if (options.problem == ProblemKind::kBppP) {
            return 1000.0;
        }
        throw std::invalid_argument(
            "the bundled Scholl set has no BPP-GP benchmark schedule; "
            "specify --time-limit explicitly");
    }
    if (family != "otto") {
        return options.time_limit_seconds;
    }

    const std::optional<std::filesystem::path> otto_path =
        relative_path_below(batch_case.instance_path,
                            instance_root / "otto");
    if (!otto_path.has_value()) {
        return options.time_limit_seconds;
    }
    const std::string benchmark_set = first_directory(*otto_path);
    if (options.problem == ProblemKind::kBppGp) {
        return 75.0;
    }
    if (options.problem == ProblemKind::kBppP) {
        if (benchmark_set == "n_0100") {
            return 1000.0;
        }
        if (benchmark_set == "n_0020" || benchmark_set == "n_0050" ||
            benchmark_set == "n_0250" || benchmark_set == "n_0500" ||
            benchmark_set == "n_0750" || benchmark_set == "n_1000") {
            return 75.0;
        }
    } else if (benchmark_set == "n_0020" ||
               benchmark_set == "n_0050" ||
               benchmark_set == "n_0050_permuted") {
        return 1000.0;
    } else if (benchmark_set == "n_0100" ||
               benchmark_set == "n_1000") {
        return 350.0;
    } else if (benchmark_set == "n_0250" ||
               benchmark_set == "n_0500" ||
               benchmark_set == "n_0750") {
        return 75.0;
    }

    throw std::invalid_argument(
        "no default benchmark time limit for " +
        std::string(to_string(options.problem)) + " set " + benchmark_set +
        "; specify --time-limit explicitly");
}

}
