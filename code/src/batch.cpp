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
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
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

[[nodiscard]] std::string utc_timestamp() {
    const std::time_t now = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &now);
#else
    gmtime_r(&now, &utc);
#endif
    std::ostringstream output;
    output << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

[[nodiscard]] std::string single_line(std::string_view text) {
    std::string result;
    result.reserve(text.size());
    for (char character : text) {
        switch (character) {
            case '\r':
                result += "\\r";
                break;
            case '\n':
                result += "\\n";
                break;
            case '\t':
                result += "\\t";
                break;
            default:
                result.push_back(character);
                break;
        }
    }
    return result;
}

class BatchLog {
public:
    BatchLog(const std::filesystem::path& output_directory,
             ProblemKind problem,
             const Config& config,
             bool strict_gurobi)
        : log_directory_(output_directory / "logs"),
          problem_(problem),
          threads_(resolve_thread_count(config.threads)),
          time_limit_seconds_(config.time_limit_seconds),
          state_limit_(config.bbr_state_limit),
          memory_limit_mb_(config.bbr_memory_limit_mb),
          strict_gurobi_(strict_gurobi) {
        std::filesystem::create_directories(log_directory_);
        event_path_ = log_directory_ / "batch-events.log";
        events_.open(event_path_, std::ios::binary | std::ios::app);
        if (!events_) {
            throw std::runtime_error("cannot open batch event log: " +
                                     event_path_.string());
        }
    }

    BatchLog(const BatchLog&) = delete;
    BatchLog& operator=(const BatchLog&) = delete;

    void batch_start(std::size_t selected, std::size_t completed) {
        std::ostringstream details;
        details << "problem=" << to_string(problem_)
                << " selected=" << selected
                << " completed=" << completed
                << " pending=" << selected - completed
                << " threads=" << threads_
                << " time_limit_seconds=" << std::setprecision(12)
                << time_limit_seconds_
                << " state_limit=" << state_limit_
                << " memory_limit_mb=" << memory_limit_mb_
                << " gurobi_required=" << (strict_gurobi_ ? 1 : 0);
        event("BATCH_START", {}, details.str());
    }

    void attempt_start(std::size_t index,
                       std::size_t selected,
                       std::string_view key,
                       const BatchCase& batch_case) {
        std::ostringstream details;
        details << "position=" << index << '/' << selected
                << " instance_file="
                << std::quoted(batch_case.instance_path.generic_string())
                << " graph_file="
                << std::quoted(batch_case.graph_path.has_value()
                                   ? batch_case.graph_path->generic_string()
                                   : std::string());
        event("START", key, details.str());
    }

    void success(std::string_view key,
                 const BatchCase& batch_case,
                 const Solution& solution) {
        std::ostringstream details;
        details << "status=" << to_string(solution.status)
                << " lower_bound=" << solution.lower_bound
                << " upper_bound=" << solution.upper_bound
                << " gap=" << std::setprecision(12)
                << solution.relative_gap
                << " time_seconds=" << solution.stats.total_seconds;
        event("SUCCESS", key, details.str());
        append_recovery_if_needed(key, batch_case, solution);
    }

    [[nodiscard]] std::filesystem::path failure(
        std::string_view key,
        const BatchCase& batch_case,
        std::string_view exception_kind,
        std::string_view message) {
        const std::string normalized_message = single_line(message);
        std::ostringstream details;
        details << "exception=" << std::quoted(std::string(exception_kind))
                << " message=" << std::quoted(normalized_message);
        event("ERROR", key, details.str());

        const std::filesystem::path path = failure_path(key);
        std::ofstream output(path, std::ios::binary | std::ios::app);
        if (!output) {
            mark_unhealthy("cannot open per-instance failure log: " +
                           path.string());
            return path;
        }
        output << "---\n"
               << "timestamp=" << utc_timestamp() << '\n'
               << "event=ERROR\n"
               << "instance_key=" << key << '\n'
               << "problem=" << to_string(problem_) << '\n'
               << "instance_file="
               << batch_case.instance_path.generic_string() << '\n'
               << "graph_file="
               << (batch_case.graph_path.has_value()
                       ? batch_case.graph_path->generic_string()
                       : std::string())
               << '\n'
               << "threads=" << threads_ << '\n'
               << "time_limit_seconds=" << std::setprecision(12)
               << time_limit_seconds_ << '\n'
               << "state_limit=" << state_limit_ << '\n'
               << "memory_limit_mb=" << memory_limit_mb_ << '\n'
               << "gurobi_required=" << (strict_gurobi_ ? 1 : 0) << '\n'
               << "exception=" << exception_kind << '\n'
               << "message=" << normalized_message << '\n';
        output.flush();
        if (!output) {
            mark_unhealthy("cannot write per-instance failure log: " +
                           path.string());
        }
        return path;
    }

    void fatal(std::string_view exception_kind,
               std::string_view message) noexcept {
        try {
            const std::string normalized_message = single_line(message);
            std::ostringstream details;
            details << "exception="
                    << std::quoted(std::string(exception_kind))
                    << " message=" << std::quoted(normalized_message);
            event("FATAL", {}, details.str());
        } catch (...) {
            mark_unhealthy("cannot record fatal batch exception");
        }
    }

    void batch_end(std::size_t failures) {
        std::ostringstream details;
        details << "failures=" << failures
                << " logging_healthy=" << (healthy_ ? 1 : 0);
        event("BATCH_END", {}, details.str());
    }

    [[nodiscard]] bool healthy() const noexcept {
        return healthy_;
    }

    [[nodiscard]] const std::filesystem::path& event_path() const noexcept {
        return event_path_;
    }

private:
    [[nodiscard]] std::filesystem::path failure_path(
        std::string_view key) const {
        return log_directory_ /
               (std::string(to_slug(problem_)) + "__" + std::string(key) +
                ".failure.log");
    }

    void event(std::string_view name,
               std::string_view key,
               std::string_view details) {
        events_ << utc_timestamp() << " event=" << name;
        if (!key.empty()) {
            events_ << " instance_key=" << key;
        }
        if (!details.empty()) {
            events_ << ' ' << details;
        }
        events_ << '\n';
        events_.flush();
        if (!events_) {
            mark_unhealthy("cannot write batch event log: " +
                           event_path_.string());
        }
    }

    void append_recovery_if_needed(std::string_view key,
                                   const BatchCase& batch_case,
                                   const Solution& solution) {
        const std::filesystem::path path = failure_path(key);
        std::error_code error;
        if (!std::filesystem::is_regular_file(path, error)) {
            return;
        }
        std::ofstream output(path, std::ios::binary | std::ios::app);
        if (!output) {
            mark_unhealthy("cannot append recovery to failure log: " +
                           path.string());
            return;
        }
        output << "---\n"
               << "timestamp=" << utc_timestamp() << '\n'
               << "event=RECOVERED\n"
               << "instance_key=" << key << '\n'
               << "instance_file="
               << batch_case.instance_path.generic_string() << '\n'
               << "status=" << to_string(solution.status) << '\n'
               << "lower_bound=" << solution.lower_bound << '\n'
               << "upper_bound=" << solution.upper_bound << '\n'
               << "gap=" << std::setprecision(12)
               << solution.relative_gap << '\n'
               << "time_seconds=" << solution.stats.total_seconds << '\n';
        output.flush();
        if (!output) {
            mark_unhealthy("cannot append recovery to failure log: " +
                           path.string());
        }
    }

    void mark_unhealthy(const std::string& message) noexcept {
        healthy_ = false;
        if (!reported_logging_error_) {
            std::cerr << "  logging error: " << message << '\n';
            reported_logging_error_ = true;
        }
    }

    std::filesystem::path log_directory_;
    std::filesystem::path event_path_;
    ProblemKind problem_;
    int threads_ = 1;
    double time_limit_seconds_ = 0.0;
    std::uint64_t state_limit_ = 0U;
    std::uint64_t memory_limit_mb_ = 0U;
    bool strict_gurobi_ = false;
    bool healthy_ = true;
    bool reported_logging_error_ = false;
    std::ofstream events_;
};

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
    std::unique_ptr<BatchLog> batch_log;
    if (!options.check_only) {
        batch_log = std::make_unique<BatchLog>(
            output_directory, options.problem, config, strict_gurobi);
    }

    try {
        if (strict_gurobi && !options.check_only) {
            verify_gurobi_runtime();
        }
        std::vector<BatchCase> cases = collect_batch_cases(
            options.problem, input_path, graph_directory,
            root / "data/instances", root / "data/bpp-gp-graphs");
        if (cases.empty()) {
            throw std::invalid_argument("no instances selected");
        }
        const std::vector<ResultReference> references =
            read_result_references(result_path);
        std::map<std::string, ResultReference> recorded_results;
        for (const ResultReference& reference : references) {
            if (!recorded_results.emplace(reference.instance_key, reference)
                     .second) {
                throw std::runtime_error(
                    "duplicate instance_key in result CSV: " +
                    reference.instance_key);
            }
        }

        std::set<std::string> output_names;
        std::set<std::string> completed;
        std::size_t complete_count = 0U;
        for (const BatchCase& batch_case : cases) {
            const std::string key = make_instance_key(
                batch_case.instance_path, batch_case.graph_path);
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
        batch_log->batch_start(cases.size(), complete_count);
        std::cout << "  event log : " << batch_log->event_path().string()
                  << '\n';

        std::size_t failures = 0U;
        for (std::size_t index = 0U; index < cases.size(); ++index) {
            const BatchCase& batch_case = cases[index];
            const std::string key = make_instance_key(
                batch_case.instance_path, batch_case.graph_path);
            const std::filesystem::path reference =
                solution_reference(options.problem, batch_case);
            if (completed.contains(key)) {
                std::cout << '[' << index + 1U << '/' << cases.size()
                          << "] skip "
                          << batch_case.instance_path.filename().string()
                          << '\n';
                continue;
            }
            std::cout << '[' << index + 1U << '/' << cases.size() << "] "
                      << batch_case.instance_path.filename().string()
                      << std::endl;
            batch_log->attempt_start(index + 1U, cases.size(), key,
                                     batch_case);
            try {
                const Instance instance = read_instance(
                    batch_case.instance_path, batch_case.graph_path,
                    to_string(options.problem));
                const Solution solution = solve(instance, config);
                write_assignment(output_directory / reference, instance,
                                 solution);
                append_result_csv(result_path, key, batch_case.instance_path,
                                  batch_case.graph_path, instance, solution,
                                  reference);
                completed.insert(key);
                std::cout << "  status=" << to_string(solution.status)
                          << " LB=" << solution.lower_bound
                          << " UB=" << solution.upper_bound
                          << " gap=" << std::setprecision(8)
                          << solution.relative_gap
                          << " time=" << solution.stats.total_seconds
                          << "s\n";
                if (solution.assignment.complete()) {
                    batch_log->success(key, batch_case, solution);
                } else {
                    const std::filesystem::path failure_path =
                        batch_log->failure(
                            key, batch_case, "IncompleteAssignment",
                            "solver returned no complete validated assignment");
                    std::cerr << "  failure log: " << failure_path.string()
                              << '\n';
                    ++failures;
                }
#if PRECPACK_HAS_GUROBI
            } catch (const GRBException& error) {
                const std::string message =
                    "code=" + std::to_string(error.getErrorCode()) +
                    " message=" + error.getMessage();
                const std::filesystem::path failure_path =
                    batch_log->failure(key, batch_case, "GRBException",
                                       message);
                std::cerr << "  Gurobi error " << error.getErrorCode() << ": "
                          << error.getMessage() << '\n'
                          << "  failure log: " << failure_path.string()
                          << '\n';
                if (strict_gurobi) {
                    throw;
                }
                ++failures;
#endif
            } catch (const std::exception& error) {
                const std::filesystem::path failure_path =
                    batch_log->failure(key, batch_case, "std::exception",
                                       error.what());
                std::cerr << "  error: " << error.what() << '\n'
                          << "  failure log: " << failure_path.string()
                          << '\n';
                ++failures;
            } catch (...) {
                const std::filesystem::path failure_path =
                    batch_log->failure(key, batch_case, "unknown",
                                       "non-standard exception");
                std::cerr << "  error: non-standard exception\n"
                          << "  failure log: " << failure_path.string()
                          << '\n';
                ++failures;
            }
        }
        batch_log->batch_end(failures);
        return failures == 0U && batch_log->healthy() ? 0 : 1;
#if PRECPACK_HAS_GUROBI
    } catch (const GRBException& error) {
        if (batch_log) {
            const std::string message =
                "code=" + std::to_string(error.getErrorCode()) +
                " message=" + error.getMessage();
            batch_log->fatal("GRBException", message);
        }
        throw;
#endif
    } catch (const std::exception& error) {
        if (batch_log) {
            batch_log->fatal("std::exception", error.what());
        }
        throw;
    } catch (...) {
        if (batch_log) {
            batch_log->fatal("unknown", "non-standard exception");
        }
        throw;
    }
}

}
