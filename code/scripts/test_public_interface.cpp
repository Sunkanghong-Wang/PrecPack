#include "precpack/build_config.hpp"
#include "precpack/batch.hpp"
#include "precpack/cli.hpp"
#include "precpack/dff.hpp"
#include "precpack/exact_arithmetic.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/output_lock.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver_profile.hpp"

#include "environment.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

class TemporaryDirectory {
public:
    TemporaryDirectory() {
        const std::filesystem::path base =
            std::filesystem::temp_directory_path();
        const auto stamp = static_cast<std::uint64_t>(
            std::chrono::steady_clock::now().time_since_epoch().count());
        for (int attempt = 0; attempt < 100; ++attempt) {
            const std::filesystem::path candidate =
                base / ("precpack-interface-test-" + std::to_string(stamp) +
                        "-" + std::to_string(attempt));
            std::error_code error;
            if (std::filesystem::create_directory(candidate, error)) {
                path_ = candidate;
                return;
            }
            if (error) {
                throw std::filesystem::filesystem_error(
                    "cannot create temporary test directory", candidate, error);
            }
        }
        throw std::runtime_error("cannot allocate a unique test directory");
    }

    ~TemporaryDirectory() {
        std::error_code error;
        static_cast<void>(std::filesystem::remove_all(path_, error));
    }

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return path_;
    }

private:
    std::filesystem::path path_;
};

class ScopedEnvironment {
public:
    ScopedEnvironment(const char* name, const char* value) : name_(name) {
        const std::optional<std::string> previous =
            precpack::internal::environment_value(name);
        if (previous.has_value()) {
            previous_ = *previous;
        }
#ifdef _WIN32
        require(_putenv_s(name, value) == 0,
                "cannot set test environment variable");
#else
        require(setenv(name, value, 1) == 0,
                "cannot set test environment variable");
#endif
    }

    ~ScopedEnvironment() {
#ifdef _WIN32
        static_cast<void>(_putenv_s(name_.c_str(),
                                    previous_.has_value()
                                        ? previous_->c_str()
                                        : ""));
#else
        if (previous_.has_value()) {
            static_cast<void>(setenv(name_.c_str(), previous_->c_str(), 1));
        } else {
            static_cast<void>(unsetenv(name_.c_str()));
        }
#endif
    }

    ScopedEnvironment(const ScopedEnvironment&) = delete;
    ScopedEnvironment& operator=(const ScopedEnvironment&) = delete;

private:
    std::string name_;
    std::optional<std::string> previous_;
};

void write_text_file(const std::filesystem::path& path,
                     std::string_view contents) {
    std::ofstream output(path, std::ios::binary);
    require(static_cast<bool>(output), "cannot create test file: " + path.string());
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    require(static_cast<bool>(output), "cannot write test file: " + path.string());
}

[[nodiscard]] std::string read_text_file(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    require(static_cast<bool>(input), "cannot read test file: " + path.string());
    std::ostringstream contents;
    contents << input.rdbuf();
    require(static_cast<bool>(input) || input.eof(),
            "cannot read test file: " + path.string());
    return contents.str();
}

precpack::CommandLineOptions parse(std::vector<std::string> arguments) {
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) {
        argv.push_back(argument.data());
    }
    return precpack::parse_command_line(static_cast<int>(argv.size()),
                                        argv.data());
}

void test_public_defaults() {
    const precpack::CommandLineOptions options = parse(
        {"precpack", "--problem", "bpp-p", "--instance", "case.txt"});
    require(options.problem == precpack::ProblemKind::kBppP,
            "BPP-P problem parsing failed");
    require(options.instance_path == "case.txt",
            "instance path parsing failed");
    require(!options.graph_path.has_value(),
            "BPP-P unexpectedly accepted a graph by default");
    require(std::abs(options.time_limit_seconds - 300.0) < 1e-12,
            "public time-limit default changed");
    require(options.memory_limit_mb == 24ULL * 1024ULL,
            "public memory-limit default changed");
    require(options.threads == 1, "public thread default changed");
    require(options.output_directory == "results",
            "public output-directory default changed");
    require(!options.batch_mode && !options.check_only,
            "single-instance mode unexpectedly enabled batch behavior");
}

void test_exact_arithmetic() {
    using precpack::exact_arithmetic::compare_nonnegative_fractions;
    require(compare_nonnegative_fractions(1, 2, 2, 3) < 0,
            "exact fraction ordering failed");

    constexpr std::uint64_t maximum =
        std::numeric_limits<std::uint64_t>::max();
    require(compare_nonnegative_fractions(maximum - 1, maximum,
                                          maximum - 2, maximum - 1) > 0,
            "overflow-safe fraction ordering failed");
    require(compare_nonnegative_fractions(maximum - 1, maximum - 1,
                                          maximum, maximum) == 0,
            "overflow-safe fraction equality failed");

    precpack::exact_arithmetic::NonnegativeRatioSum sum(7);
    sum.add(6);
    sum.add(8);
    sum.add(1);
    require(sum.ceil_to_int() == 3, "exact ratio accumulation failed");

    require(precpack::exact_arithmetic::ceil_ratio_to_int(
                1, 2, "ceil test") == 1 &&
                precpack::exact_arithmetic::ceil_ratio_to_int(
                    -1, 2, "ceil test") == 0,
            "signed ceiling division failed");
    require(precpack::exact_arithmetic::ceil_nonnegative_product_ratio(
                11, 3, 5, "product ratio test") == 7,
            "exact product-ratio ceiling failed");

    bool detected_overflow = false;
    try {
        static_cast<void>(precpack::exact_arithmetic::checked_multiply(
            std::numeric_limits<std::int64_t>::max(), 2,
            "expected multiplication overflow"));
    } catch (const std::overflow_error&) {
        detected_overflow = true;
    }
    require(detected_overflow, "signed multiplication overflow was missed");

    const int maximum_int = std::numeric_limits<int>::max();
    const precpack::DffTransformSet transforms =
        precpack::build_complete_dff_transforms(
            {maximum_int, maximum_int - 1, maximum_int / 2},
            maximum_int, false);
    require(!transforms.capacities.empty(),
            "large-integer DFF construction failed");
}

void test_public_validation() {
    bool rejected_removed_switch = false;
    try {
        static_cast<void>(parse({"precpack", "--problem", "bpp-p",
                                 "--instance", "case.txt",
                                 "--bbr-paper-queue-order"}));
    } catch (const std::invalid_argument&) {
        rejected_removed_switch = true;
    }
    require(rejected_removed_switch,
            "a removed algorithm switch remains publicly accepted");

    bool required_graph = false;
    try {
        static_cast<void>(parse({"precpack", "--problem", "bpp-gp",
                                 "--instance", "case.txt"}));
    } catch (const std::invalid_argument&) {
        required_graph = true;
    }
    require(required_graph, "BPP-GP did not require its GRAPH file");

    const precpack::CommandLineOptions parallel = parse(
        {"precpack", "--problem", "bpp-p", "--instance", "case.txt",
         "--threads", "-1"});
    require(parallel.threads == -1, "automatic thread selection was not parsed");
    bool rejected_zero_threads = false;
    try {
        static_cast<void>(parse({"precpack", "--problem", "bpp-p",
                                 "--instance", "case.txt", "--threads", "0"}));
    } catch (const std::invalid_argument&) {
        rejected_zero_threads = true;
    }
    require(rejected_zero_threads, "zero worker threads were accepted");

    const precpack::CommandLineOptions batch = parse(
        {"precpack", "--batch", "--problem", "bpp-gp", "--input",
         "items", "--graph-dir", "graphs", "--check-only"});
    require(batch.batch_mode && batch.check_only &&
                batch.input_path == std::filesystem::path("items") &&
                batch.graph_directory == std::filesystem::path("graphs") &&
                batch.output_directory ==
                    std::filesystem::path("results/bpp-gp"),
            "batch command-line parsing failed");

    bool rejected_mixed_inputs = false;
    try {
        static_cast<void>(parse(
            {"precpack", "--batch", "--problem", "bpp-p", "--instance",
             "case.txt"}));
    } catch (const std::invalid_argument&) {
        rejected_mixed_inputs = true;
    }
    require(rejected_mixed_inputs,
            "batch mode accepted the single-instance interface");

    bool rejected_empty_batch_input = false;
    try {
        static_cast<void>(parse(
            {"precpack", "--batch", "--problem", "bpp-p", "--input", ""}));
    } catch (const std::invalid_argument&) {
        rejected_empty_batch_input = true;
    }
    require(rejected_empty_batch_input,
            "batch mode accepted an empty input path");
}

void test_batch_pairing() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path& root = temporary_directory.path();
    const std::filesystem::path instance_root = root / "instances";
    const std::filesystem::path graph_root = root / "graphs";
    const std::filesystem::path instance =
        instance_root / "otto" / "n_0020" / "case.txt";
    const std::filesystem::path graph_01 =
        graph_root / "separation-01" / "n_0020" / "case.graph";
    const std::filesystem::path graph_03 =
        graph_root / "separation-03" / "n_0020" / "case.graph";
    std::filesystem::create_directories(instance.parent_path());
    std::filesystem::create_directories(graph_01.parent_path());
    std::filesystem::create_directories(graph_03.parent_path());
    write_text_file(instance, "instance");
    write_text_file(graph_01, "graph");
    write_text_file(graph_03, "graph");

    const auto canonical = [](const std::filesystem::path& path) {
        return std::filesystem::weakly_canonical(path);
    };
    const std::vector<precpack::BatchCase> file_pair =
        precpack::collect_batch_cases(
            precpack::ProblemKind::kBppGp, instance, graph_01,
            instance_root, graph_root);
    require(file_pair.size() == 1U &&
                file_pair.front().instance_path == canonical(instance) &&
                file_pair.front().graph_path == canonical(graph_01),
            "single BPP-GP instance/graph pairing failed");

    const std::vector<precpack::BatchCase> file_with_graph_set =
        precpack::collect_batch_cases(
            precpack::ProblemKind::kBppGp, instance, graph_root,
            instance_root, graph_root);
    require(file_with_graph_set.size() == 2U,
            "single BPP-GP instance/graph-set pairing failed");

    const std::vector<precpack::BatchCase> size_pair =
        precpack::collect_batch_cases(
            precpack::ProblemKind::kBppGp,
            instance_root / "otto" / "n_0020",
            graph_root / "separation-01" / "n_0020", instance_root,
            graph_root);
    require(size_pair.size() == 1U &&
                size_pair.front().graph_path == canonical(graph_01),
            "same-size BPP-GP directory pairing failed");

    const std::vector<precpack::BatchCase> set_with_graph_file =
        precpack::collect_batch_cases(
            precpack::ProblemKind::kBppGp, instance_root / "otto", graph_01,
            instance_root, graph_root);
    require(set_with_graph_file.size() == 1U &&
                set_with_graph_file.front().instance_path ==
                    canonical(instance),
            "BPP-GP instance-set/single-graph pairing failed");

    const std::vector<precpack::BatchCase> default_pairs =
        precpack::collect_batch_cases(precpack::ProblemKind::kBppGp,
                                      std::nullopt, std::nullopt,
                                      instance_root, graph_root);
    require(default_pairs.size() == 2U &&
                default_pairs[0].graph_path == canonical(graph_01) &&
                default_pairs[1].graph_path == canonical(graph_03),
            "default BPP-GP benchmark-set pairing failed");

    bool rejected_mismatch = false;
    const std::filesystem::path other_graph =
        graph_root / "separation-01" / "n_0020" / "other.graph";
    write_text_file(other_graph, "graph");
    try {
        static_cast<void>(precpack::collect_batch_cases(
            precpack::ProblemKind::kBppGp, instance, other_graph,
            instance_root, graph_root));
    } catch (const std::invalid_argument&) {
        rejected_mismatch = true;
    }
    require(rejected_mismatch,
            "mismatched BPP-GP instance and graph were accepted");
}

void test_batch_resume_profile() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path instance_path =
        temporary_directory.path() / "items" / "case.txt";
    const std::filesystem::path output_directory =
        temporary_directory.path() / "results";
    std::filesystem::create_directories(instance_path.parent_path());
    write_text_file(instance_path, "instance");
    const std::filesystem::path canonical_instance =
        std::filesystem::weakly_canonical(instance_path);
    const std::string key =
        precpack::make_instance_key(canonical_instance, std::nullopt);
    const std::filesystem::path solution_reference =
        std::filesystem::path("solutions") / ("bpp-p__" + key + ".sol");
    std::filesystem::create_directories(
        (output_directory / solution_reference).parent_path());
    write_text_file(output_directory / solution_reference, "Bin 1: 1\n");

    precpack::Instance instance;
    instance.problem_type = "BPP-P";
    instance.capacity = 10;
    instance.items = {{0, 5}};
    precpack::Solution solution;
    solution.status = precpack::SolveStatus::kOptimal;
    solution.optimal = true;
    solution.lower_bound = 1;
    solution.upper_bound = 1;
    solution.threads = 1;
    solution.bbr_stats.time_limit_seconds = 60.0;
    solution.bbr_stats.memory_limit_bytes = 512ULL * 1024ULL * 1024ULL;
    precpack::append_result_csv(
        output_directory / "BPP-P_Results.csv", key, canonical_instance,
        std::nullopt, instance, solution, solution_reference);

    precpack::CommandLineOptions options;
    options.batch_mode = true;
    options.check_only = true;
    options.problem = precpack::ProblemKind::kBppP;
    options.input_path = canonical_instance;
    options.output_directory = output_directory;
    options.time_limit_seconds = 60.0;
    options.memory_limit_mb = 512U;
    options.threads = 1;
    require(precpack::run_batch(options) == 0,
            "a compatible completed batch result was not resumed");

    options.threads = 2;
    bool rejected_profile_mismatch = false;
    try {
        static_cast<void>(precpack::run_batch(options));
    } catch (const std::runtime_error&) {
        rejected_profile_mismatch = true;
    }
    require(rejected_profile_mismatch,
            "batch resume accepted a mismatched thread profile");

    options.threads = 1;
    bool rejected_gurobi_profile = false;
    {
        ScopedEnvironment strict_gurobi(
            "PRECPACK_REQUIRE_GUROBI_RUNTIME", "1");
        try {
            static_cast<void>(precpack::run_batch(options));
        } catch (const std::runtime_error&) {
            rejected_gurobi_profile = true;
        }
    }
    require(rejected_gurobi_profile,
            "batch resume mixed strict and optional Gurobi profiles");
}

void test_batch_caller_directory() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path caller =
        temporary_directory.path() / "caller directory";
    std::filesystem::create_directories(caller);
    write_text_file(caller / "local-instance.txt", "instance");

    precpack::CommandLineOptions options;
    options.batch_mode = true;
    options.check_only = true;
    options.problem = precpack::ProblemKind::kBppP;
    options.input_path = "local-instance.txt";
    options.output_directory = "relative output";

    const std::string caller_string = caller.string();
    ScopedEnvironment caller_environment(
        "PRECPACK_CALLER_DIRECTORY", caller_string.c_str());
    require(precpack::run_batch(options) == 0,
            "batch paths were not resolved from the launcher caller");
    require(!std::filesystem::exists(caller / "relative output"),
            "check-only batch unexpectedly created its output directory");
}

void test_batch_failure_logging() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path instance_path =
        temporary_directory.path() / "invalid-instance.txt";
    const std::filesystem::path output_directory =
        temporary_directory.path() / "results";
    write_text_file(instance_path, "this is not a PrecPack instance\n");

    precpack::CommandLineOptions options;
    options.batch_mode = true;
    options.problem = precpack::ProblemKind::kBppP;
    options.input_path = instance_path;
    options.output_directory = output_directory;
    options.time_limit_seconds = 1.0;
    options.memory_limit_mb = 64U;
    options.threads = 1;
    require(precpack::run_batch(options) == 1,
            "an invalid batch instance did not report failure");

    const std::filesystem::path canonical_instance =
        std::filesystem::weakly_canonical(instance_path);
    const std::string key =
        precpack::make_instance_key(canonical_instance, std::nullopt);
    const std::filesystem::path event_path =
        output_directory / "logs" / "batch-events.log";
    const std::filesystem::path failure_path =
        output_directory / "logs" /
        ("bpp-p__" + key + ".failure.log");
    const std::string events = read_text_file(event_path);
    require(events.find("event=BATCH_START") != std::string::npos &&
                events.find("event=START instance_key=" + key) !=
                    std::string::npos &&
                events.find("event=ERROR instance_key=" + key) !=
                    std::string::npos &&
                events.find("event=BATCH_END failures=1") !=
                    std::string::npos,
            "batch event log does not preserve the failed attempt");
    const std::string failure = read_text_file(failure_path);
    require(failure.find("event=ERROR") != std::string::npos &&
                failure.find("instance_key=" + key) != std::string::npos &&
                failure.find("threads=1") != std::string::npos &&
                failure.find("exception=std::exception") !=
                    std::string::npos &&
                failure.find("message=") != std::string::npos,
            "per-instance failure log is incomplete");
}

void test_instance_file_validation() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path& data = temporary_directory.path();
    write_text_file(data / "simple.txt", R"(<number of tasks>
4

<cycle time>
10

<order strength>
0.6

<task times>
1 6
2 4
3 6
4 4

<precedence relations>

<end>
)");
    write_text_file(data / "generalized.graph", "1\n1 4 2\n");
    write_text_file(data / "invalid_graph_missing_count.graph", "1 4 2\n");
    write_text_file(data / "invalid_graph_count_mismatch.graph", "2\n1 4 2\n");
    write_text_file(data / "invalid_graph_extra_column.graph", "1\n1 4 2 9\n");
    write_text_file(data / "invalid_task_count.txt", R"(<number of tasks>
4 tasks
<cycle time>
10
<task times>
1 6
2 4
3 6
4 4
<end>
)");
    write_text_file(data / "invalid_cycle_time.txt", R"(<number of tasks>
4
<cycle time>
10 seconds
<task times>
1 6
2 4
3 6
4 4
<end>
)");
    write_text_file(data / "invalid_order_strength.txt", R"(<number of tasks>
4
<cycle time>
10
<order strength>
0.6 extra
<task times>
1 6
2 4
3 6
4 4
<end>
)");
    write_text_file(data / "duplicate_task_id.txt", R"(<number of tasks>
4
<cycle time>
10
<task times>
1 6
1 4
3 6
4 4
<end>
)");

    const precpack::Instance instance = precpack::read_instance(
        data / "simple.txt", data / "generalized.graph", "BPP-GP-01", 7);
    require(instance.size() == 4, "parser: wrong item count");
    require(instance.capacity == 10, "parser: wrong capacity");
    require(instance.order_strength == "0.6",
            "parser: canonical order-strength metadata mismatch");
    require(instance.arcs.size() == 1U,
            "graph override did not replace ALB arcs");
    require(instance.arcs.front().from == 0 && instance.arcs.front().to == 3 &&
                instance.arcs.front().separation == 2,
            "parser: generalized arc mismatch");
    require(instance.front[3] == 2 && instance.back[0] == 2,
            "graph preprocessing: position bounds mismatch");

    for (const char* filename : {
             "invalid_graph_missing_count.graph",
             "invalid_graph_count_mismatch.graph",
             "invalid_graph_extra_column.graph",
         }) {
        bool rejected = false;
        try {
            static_cast<void>(precpack::read_instance(
                data / "simple.txt", data / filename, "BPP-GP"));
        } catch (const std::exception&) {
            rejected = true;
        }
        require(rejected, std::string("malformed graph was accepted: ") +
                              filename);
    }

    for (const char* filename : {
             "invalid_task_count.txt",
             "invalid_cycle_time.txt",
             "invalid_order_strength.txt",
             "duplicate_task_id.txt",
         }) {
        bool rejected = false;
        try {
            static_cast<void>(precpack::read_instance(
                data / filename, std::nullopt, "BPP-P"));
        } catch (const std::exception&) {
            rejected = true;
        }
        require(rejected, std::string("malformed ALB file was accepted: ") +
                              filename);
    }
}

void test_bpp_profile() {
    const precpack::Config config = precpack::make_solver_config(
        precpack::ProblemKind::kBppGp, 17.0, 512);
    require(config.seed == 1,
            "fixed reproducibility controls changed");
    require(config.threads == 1, "serial profile default changed");
    require(config.time_limit_seconds == 17.0 &&
                config.bbr_memory_limit_mb == 512,
            "resource limits were not forwarded");
    require(config.bbr_enable_jackson &&
                config.bbr_enable_generalized_item_dominance,
            "BPP-GP item-dominance profile changed");
    require(config.bbr_enable_paper_queue_order &&
                config.bbr_enable_complete_dff,
            "BPP-GP 2016 queue/DFF profile changed");
    require(config.bbr_enable_root_strengthening &&
                config.bbr_root_cg_time_limit_seconds == 5.0,
            "BPP-GP adaptive price-and-switch profile changed");
    require(!config.bbr_enable_binlb,
            "BPP-GP enabled an excluded component");
}

void test_parallel_profile() {
    const precpack::Config config = precpack::make_solver_config(
        precpack::ProblemKind::kBppGp, 17.0, 512, 8);
    require(config.threads == 8,
            "explicit shared-memory worker count was not forwarded");
    require(precpack::resolve_thread_count(8) == 8 &&
                precpack::resolve_thread_count(-1) >= 1,
            "thread-count resolution is invalid");
}

void test_salbp_profile() {
    const precpack::Config config = precpack::make_solver_config(
        precpack::ProblemKind::kSalbpI, 23.0, 1024);
    require(config.bbr_enable_jackson &&
                !config.bbr_enable_generalized_item_dominance,
            "SALBP-I item-dominance profile changed");
    require(config.bbr_enable_bbr12_mhh &&
                config.bbr_enable_bbr12_mhh_portfolio &&
                config.bbr12_mhh_portfolio_max_items == 200 &&
                config.bbr12_mhh_full_load_limit == 1000,
            "SALBP-I bounded MHH profile changed");
    require(config.bbr_enable_binlb &&
                config.bbr_binlb_call_time_limit_seconds == 1.0 &&
                config.bbr_binlb_total_time_limit_seconds == 0.1 &&
                config.bbr_binlb_node_limit == 1'000'000U &&
                config.bbr_binlb_load_limit == 50U,
            "SALBP-I bounded BINLB profile changed");
    require(config.bbr_enable_complete_dff &&
                config.bbr_dff_transform_limit == 15,
            "SALBP-I complete DFF profile changed");
    require(!config.bbr_enable_initial_bdp &&
                !config.bbr_enable_paper_queue_order &&
                !config.bbr_enable_closure_bound &&
                !config.bbr_enable_root_strengthening,
            "SALBP-I enabled an excluded component");
}

void test_output_schema() {
    require(precpack::make_instance_key("external-a/case.txt", std::nullopt) ==
                "case__3a95170ceb4e4d41",
            "instance-key path fingerprint changed");
    require(precpack::make_instance_key("external-a/case.txt", std::nullopt) !=
                precpack::make_instance_key("external-b/case.txt", std::nullopt),
            "same-stem instances from different directories share a key");
    require(precpack::make_instance_key(
                "data/items/case.txt",
                std::filesystem::path("data/graphs/separation-01/case.graph")) !=
                precpack::make_instance_key(
                    "data/items/case.txt",
                    std::filesystem::path("data/graphs/separation-03/case.graph")),
            "different BPP-GP graphs share an instance key");

    TemporaryDirectory temporary_directory;
    const std::filesystem::path csv_path =
        temporary_directory.path() / "BPP-P_Results.csv";
    const std::filesystem::path assignment_path =
        temporary_directory.path() / "solutions" / "case.sol";

    precpack::Instance instance;
    instance.problem_type = "BPP-P";
    instance.capacity = 10;
    instance.items = {{0, 6}, {1, 4}, {2, 7}};

    precpack::Solution solution;
    solution.status = precpack::SolveStatus::kOptimal;
    solution.optimal = true;
    solution.gurobi_runtime_required = true;
    solution.lower_bound = 2;
    solution.upper_bound = 2;
    solution.relative_gap = 0.0;
    solution.stats.total_seconds = 1.25;
    solution.threads = 4;
    solution.bbr_stats.time_limit_seconds = 60.0;
    solution.bbr_stats.memory_limit_bytes = 512ULL * 1024ULL * 1024ULL;
    solution.bbr_stats.peak_memory_bytes = 123'456ULL;
    solution.assignment.bin_of_item = {0, 0, 1};
    solution.assignment.bin_count = 2;

    precpack::write_assignment(assignment_path, instance, solution);
    std::filesystem::path temporary_assignment_path = assignment_path;
    temporary_assignment_path += ".tmp";
    require(!std::filesystem::exists(temporary_assignment_path),
            "completed assignment write left a temporary file");
    precpack::append_result_csv(
        csv_path, "case,1", "data/case,1.txt", std::nullopt, instance,
        solution, "solutions/case.sol");

    const std::string expected_csv =
        "instance_key,problem,instance_file,graph_file,n,capacity,status,"
        "lower_bound,upper_bound,gap,time_seconds,time_limit_seconds,threads,"
        "memory_limit_mb,bbr_peak_memory_bytes,gurobi_enabled,"
        "gurobi_required,solution_file\n"
        "\"case,1\",BPP-P,\"data/case,1.txt\",,3,10,OPTIMAL,2,2,0,1.25,"
        "60,4,512,123456," +
        std::string(precpack::kHasGurobiSupport ? "1" : "0") +
        ",1,solutions/case.sol\n";
    require(read_text_file(csv_path) == expected_csv,
            "result CSV schema or serialization changed");
    const std::vector<precpack::ResultReference> references =
        precpack::read_result_references(csv_path);
    require(references.size() == 1U &&
                references.front().instance_key == "case,1" &&
                references.front().problem == "BPP-P" &&
                std::abs(references.front().time_limit_seconds - 60.0) <
                    1e-12 &&
                references.front().threads == 4 &&
                references.front().memory_limit_mb == 512U &&
                references.front().gurobi_enabled ==
                    precpack::kHasGurobiSupport &&
                references.front().gurobi_required &&
                references.front().solution_file == "solutions/case.sol",
            "result CSV references were not parsed correctly");

    const std::filesystem::path malformed_result_path =
        temporary_directory.path() / "malformed-result.csv";
    std::string malformed_result = expected_csv;
    const std::size_t time_field = malformed_result.find(",60,4,512,");
    require(time_field != std::string::npos,
            "test result row no longer contains the expected profile fields");
    malformed_result.replace(time_field, 4U, ",,");
    write_text_file(malformed_result_path, malformed_result);
    bool rejected_malformed_result = false;
    try {
        static_cast<void>(
            precpack::read_result_references(malformed_result_path));
    } catch (const std::runtime_error&) {
        rejected_malformed_result = true;
    }
    require(rejected_malformed_result,
            "a result row with an empty time limit was accepted");
    require(read_text_file(assignment_path) ==
                "Bin 1: 1 2\nBin 2: 3\n",
            "assignment file contains redundant metadata or changed format");

    bool rejected_duplicate_key = false;
    try {
        precpack::require_unused_instance_key(csv_path, "case,1");
    } catch (const std::runtime_error&) {
        rejected_duplicate_key = true;
    }
    require(rejected_duplicate_key,
            "an existing single-instance result key was accepted");
    precpack::require_unused_instance_key(csv_path, "unused-case");

    precpack::Solution replacement = solution;
    replacement.assignment.bin_of_item = {0, 1, 2};
    replacement.assignment.bin_count = 3;
    precpack::write_assignment(assignment_path, instance, replacement);
    require(read_text_file(assignment_path) ==
                "Bin 1: 1\nBin 2: 2\nBin 3: 3\n" &&
                !std::filesystem::exists(temporary_assignment_path),
            "atomic assignment replacement failed");

    const std::filesystem::path incompatible_path =
        temporary_directory.path() / "incompatible.csv";
    write_text_file(incompatible_path, "different_header\n");
    bool rejected_incompatible_header = false;
    try {
        precpack::append_result_csv(
            incompatible_path, "case", "data/case.txt", std::nullopt,
            instance, solution, "solutions/case.sol");
    } catch (const std::runtime_error&) {
        rejected_incompatible_header = true;
    }
    require(rejected_incompatible_header,
            "an incompatible result CSV schema was silently appended");
}

void test_output_lock() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path output_directory =
        temporary_directory.path() / "results";
    {
        precpack::OutputLock first(output_directory);
        bool rejected_second = false;
        try {
            precpack::OutputLock second(output_directory);
        } catch (const std::runtime_error&) {
            rejected_second = true;
        }
        require(rejected_second,
                "two writers acquired the same output-directory lock");
    }
    precpack::OutputLock reacquired(output_directory);
    require(std::filesystem::is_regular_file(
                output_directory / ".precpack.lock"),
            "output lock file was not retained for race-free reuse");
}

}

int main() {
    try {
        test_public_defaults();
        test_exact_arithmetic();
        test_public_validation();
        test_batch_pairing();
        test_batch_resume_profile();
        test_batch_caller_directory();
        test_batch_failure_logging();
        test_instance_file_validation();
        test_bpp_profile();
        test_parallel_profile();
        test_salbp_profile();
        test_output_lock();
        test_output_schema();
        std::cout << "Public interface tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Public interface test failed: " << error.what() << '\n';
        return 1;
    }
}
