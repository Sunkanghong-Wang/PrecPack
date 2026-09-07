#include "precpack/build_config.hpp"
#include "precpack/batch.hpp"
#include "precpack/cli.hpp"
#include "precpack/dff.hpp"
#include "precpack/exact_arithmetic.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/output_lock.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver_profile.hpp"

#include "batch_schedule.hpp"
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
    require(!options.time_limit_was_set,
            "the default time limit was marked as an explicit override");
    require(options.memory_limit_mb == 24ULL * 1024ULL,
            "public memory-limit default changed");
    require(options.output_directory ==
                std::filesystem::path("results/bpp-p"),
            "public output-directory default changed");
    require(!options.batch_mode,
            "single-instance mode unexpectedly enabled batch behavior");

    const precpack::CommandLineOptions explicit_limit = parse(
        {"precpack", "--problem", "bpp-p", "--instance", "case.txt",
         "--time-limit", "42"});
    require(explicit_limit.time_limit_was_set &&
                std::abs(explicit_limit.time_limit_seconds - 42.0) < 1e-12,
            "an explicit time-limit override was not recorded");
    require(precpack::make_command_line_solver_config(options, 42.0)
                .bbr_enable_root_strengthening,
            "the default BPP-P profile unexpectedly disabled root "
            "strengthening");
}

void test_cli_memory_units() {
    const precpack::CommandLineOptions options = parse(
        {"precpack", "--problem", "bpp-p", "--instance", "case.txt",
         "--memory-limit-mb", "512"});
    require(options.memory_limit_mb == 512U &&
                precpack::make_command_line_solver_config(options, 300.0)
                        .bbr_memory_limit_mb == 512U,
            "the existing memory-limit option changed its value or meaning");

    std::ostringstream output;
    precpack::print_help(output, "precpack");
    const std::string help = output.str();
    require(help.find("--memory-limit-mb MiB") != std::string::npos,
            "CLI help must label the memory-limit argument in MiB");
    require(help.find("24576 MiB = 24 GiB") != std::string::npos,
            "CLI help must state the default memory-limit conversion");
    require(help.find("1 MiB = 2^20 bytes; 1 GiB = 2^30 bytes") !=
                std::string::npos,
            "CLI help must define the binary memory units");
}

void test_bundled_benchmark_time_schedule() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path instance_root =
        temporary_directory.path() / "data" / "instances";

    const auto make_case = [&](const std::filesystem::path& relative) {
        const std::filesystem::path path =
            instance_root / relative / "case.txt";
        std::filesystem::create_directories(path.parent_path());
        write_text_file(path, "instance\n");
        return precpack::BatchCase{
            std::filesystem::weakly_canonical(path), std::nullopt};
    };
    const precpack::BatchCase scholl = make_case("scholl/Bowman");
    const precpack::BatchCase otto20 = make_case("otto/n_0020");
    const precpack::BatchCase otto50_permuted =
        make_case("otto/n_0050_permuted");
    const precpack::BatchCase otto100 = make_case("otto/n_0100");
    const precpack::BatchCase otto250 = make_case("otto/n_0250");
    const precpack::BatchCase otto1000 = make_case("otto/n_1000");
    const precpack::BatchCase external = [&] {
        const std::filesystem::path path =
            temporary_directory.path() / "external" / "case.txt";
        std::filesystem::create_directories(path.parent_path());
        write_text_file(path, "instance\n");
        return precpack::BatchCase{
            std::filesystem::weakly_canonical(path), std::nullopt};
    }();

    const auto limit = [&](precpack::ProblemKind problem,
                           const precpack::BatchCase& batch_case) {
        precpack::CommandLineOptions options;
        options.problem = problem;
        return precpack::internal::effective_batch_time_limit(
            options, batch_case, instance_root);
    };

    require(limit(precpack::ProblemKind::kSalbpI, scholl) == 350.0 &&
                limit(precpack::ProblemKind::kSalbpI, otto20) == 1000.0 &&
                limit(precpack::ProblemKind::kSalbpI,
                      otto50_permuted) == 1000.0 &&
                limit(precpack::ProblemKind::kSalbpI, otto100) == 350.0 &&
                limit(precpack::ProblemKind::kSalbpI, otto250) == 75.0 &&
                limit(precpack::ProblemKind::kSalbpI, otto1000) == 350.0,
            "SALBP-I bundled benchmark time schedule changed");
    require(limit(precpack::ProblemKind::kBppP, scholl) == 1000.0 &&
                limit(precpack::ProblemKind::kBppP, otto20) == 75.0 &&
                limit(precpack::ProblemKind::kBppP, otto100) == 1000.0 &&
                limit(precpack::ProblemKind::kBppP, otto1000) == 75.0,
            "BPP-P bundled benchmark time schedule changed");
    require(limit(precpack::ProblemKind::kBppGp, otto20) == 75.0 &&
                limit(precpack::ProblemKind::kBppGp, otto1000) == 75.0,
            "BPP-GP bundled benchmark time schedule changed");
    require(limit(precpack::ProblemKind::kBppP, external) == 300.0,
            "custom batch input did not retain the product default");

    precpack::CommandLineOptions override_options;
    override_options.problem = precpack::ProblemKind::kBppP;
    override_options.time_limit_seconds = 42.0;
    override_options.time_limit_was_set = true;
    require(precpack::internal::effective_batch_time_limit(
                override_options, otto100, instance_root) == 42.0,
            "explicit time limit did not override the benchmark schedule");
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

    bool rejected_removed_check_only = false;
    try {
        static_cast<void>(parse({"precpack", "--batch", "--problem",
                                 "bpp-p", "--check-only"}));
    } catch (const std::invalid_argument&) {
        rejected_removed_check_only = true;
    }
    require(rejected_removed_check_only,
            "the removed --check-only option remains publicly accepted");

    bool required_graph = false;
    try {
        static_cast<void>(parse({"precpack", "--problem", "bpp-gp",
                                 "--instance", "case.txt"}));
    } catch (const std::invalid_argument&) {
        required_graph = true;
    }
    require(required_graph, "BPP-GP did not require its GRAPH file");

    const precpack::CommandLineOptions batch = parse(
        {"precpack", "--batch", "--problem", "bpp-gp", "--input",
         "items", "--graph-dir", "graphs"});
    require(batch.batch_mode &&
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

void test_batch_resume_by_solution() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path instance_path =
        temporary_directory.path() / "items" / "case.txt";
    const std::filesystem::path output_directory =
        temporary_directory.path() / "results";
    std::filesystem::create_directories(instance_path.parent_path());
    write_text_file(instance_path, "instance");
    const std::filesystem::path canonical_instance =
        std::filesystem::weakly_canonical(instance_path);
    const std::filesystem::path solution_reference =
        std::filesystem::path("solutions") / "items" / "case.sol";
    std::filesystem::create_directories(
        (output_directory / solution_reference).parent_path());
    write_text_file(output_directory / solution_reference, "Bin 1: 1\n");
    const std::filesystem::path result_path =
        output_directory / "BPP-P_Results.csv";
    write_text_file(result_path, "old_schema\n");

    precpack::CommandLineOptions options;
    options.batch_mode = true;
    options.problem = precpack::ProblemKind::kBppP;
    options.input_path = canonical_instance;
    options.output_directory = output_directory;
    options.time_limit_seconds = 60.0;
    options.memory_limit_mb = 512U;
    require(precpack::run_batch(options) == 0,
            "an existing solution was not skipped");
    require(read_text_file(result_path) == "old_schema\n",
            "solution-file resume inspected or changed an unused CSV");

    options.time_limit_seconds = 1.0;
    options.memory_limit_mb = 64U;
    require(precpack::run_batch(options) == 0,
            "resource settings prevented solution-file resume");
}

void test_batch_caller_directory() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path caller =
        temporary_directory.path() / "caller directory";
    std::filesystem::create_directories(caller);
    write_text_file(caller / "local-instance.txt", "instance");
    const std::filesystem::path output_directory = caller / "relative output";
    const std::filesystem::path existing_solution =
        output_directory / "solutions" / caller.filename() /
        "local-instance.sol";
    std::filesystem::create_directories(existing_solution.parent_path());
    write_text_file(existing_solution, "Bin 1: 1\n");

    precpack::CommandLineOptions options;
    options.batch_mode = true;
    options.problem = precpack::ProblemKind::kBppP;
    options.input_path = "local-instance.txt";
    options.output_directory = "relative output";

    const std::string caller_string = caller.string();
    ScopedEnvironment caller_environment(
        "PRECPACK_CALLER_DIRECTORY", caller_string.c_str());
    require(precpack::run_batch(options) == 0,
            "batch paths were not resolved from the launcher caller");
    require(!std::filesystem::exists(output_directory / "BPP-P_Results.csv"),
            "a fully resumed batch unexpectedly created a result CSV");
    require(!std::filesystem::exists(output_directory / "errors.log"),
            "a batch without errors unexpectedly created errors.log");
}

void test_batch_failure_output() {
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
    require(precpack::run_batch(options) == 1,
            "an invalid batch instance did not report failure");

    const std::filesystem::path error_log = output_directory / "errors.log";
    require(std::filesystem::is_regular_file(error_log),
            "batch failure did not create errors.log");
    const std::string error_text = read_text_file(error_log);
    require(error_text.find("invalid-instance.txt") != std::string::npos &&
                error_text.find("error") != std::string::npos,
            "batch failure log omitted the instance or diagnostic");
    require(!std::filesystem::exists(output_directory / "logs"),
            "batch failure created a routine logs directory");
    require(!std::filesystem::exists(output_directory / ".precpack.lock"),
            "batch failure created a standalone lock file");
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
    require(config.time_limit_seconds == 17.0 &&
                config.bbr_memory_limit_mb == 512,
            "resource limits were not forwarded");
    require(config.bbr_enable_jackson &&
                config.bbr_enable_generalized_item_dominance,
            "BPP-GP item-dominance profile changed");
    require(config.bbr_enable_paper_queue_order &&
                config.bbr_enable_complete_dff,
            "BPP-GP 2016 queue/DFF profile changed");
    require(config.bbr_enable_root_strengthening,
            "BPP-GP position-free root profile changed");
    require(config.bbr_enable_binlb &&
                config.bbr_enable_conflict_binlb &&
                config.bbr_binlb_total_time_limit_seconds == 0.1 &&
                config.bbr_conflict_binlb_call_time_limit_seconds == 0.005 &&
                config.bbr_conflict_binlb_node_limit == 50'000U,
            "BPP-GP bounded conflict-aware BINLB profile changed");
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
                config.bbr_enable_conflict_binlb &&
                config.bbr_binlb_call_time_limit_seconds == 1.0 &&
                config.bbr_binlb_total_time_limit_seconds == 0.1 &&
                config.bbr_binlb_node_limit == 1'000'000U &&
                config.bbr_binlb_load_limit == 50U &&
                config.bbr_conflict_binlb_call_time_limit_seconds == 0.005 &&
                config.bbr_conflict_binlb_node_limit == 50'000U,
            "SALBP-I bounded BINLB profile changed");
    require(config.bbr_enable_complete_dff &&
                config.bbr_dff_transform_limit == 15,
            "SALBP-I complete DFF profile changed");
    require(config.bbr_enable_closure_bound,
            "SALBP-I residual closure bound was disabled");
    require(!config.bbr_enable_initial_bdp &&
                !config.bbr_enable_paper_queue_order,
            "SALBP-I enabled an excluded component");
    require(config.bbr_enable_root_strengthening,
            "SALBP-I universal position-free root bound was disabled");
}

void test_output_schema() {
    require(precpack::make_instance_set(
                "data/instances/otto/n_0020/case.txt", std::nullopt) ==
                "otto/n_0020",
            "Otto instance-set name changed");
    require(precpack::make_instance_set(
                "data/instances/otto/n_0020/case.txt",
                std::filesystem::path(
                    "data/bpp-gp-graphs/separation-03/n_0020/case.graph")) ==
                "separation-03/n_0020",
            "BPP-GP instance-set name changed");
    require(precpack::make_instance_set(
                "external-a/case.txt", std::nullopt) == "external-a",
            "external instance-set name changed");

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
    solution.bbr_stats.time_limit_seconds = 60.0;
    solution.bbr_stats.memory_limit_bytes = 512ULL * 1024ULL * 1024ULL;
    solution.bbr_stats.peak_memory_bytes = 123'456ULL;
    solution.bbr_stats.states_created = 987'654ULL;
    solution.assignment.bin_of_item = {0, 0, 1};
    solution.assignment.bin_count = 2;

    precpack::write_assignment(assignment_path, instance, solution);
    std::filesystem::path temporary_assignment_path = assignment_path;
    temporary_assignment_path += ".tmp";
    require(!std::filesystem::exists(temporary_assignment_path),
            "completed assignment write left a temporary file");
    precpack::append_result_csv(
        csv_path, "otto/n_0020", "case,1", instance, solution);

    const std::string expected_csv =
        "instance_set,instance,n,time_limit_seconds,memory_limit_mb,status,"
        "opt,lower_bound,upper_bound,time_seconds,"
        "bbr_peak_memory_bytes,bbr_states_created\n"
        "otto/n_0020,\"case,1\",3,60,512,OPTIMAL,1,2,2,1.25,"
        "123456,987654\n";
    require(read_text_file(csv_path) == expected_csv,
            "result CSV schema or serialization changed");
    const std::filesystem::path limited_csv_path =
        temporary_directory.path() / "limited.csv";
    precpack::Solution limited = solution;
    limited.status = precpack::SolveStatus::kTimeLimit;
    limited.optimal = false;
    limited.lower_bound = 1;
    precpack::append_result_csv(
        limited_csv_path, "otto/n_0020", "case-2", instance, limited);
    require(read_text_file(limited_csv_path).find(
                ",TIME_LIMIT,0,1,2,") != std::string::npos,
            "a limited result was not serialized with opt=0");
    require(read_text_file(assignment_path) ==
                "Bin 1: 1 2\nBin 2: 3\n",
            "assignment file contains redundant metadata or changed format");
    require(precpack::has_nonempty_solution(assignment_path),
            "nonempty solution was not detected");

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
            incompatible_path, "otto/n_0020", "case", instance, solution);
    } catch (const std::runtime_error&) {
        rejected_incompatible_header = true;
    }
    require(rejected_incompatible_header,
            "an incompatible result CSV schema was silently appended");
}

void test_output_lock() {
    TemporaryDirectory temporary_directory;
    const std::filesystem::path result_path =
        temporary_directory.path() / "results" / "BPP-P_Results.csv";
    {
        precpack::OutputLock first(result_path);
        std::ofstream output(result_path, std::ios::app | std::ios::binary);
        require(static_cast<bool>(output),
                "the result lock blocked its owning writer");
        output << "test\n";
        output.flush();
        require(static_cast<bool>(output),
                "the result lock blocked a result write");
        output.close();
        require(static_cast<bool>(output),
                "the locked result CSV could not be closed");
        bool rejected_second = false;
        try {
            precpack::OutputLock second(result_path);
        } catch (const std::runtime_error&) {
            rejected_second = true;
        }
        require(rejected_second,
                "two writers acquired the same result-CSV lock");
    }
    precpack::OutputLock reacquired(result_path);
    require(std::filesystem::is_regular_file(result_path) &&
                read_text_file(result_path) == "test\n" &&
                !std::filesystem::exists(
                    result_path.parent_path() / ".precpack.lock"),
            "result locking created a standalone lock file");
}

}

int main() {
    try {
        test_public_defaults();
        test_cli_memory_units();
        test_bundled_benchmark_time_schedule();
        test_exact_arithmetic();
        test_public_validation();
        test_batch_pairing();
        test_batch_resume_by_solution();
        test_batch_caller_directory();
        test_batch_failure_output();
        test_instance_file_validation();
        test_bpp_profile();
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
