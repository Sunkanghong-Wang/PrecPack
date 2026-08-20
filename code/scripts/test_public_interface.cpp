#include "precpack/build_config.hpp"
#include "precpack/cli.hpp"
#include "precpack/instance_io.hpp"
#include "precpack/result_io.hpp"
#include "precpack/solver_profile.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
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
}

void test_public_validation() {
    bool rejected_legacy_switch = false;
    try {
        static_cast<void>(parse({"precpack", "--problem", "bpp-p",
                                 "--instance", "case.txt",
                                 "--bbr-paper-queue-order"}));
    } catch (const std::invalid_argument&) {
        rejected_legacy_switch = true;
    }
    require(rejected_legacy_switch,
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
}

void test_bpp_profile() {
    const precpack::Config config = precpack::make_solver_config(
        precpack::ProblemKind::kBppGp, 17.0, 512);
    require(config.exact_method == precpack::ExactMethod::kBbr,
            "BPP-GP is not using unified BBR");
    require(config.seed == 1 && config.bbr_state_limit == 60'000'000ULL,
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
    require(config.bbr_root_cg_mode ==
                precpack::BbrRootCgMode::kPriceAndSwitch &&
                config.bbr_root_cg_time_limit_seconds == 5.0,
            "BPP-GP adaptive price-and-switch profile changed");
    require(!config.enable_initial_alns &&
                !config.bbr_enable_initial_alns &&
                !config.bbr_enable_binlb &&
                !config.bbr_enable_conflict_binlb,
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
    require(config.exact_method == precpack::ExactMethod::kBbr &&
                config.bbr_heuristic_load_limit == 0,
            "SALBP-I is not using exact BBR");
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
                config.bbr_root_cg_mode == precpack::BbrRootCgMode::kNone,
            "SALBP-I enabled an excluded component");
}

void test_output_schema() {
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
    solution.lower_bound = 2;
    solution.upper_bound = 2;
    solution.relative_gap = 0.0;
    solution.stats.total_seconds = 1.25;
    solution.threads = 4;
    solution.bbr_stats.time_limit_seconds = 60.0;
    solution.bbr_stats.configured_state_limit = 60'000'000ULL;
    solution.bbr_stats.memory_limit_bytes = 512ULL * 1024ULL * 1024ULL;
    solution.bbr_stats.peak_memory_bytes = 123'456ULL;
    solution.assignment.bin_of_item = {0, 0, 1};
    solution.assignment.bin_count = 2;

    precpack::write_assignment(assignment_path, instance, solution);
    precpack::append_result_csv(
        csv_path, "case,1", "data/case,1.txt", std::nullopt, instance,
        solution, "solutions/case.sol");

    const std::string expected_csv =
        "instance_key,problem,instance_file,graph_file,n,capacity,status,"
        "lower_bound,upper_bound,gap,time_seconds,time_limit_seconds,threads,"
        "state_limit,memory_limit_mb,bbr_peak_memory_bytes,gurobi_enabled,"
        "solution_file\n"
        "\"case,1\",BPP-P,\"data/case,1.txt\",,3,10,OPTIMAL,2,2,0,1.25,"
        "60,4,60000000,512,123456," +
        std::string(precpack::kHasGurobiSupport ? "1" : "0") +
        ",solutions/case.sol\n";
    require(read_text_file(csv_path) == expected_csv,
            "result CSV schema or serialization changed");
    require(read_text_file(assignment_path) ==
                "Bin 1: 1 2\nBin 2: 3\n",
            "assignment file contains redundant metadata or changed format");

    const std::filesystem::path incompatible_path =
        temporary_directory.path() / "incompatible.csv";
    write_text_file(incompatible_path, "legacy_header\n");
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

}

int main() {
    try {
        test_public_defaults();
        test_public_validation();
        test_instance_file_validation();
        test_bpp_profile();
        test_parallel_profile();
        test_salbp_profile();
        test_output_schema();
        std::cout << "Public interface tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Public interface test failed: " << error.what() << '\n';
        return 1;
    }
}
