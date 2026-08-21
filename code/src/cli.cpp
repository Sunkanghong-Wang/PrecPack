#include "precpack/cli.hpp"

#include "precpack/build_config.hpp"

#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>

namespace precpack {
namespace {

[[noreturn]] void usage_error(const std::string& message) {
    throw std::invalid_argument(message + "\nUse --help for command syntax.");
}

[[nodiscard]] std::string require_value(int& index, int argc,
                                        char* const argv[]) {
    if (index + 1 >= argc) {
        usage_error(std::string("missing value after ") + argv[index]);
    }
    return argv[++index];
}

[[nodiscard]] double parse_positive_double(const std::string& text,
                                           std::string_view option) {
    std::size_t parsed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &parsed);
    } catch (const std::exception&) {
        usage_error(std::string(option) +
                    " requires a positive number: " + text);
    }
    if (parsed != text.size() || !std::isfinite(value) || value <= 0.0) {
        usage_error(std::string(option) +
                    " requires a positive number: " + text);
    }
    return value;
}

[[nodiscard]] std::uint64_t parse_positive_integer(const std::string& text,
                                                   std::string_view option) {
    if (text.empty() || text.front() == '-' || text.front() == '+') {
        usage_error(std::string(option) +
                    " requires a positive integer: " + text);
    }
    std::size_t parsed = 0;
    std::uint64_t value = 0;
    try {
        value = std::stoull(text, &parsed);
    } catch (const std::exception&) {
        usage_error(std::string(option) +
                    " requires a positive integer: " + text);
    }
    if (parsed != text.size() || value == 0U ||
        value >
            std::numeric_limits<std::uint64_t>::max() / (1024ULL * 1024ULL)) {
        usage_error(std::string(option) +
                    " requires a positive integer: " + text);
    }
    return value;
}

[[nodiscard]] int parse_threads(const std::string& text) {
    if (text == "-1") {
        return -1;
    }
    const std::uint64_t value = parse_positive_integer(text, "--threads");
    if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        usage_error("--threads is too large: " + text);
    }
    return static_cast<int>(value);
}

}

CommandLineOptions parse_command_line(int argc, char* const argv[]) {
    CommandLineOptions options;
    bool problem_was_set = false;
    bool output_was_set = false;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            options.show_help = true;
            return options;
        }
        if (argument == "--batch") {
            options.batch_mode = true;
        } else if (argument == "--check-only") {
            options.check_only = true;
        } else if (argument == "--problem") {
            options.problem =
                parse_problem_kind(require_value(index, argc, argv));
            problem_was_set = true;
        } else if (argument == "--instance") {
            options.instance_path = require_value(index, argc, argv);
        } else if (argument == "--graph") {
            options.graph_path = require_value(index, argc, argv);
        } else if (argument == "--input") {
            options.input_path = require_value(index, argc, argv);
        } else if (argument == "--graph-dir") {
            options.graph_directory = require_value(index, argc, argv);
        } else if (argument == "--time-limit") {
            options.time_limit_seconds = parse_positive_double(
                require_value(index, argc, argv), "--time-limit");
        } else if (argument == "--memory-limit-mb") {
            options.memory_limit_mb = parse_positive_integer(
                require_value(index, argc, argv), "--memory-limit-mb");
        } else if (argument == "--threads") {
            options.threads = parse_threads(require_value(index, argc, argv));
        } else if (argument == "--output-dir") {
            options.output_directory = require_value(index, argc, argv);
            output_was_set = true;
        } else {
            usage_error("unknown argument: " + argument);
        }
    }

    if (!problem_was_set) {
        usage_error("--problem is required");
    }
    if (options.output_directory.empty()) {
        usage_error("--output-dir cannot be empty");
    }
    if (options.batch_mode) {
        if (!options.instance_path.empty() || options.graph_path.has_value()) {
            usage_error("--instance and --graph are not valid with --batch");
        }
        if (options.input_path.has_value() && options.input_path->empty()) {
            usage_error("--input cannot be empty");
        }
        if (options.graph_directory.has_value() &&
            options.graph_directory->empty()) {
            usage_error("--graph-dir cannot be empty");
        }
        if (options.problem != ProblemKind::kBppGp &&
            options.graph_directory.has_value()) {
            usage_error("--graph-dir is only valid for bpp-gp");
        }
        if (!output_was_set) {
            options.output_directory =
                std::filesystem::path("results") / to_slug(options.problem);
        }
    } else {
        if (options.input_path.has_value() ||
            options.graph_directory.has_value()) {
            usage_error("--input and --graph-dir require --batch");
        }
        if (options.check_only) {
            usage_error("--check-only requires --batch");
        }
        if (options.instance_path.empty()) {
            usage_error("--instance is required");
        }
        if (options.problem == ProblemKind::kBppGp &&
            !options.graph_path.has_value()) {
            usage_error("--graph is required for bpp-gp");
        }
        if (options.graph_path.has_value() && options.graph_path->empty()) {
            usage_error("--graph cannot be empty");
        }
        if (options.problem != ProblemKind::kBppGp &&
            options.graph_path.has_value()) {
            usage_error("--graph is only valid for bpp-gp");
        }
    }
    return options;
}

void print_help(std::ostream& output, std::string_view executable) {
    output
        << "PrecPack exact solver\n\n"
        << "Usage:\n  " << executable
        << " --problem TYPE --instance FILE [options]\n  " << executable
        << " --batch --problem TYPE [batch options]\n\n"
        << "Required:\n"
        << "  --problem TYPE          salbp-i, bpp-p, or bpp-gp\n\n"
        << "Single-instance input:\n"
        << "  --instance FILE         ALB-format .txt instance file\n"
        << "  --graph FILE            Labeled .graph file (bpp-gp only)\n\n"
        << "Batch options:\n"
        << "  --batch                 Run selected instances sequentially\n"
        << "  --input PATH            .txt file or directory; bundled data by default\n"
        << "  --graph-dir PATH        .graph file or directory (bpp-gp only)\n"
        << "  --check-only            Validate selection and existing output only\n\n"
        << "Resource and output options:\n"
        << "  --time-limit SECONDS    Global limit (default: 300)\n"
        << "  --memory-limit-mb MB    BBR memory cap (default: 24576)\n"
        << "  --threads N             BBR workers; -1 uses available CPUs "
           "(default: 1)\n"
        << "  --output-dir DIR        Results directory (default: results;\n"
        << "                          results/TYPE in batch mode)\n"
        << "  --help                  Show this message\n\n"
        << "Build capabilities:\n"
        << "  Optional Gurobi root strengthening: "
        << (kHasGurobiSupport ? "enabled" : "disabled") << "\n\n"
        << "Algorithm choices are fixed by --problem. Parallelism is optional "
           "and applies only to exact BBR.\n";
}

}
