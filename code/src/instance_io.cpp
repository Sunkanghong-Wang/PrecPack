#include "precpack/instance_io.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace precpack {
namespace {

[[nodiscard]] std::string trim(std::string value) {
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char c) {
        return std::isspace(c) != 0;
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) {
                          return std::isspace(c) != 0;
                      }).base();
    if (first >= last) {
        return {};
    }
    return std::string(first, last);
}

[[nodiscard]] int default_separation(const std::string& problem_type) {
    return problem_type == "SALBP-I" ? 0 : 1;
}

[[nodiscard]] std::vector<int> parse_integers(std::string line) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::replace(line.begin(), line.end(), ';', ' ');
    std::istringstream input(line);
    std::vector<int> values;
    int value = 0;
    while (input >> value) {
        values.push_back(value);
    }
    if (!input.eof()) {
        throw std::invalid_argument("invalid integer row: " + line);
    }
    return values;
}

[[nodiscard]] int parse_integer_scalar(
    const std::string& value,
    std::string_view field,
    const std::filesystem::path& path) {
    int result = 0;
    const auto [end, error] = std::from_chars(
        value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::runtime_error(
            "invalid " + std::string(field) + " in " + path.string());
    }
    return result;
}

[[nodiscard]] double parse_double_scalar(
    const std::string& value,
    std::string_view field,
    const std::filesystem::path& path) {
    std::size_t parsed = 0U;
    double result = 0.0;
    try {
        result = std::stod(value, &parsed);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "invalid " + std::string(field) + " in " + path.string());
    }
    if (parsed != value.size() || !std::isfinite(result)) {
        throw std::runtime_error(
            "invalid " + std::string(field) + " in " + path.string());
    }
    return result;
}

void read_graph_file(const std::filesystem::path& path,
                     const std::string& problem_type,
                     std::vector<Arc>& arcs) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open graph file: " + path.string());
    }
    arcs.clear();
    std::string line;
    bool first_nonempty_line = true;
    std::optional<std::size_t> declared_arc_count;
    while (std::getline(input, line)) {
        line = trim(std::move(line));
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const std::vector<int> values = parse_integers(line);
        if (values.empty()) {
            continue;
        }
        if (first_nonempty_line) {
            if (values.size() != 1U) {
                throw std::runtime_error(
                    "the first graph row must be one arc count: " +
                    path.string());
            }
            if (values.front() < 0) {
                throw std::runtime_error("negative arc count in " +
                                         path.string());
            }
            declared_arc_count = static_cast<std::size_t>(values.front());
            arcs.reserve(*declared_arc_count);
            first_nonempty_line = false;
            continue;
        }
        first_nonempty_line = false;
        if (values.size() != 3U) {
            throw std::runtime_error("invalid graph line in " + path.string() + ": " +
                                     line);
        }
        int separation_value = default_separation(problem_type);
        if (problem_type != "SALBP-I" && problem_type != "BPP-P") {
            separation_value = values[2];
        }
        if (separation_value < 0) {
            throw std::runtime_error("negative separation in " + path.string() +
                                     ": " + line);
        }
        arcs.push_back(Arc{values[0] - 1, values[1] - 1, separation_value});
    }
    if (!declared_arc_count.has_value()) {
        throw std::runtime_error("missing arc count in " + path.string());
    }
    if (arcs.size() != *declared_arc_count) {
        throw std::runtime_error(
            "declared arc count does not match graph contents: " +
            path.string());
    }
}

}

Instance read_instance(const std::filesystem::path& alb_path,
                       const std::optional<std::filesystem::path>& graph_path,
                       const std::string& problem_type,
                       int instance_id) {
    std::ifstream input(alb_path);
    if (!input) {
        throw std::runtime_error("cannot open ALB file: " + alb_path.string());
    }

    Instance instance;
    instance.problem_type = problem_type;
    instance.id = instance_id;
    std::string section;
    std::vector<int> weights;
    std::vector<unsigned char> weight_seen;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(std::move(line));
        if (line.empty()) {
            continue;
        }
        if (line.front() == '<' && line.back() == '>') {
            section = line;
            continue;
        }
        if (line.front() == '#') {
            continue;
        }

        if (section == "<number of tasks>") {
            const int n = parse_integer_scalar(
                line, "task count", alb_path);
            if (n <= 0) {
                throw std::runtime_error("invalid task count in " + alb_path.string());
            }
            weights.assign(static_cast<std::size_t>(n), -1);
            weight_seen.assign(static_cast<std::size_t>(n), 0U);
        } else if (section == "<cycle time>") {
            instance.capacity = parse_integer_scalar(
                line, "cycle time", alb_path);
        } else if (section == "<order strength>") {
            std::replace(line.begin(), line.end(), ',', '.');
            const double value = parse_double_scalar(
                line, "order strength", alb_path);
            constexpr std::array<double, 3> canonical{0.2, 0.6, 0.9};
            const double closest = *std::min_element(
                canonical.begin(), canonical.end(), [&](double lhs, double rhs) {
                    return std::abs(lhs - value) < std::abs(rhs - value);
                });
            std::ostringstream normalized;
            normalized << std::fixed << std::setprecision(1) << closest;
            instance.order_strength = normalized.str();
        } else if (section == "<task times>") {
            const auto values = parse_integers(line);
            if (values.size() != 2U || weights.empty()) {
                throw std::runtime_error("invalid task-time line in " + alb_path.string());
            }
            const int index = values[0] - 1;
            if (index < 0 || index >= static_cast<int>(weights.size())) {
                throw std::runtime_error("task id out of range in " + alb_path.string());
            }
            if (weight_seen[static_cast<std::size_t>(index)] != 0U) {
                throw std::runtime_error(
                    "duplicate task id in " + alb_path.string());
            }
            weight_seen[static_cast<std::size_t>(index)] = 1U;
            weights[static_cast<std::size_t>(index)] = values[1];
        } else if (section == "<precedence relations>") {
            const auto values = parse_integers(line);
            if (values.size() != 2U) {
                throw std::runtime_error("invalid precedence line in " + alb_path.string());
            }
            instance.arcs.push_back(
                Arc{values[0] - 1, values[1] - 1, default_separation(problem_type)});
        }
    }

    if (weights.empty()) {
        throw std::runtime_error("missing <number of tasks> in " + alb_path.string());
    }
    instance.items.reserve(weights.size());
    for (std::size_t i = 0; i < weights.size(); ++i) {
        if (weights[i] <= 0) {
            throw std::runtime_error("missing or invalid weight for task " +
                                     std::to_string(i + 1));
        }
        instance.items.push_back(Item{static_cast<int>(i), weights[i]});
    }

    if (graph_path.has_value()) {
        read_graph_file(*graph_path, problem_type, instance.arcs);
    }
    instance.initialize();
    return instance;
}

}
