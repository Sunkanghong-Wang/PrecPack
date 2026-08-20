#include "precpack/result_io.hpp"

#include "precpack/build_config.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace precpack {
namespace {

inline constexpr const char* kResultHeader =
    "instance_key,problem,instance_file,graph_file,n,capacity,status,"
    "lower_bound,upper_bound,gap,time_seconds,time_limit_seconds,threads,"
    "state_limit,memory_limit_mb,bbr_peak_memory_bytes,gurobi_enabled,"
    "solution_file";

[[nodiscard]] std::string csv_field(std::string_view value) {
    if (value.find_first_of(",\"\r\n") == std::string_view::npos) {
        return std::string(value);
    }
    std::string escaped = "\"";
    for (const char character : value) {
        if (character == '\"') {
            escaped += '\"';
        }
        escaped += character;
    }
    escaped += '\"';
    return escaped;
}

[[nodiscard]] std::string recorded_path(const std::filesystem::path& path) {
    if (path.empty()) {
        return {};
    }
    const std::filesystem::path normalized = path.lexically_normal();
    if (!normalized.is_absolute()) {
        return normalized.generic_string();
    }

    std::error_code error;
    const std::filesystem::path current =
        std::filesystem::current_path(error).lexically_normal();
    if (!error) {
        const std::filesystem::path relative =
            normalized.lexically_relative(current);
        const auto first = relative.begin();
        if (!relative.empty() && first != relative.end() && *first != "..") {
            return relative.generic_string();
        }
    }
    return normalized.generic_string();
}

[[nodiscard]] std::string sanitize_filename(std::string value) {
    for (char& character : value) {
        const bool alphanumeric =
            (character >= '0' && character <= '9') ||
            (character >= 'A' && character <= 'Z') ||
            (character >= 'a' && character <= 'z');
        if (!alphanumeric && character != '-' && character != '_' &&
            character != '.') {
            character = '_';
        }
    }
    return value.empty() ? "instance" : value;
}

void update_path_hash(std::uint64_t& hash, std::string_view value) noexcept {
    constexpr std::uint64_t kFnvPrime = 1'099'511'628'211ULL;
    for (const unsigned char byte : value) {
        hash ^= byte;
        hash *= kFnvPrime;
    }
}

void require_compatible_header(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot read result CSV: " + path.string());
    }
    std::string header;
    std::getline(input, header);
    if (!header.empty() && header.back() == '\r') {
        header.pop_back();
    }
    if (header != kResultHeader) {
        throw std::runtime_error(
            "existing result CSV has an incompatible header: " +
            path.string());
    }
}

}

std::string make_instance_key(
    const std::filesystem::path& instance_path,
    const std::optional<std::filesystem::path>& graph_path) {
    constexpr std::uint64_t kFnvOffsetBasis = 14'695'981'039'346'656'037ULL;
    std::uint64_t hash = kFnvOffsetBasis;
    update_path_hash(hash, recorded_path(instance_path));
    update_path_hash(hash, std::string_view{"\0", 1});
    if (graph_path.has_value()) {
        update_path_hash(hash, recorded_path(*graph_path));
    }

    std::ostringstream key;
    key << sanitize_filename(instance_path.stem().string()) << "__"
        << std::hex << std::setfill('0') << std::setw(16) << hash;
    return key.str();
}

void append_result_csv(const std::filesystem::path& path,
                       std::string_view instance_key,
                       const std::filesystem::path& instance_path,
                       const std::optional<std::filesystem::path>& graph_path,
                       const Instance& instance,
                       const Solution& solution,
                       const std::filesystem::path& assignment_path) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const bool needs_header = !std::filesystem::exists(path) ||
                              std::filesystem::file_size(path) == 0U;
    if (!needs_header) {
        require_compatible_header(path);
    }
    std::ofstream output(path, std::ios::app | std::ios::binary);
    if (!output) {
        throw std::runtime_error("cannot open result CSV: " + path.string());
    }
    if (needs_header) {
        output << kResultHeader << '\n';
    }

    constexpr std::uint64_t kMegabyte = 1024ULL * 1024ULL;
    const std::uint64_t memory_limit_mb =
        solution.bbr_stats.memory_limit_bytes / kMegabyte;
    output << std::setprecision(12) << csv_field(instance_key) << ','
           << csv_field(instance.problem_type) << ','
           << csv_field(recorded_path(instance_path)) << ','
           << csv_field(graph_path.has_value()
                            ? recorded_path(*graph_path)
                            : std::string{})
           << ',' << instance.size() << ',' << instance.capacity << ','
           << to_string(solution.status) << ',' << solution.lower_bound << ','
           << solution.upper_bound << ',' << solution.relative_gap << ','
           << solution.stats.total_seconds << ','
           << solution.bbr_stats.time_limit_seconds << ',' << solution.threads
           << ',' << solution.bbr_stats.configured_state_limit << ','
           << memory_limit_mb << ','
           << solution.bbr_stats.peak_memory_bytes << ','
           << (kHasGurobiSupport ? 1 : 0) << ','
           << csv_field(assignment_path.generic_string()) << '\n';
}

void write_assignment(const std::filesystem::path& path,
                      const Instance& instance,
                      const Solution& solution) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("cannot open assignment file: " +
                                 path.string());
    }

    for (int bin = 0; bin < solution.assignment.bin_count; ++bin) {
        output << "Bin " << bin + 1 << ':';
        for (int item = 0; item < instance.size(); ++item) {
            if (solution.assignment.bin_of_item[static_cast<std::size_t>(item)] !=
                bin) {
                continue;
            }
            const Item& data = instance.items[static_cast<std::size_t>(item)];
            output << ' ' << data.original_index + 1;
        }
        output << '\n';
    }
}

}
