#include "precpack/result_io.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace precpack {
namespace {

inline constexpr const char* kResultHeader =
    "instance_set,instance,n,time_limit_seconds,memory_limit_mb,status,opt,"
    "lower_bound,upper_bound,time_seconds,"
    "bbr_peak_memory_bytes,bbr_states_created";

class TemporaryFileGuard {
public:
    explicit TemporaryFileGuard(std::filesystem::path path)
        : path_(std::move(path)) {}

    ~TemporaryFileGuard() {
        if (path_.empty()) {
            return;
        }
        std::error_code error;
        std::filesystem::remove(path_, error);
    }

    void release() noexcept {
        path_.clear();
    }

private:
    std::filesystem::path path_;
};

void finish_output(std::ofstream& output,
                   const std::filesystem::path& path,
                   std::string_view description) {
    output.flush();
    if (!output) {
        throw std::runtime_error(
            "cannot flush " + std::string(description) + ": " +
            path.string());
    }
    output.close();
    if (!output) {
        throw std::runtime_error(
            "cannot close " + std::string(description) + ": " +
            path.string());
    }
}

void replace_file(const std::filesystem::path& source,
                  const std::filesystem::path& destination) {
#if defined(_WIN32)
    if (MoveFileExW(source.c_str(), destination.c_str(),
                    MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) == 0) {
        throw std::system_error(
            static_cast<int>(GetLastError()), std::system_category(),
            "cannot replace assignment file: " + destination.string());
    }
#else
    std::error_code error;
    std::filesystem::rename(source, destination, error);
    if (error) {
        throw std::system_error(
            error, "cannot replace assignment file: " +
                       destination.string());
    }
#endif
}

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

[[nodiscard]] std::filesystem::path parent_below(
    const std::filesystem::path& file,
    std::string_view anchor) {
    const std::filesystem::path parent =
        file.lexically_normal().parent_path();
    bool below_anchor = false;
    std::filesystem::path result;
    for (const std::filesystem::path& component : parent) {
        if (below_anchor) {
            result /= component;
        } else if (component == std::filesystem::path(anchor)) {
            below_anchor = true;
        }
    }
    return below_anchor ? result : std::filesystem::path{};
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

std::string make_instance_set(
    const std::filesystem::path& instance_path,
    const std::optional<std::filesystem::path>& graph_path) {
    std::filesystem::path result;
    if (graph_path.has_value()) {
        result = parent_below(*graph_path, "bpp-gp-graphs");
        if (result.empty()) {
            result = graph_path->parent_path().filename();
        }
    } else {
        result = parent_below(instance_path, "instances");
        if (result.empty()) {
            result = instance_path.parent_path().filename();
        }
    }
    return result.empty() ? "external" : result.generic_string();
}

void append_result_csv(const std::filesystem::path& path,
                       std::string_view instance_set,
                       std::string_view instance_name,
                       const Instance& instance,
                       const Solution& solution) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    validate_result_csv(path);
    const bool needs_header = !std::filesystem::exists(path) ||
                              std::filesystem::file_size(path) == 0U;
    std::ostringstream record;
    if (needs_header) {
        record << kResultHeader << '\n';
    }

    constexpr std::uint64_t kMegabyte = 1024ULL * 1024ULL;
    const std::uint64_t memory_limit_mb =
        solution.bbr_stats.memory_limit_bytes / kMegabyte;
    record << std::setprecision(12) << csv_field(instance_set) << ','
           << csv_field(instance_name) << ',' << instance.size() << ','
           << solution.bbr_stats.time_limit_seconds << ',' << memory_limit_mb
           << ',' << to_string(solution.status)
           << ',' << (solution.status == SolveStatus::kOptimal ? 1 : 0) << ','
           << solution.lower_bound << ',' << solution.upper_bound << ','
           << solution.stats.total_seconds << ','
           << solution.bbr_stats.peak_memory_bytes << ','
           << solution.bbr_stats.states_created << '\n';

    const std::string serialized = record.str();
    std::ofstream output(path, std::ios::app | std::ios::binary);
    if (!output) {
        throw std::runtime_error("cannot open result CSV: " + path.string());
    }
    output.write(serialized.data(),
                 static_cast<std::streamsize>(serialized.size()));
    if (!output) {
        throw std::runtime_error("cannot write result CSV: " + path.string());
    }
    finish_output(output, path, "result CSV");
}

void validate_result_csv(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        return;
    }
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("result CSV is not a file: " + path.string());
    }
    if (std::filesystem::file_size(path) > 0U) {
        require_compatible_header(path);
    }
}

void write_assignment(const std::filesystem::path& path,
                      const Instance& instance,
                      const Solution& solution) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::filesystem::path temporary_path = path;
    temporary_path += ".tmp";
    TemporaryFileGuard temporary_file(temporary_path);
    {
        std::ofstream output(
            temporary_path, std::ios::binary | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("cannot open assignment file: " +
                                     temporary_path.string());
        }

        for (int bin = 0; bin < solution.assignment.bin_count; ++bin) {
            output << "Bin " << bin + 1 << ':';
            for (int item = 0; item < instance.size(); ++item) {
                if (solution.assignment.bin_of_item[
                        static_cast<std::size_t>(item)] != bin) {
                    continue;
                }
                const Item& data =
                    instance.items[static_cast<std::size_t>(item)];
                output << ' ' << data.original_index + 1;
            }
            output << '\n';
        }
        if (!output) {
            throw std::runtime_error(
                "cannot write assignment file: " + temporary_path.string());
        }
        finish_output(output, temporary_path, "assignment file");
    }
    replace_file(temporary_path, path);
    temporary_file.release();
}

bool has_nonempty_solution(const std::filesystem::path& path) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error) || error) {
        return false;
    }
    return std::filesystem::file_size(path, error) > 0U && !error;
}

}
