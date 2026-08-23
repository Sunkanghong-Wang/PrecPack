#include "precpack/result_io.hpp"

#include "precpack/build_config.hpp"

#include "environment.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
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
    "instance_key,problem,instance_file,graph_file,n,capacity,status,"
    "lower_bound,upper_bound,gap,time_seconds,time_limit_seconds,threads,"
    "memory_limit_mb,bbr_peak_memory_bytes,gurobi_enabled,"
    "gurobi_required,solution_file";
inline constexpr std::size_t kResultColumnCount = 18U;

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

[[nodiscard]] std::string recorded_path(const std::filesystem::path& path) {
    if (path.empty()) {
        return {};
    }
    const std::filesystem::path normalized = path.lexically_normal();
    if (!normalized.is_absolute()) {
        return normalized.generic_string();
    }

    std::error_code error;
    std::filesystem::path current;
    const std::optional<std::string> repository =
        internal::environment_value("PRECPACK_REPOSITORY_ROOT");
    if (repository.has_value()) {
        current =
            std::filesystem::absolute(*repository, error).lexically_normal();
    } else {
        current = std::filesystem::current_path(error).lexically_normal();
    }
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

bool read_csv_record(std::istream& input, std::vector<std::string>& fields) {
    fields.clear();
    std::string field;
    bool quoted = false;
    bool quote_closed = false;
    bool consumed = false;
    char character = '\0';
    while (input.get(character)) {
        consumed = true;
        if (quoted) {
            if (character != '"') {
                field += character;
                continue;
            }
            if (input.peek() == '"') {
                static_cast<void>(input.get(character));
                field += '"';
                continue;
            }
            quoted = false;
            quote_closed = true;
            continue;
        }
        if (character == ',') {
            fields.push_back(std::move(field));
            field.clear();
            quote_closed = false;
        } else if (character == '\n' || character == '\r') {
            if (character == '\r' && input.peek() == '\n') {
                static_cast<void>(input.get(character));
            }
            fields.push_back(std::move(field));
            return true;
        } else if (character == '"') {
            if (!field.empty() || quote_closed) {
                throw std::runtime_error("malformed quoted CSV field");
            }
            quoted = true;
        } else {
            if (quote_closed) {
                throw std::runtime_error(
                    "unexpected characters after a quoted CSV field");
            }
            field += character;
        }
    }
    if (quoted) {
        throw std::runtime_error("unterminated quoted CSV field");
    }
    if (!consumed && field.empty() && fields.empty()) {
        return false;
    }
    fields.push_back(std::move(field));
    return true;
}

[[nodiscard]] std::uint64_t parse_unsigned_csv_field(
    const std::string& value,
    std::string_view column,
    std::size_t row,
    const std::filesystem::path& path) {
    std::size_t parsed = 0U;
    std::uint64_t result = 0U;
    bool converted = true;
    try {
        result = std::stoull(value, &parsed);
    } catch (const std::exception&) {
        converted = false;
    }
    if (!converted || value.empty() || value.front() == '-' ||
        value.front() == '+' || parsed != value.size()) {
        throw std::runtime_error(
            "invalid " + std::string(column) + " in result CSV row " +
            std::to_string(row) + ": " + path.string());
    }
    return result;
}

[[nodiscard]] double parse_double_csv_field(
    const std::string& value,
    std::string_view column,
    std::size_t row,
    const std::filesystem::path& path) {
    std::size_t parsed = 0U;
    double result = 0.0;
    bool converted = true;
    try {
        result = std::stod(value, &parsed);
    } catch (const std::exception&) {
        converted = false;
    }
    if (!converted || value.empty() || parsed != value.size() ||
        !std::isfinite(result)) {
        throw std::runtime_error(
            "invalid " + std::string(column) + " in result CSV row " +
            std::to_string(row) + ": " + path.string());
    }
    return result;
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
    std::ostringstream record;
    if (needs_header) {
        record << kResultHeader << '\n';
    }

    constexpr std::uint64_t kMegabyte = 1024ULL * 1024ULL;
    const std::uint64_t memory_limit_mb =
        solution.bbr_stats.memory_limit_bytes / kMegabyte;
    record << std::setprecision(12) << csv_field(instance_key) << ','
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
           << ',' << memory_limit_mb << ','
           << solution.bbr_stats.peak_memory_bytes << ','
           << (kHasGurobiSupport ? 1 : 0) << ','
           << (solution.gurobi_runtime_required ? 1 : 0) << ','
           << csv_field(assignment_path.generic_string()) << '\n';

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

std::vector<ResultReference> read_result_references(
    const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        return {};
    }
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("result CSV is not a file: " + path.string());
    }
    if (std::filesystem::file_size(path) == 0U) {
        return {};
    }
    require_compatible_header(path);

    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot read result CSV: " + path.string());
    }
    std::string header;
    std::getline(input, header);

    std::vector<ResultReference> references;
    std::vector<std::string> fields;
    std::size_t row = 1U;
    while (read_csv_record(input, fields)) {
        ++row;
        if (fields.size() != kResultColumnCount || fields.front().empty() ||
            fields.back().empty()) {
            throw std::runtime_error(
                "malformed result CSV row " + std::to_string(row) + ": " +
                path.string());
        }
        const std::uint64_t threads = parse_unsigned_csv_field(
            fields[12], "threads", row, path);
        if (threads == 0U ||
            threads > static_cast<std::uint64_t>(
                          std::numeric_limits<int>::max())) {
            throw std::runtime_error(
                "invalid threads in result CSV row " + std::to_string(row) +
                ": " + path.string());
        }
        const std::uint64_t gurobi = parse_unsigned_csv_field(
            fields[15], "gurobi_enabled", row, path);
        if (gurobi > 1U) {
            throw std::runtime_error(
                "invalid gurobi_enabled in result CSV row " +
                std::to_string(row) + ": " + path.string());
        }
        const std::uint64_t gurobi_required = parse_unsigned_csv_field(
            fields[16], "gurobi_required", row, path);
        if (gurobi_required > 1U) {
            throw std::runtime_error(
                "invalid gurobi_required in result CSV row " +
                std::to_string(row) + ": " + path.string());
        }
        references.push_back(ResultReference{
            fields[0],
            fields[1],
            parse_double_csv_field(fields[11], "time_limit_seconds", row,
                                   path),
            static_cast<int>(threads),
            parse_unsigned_csv_field(fields[13], "memory_limit_mb", row,
                                     path),
            gurobi == 1U,
            gurobi_required == 1U,
            fields[17],
        });
    }
    return references;
}

void require_unused_instance_key(const std::filesystem::path& path,
                                 std::string_view instance_key) {
    for (const ResultReference& reference : read_result_references(path)) {
        if (reference.instance_key == instance_key) {
            throw std::runtime_error(
                "result CSV already contains instance_key=" +
                std::string(instance_key) +
                "; use a different output directory or batch mode");
        }
    }
}

}
