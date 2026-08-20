#pragma once

#include "precpack/types.hpp"

#include <cstdint>
#include <string_view>

namespace precpack {

enum class ProblemKind {
    kSalbpI,
    kBppP,
    kBppGp,
};

[[nodiscard]] ProblemKind parse_problem_kind(std::string_view value);
[[nodiscard]] const char* to_string(ProblemKind problem) noexcept;

[[nodiscard]] Config make_solver_config(
    ProblemKind problem,
    double time_limit_seconds,
    std::uint64_t memory_limit_mb,
    int threads = 1);

[[nodiscard]] int resolve_thread_count(int requested_threads) noexcept;

}
