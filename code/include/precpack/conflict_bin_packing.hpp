#pragma once

#include "precpack/types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace precpack {

class BinPackingBound;

struct ConflictBinPackingLimits {
    double call_time_limit_seconds = 1.0;
    double total_time_limit_seconds = 1.0;
    std::uint64_t search_node_limit = 250'000U;
    std::uint64_t maximal_load_limit_per_state = 50U;
    std::uint64_t memo_entry_limit = 200'000U;
    int maximum_item_count = 400;
};

struct ConflictBinPackingResult {
    bool attempted = false;
    bool completed = false;
    bool timed_out = false;
    bool node_limited = false;
    bool load_limited = false;
    bool item_limited = false;
    bool total_budget_exhausted = false;
    int lower_bound = 0;
    int optimum = 0;  // Valid only when completed is true.
    std::uint64_t search_nodes = 0;
    std::uint64_t maximal_loads = 0;
    std::uint64_t memo_hits = 0;
    std::uint64_t ordinary_memo_hits = 0;
    double seconds = 0.0;
};

class ConflictBinPackingBound {
public:
    ConflictBinPackingBound(const Instance& instance,
                            ConflictBinPackingLimits limits,
                            BinPackingBound* ordinary_bound = nullptr);
    ~ConflictBinPackingBound();

    ConflictBinPackingBound(const ConflictBinPackingBound&) = delete;
    ConflictBinPackingBound& operator=(const ConflictBinPackingBound&) = delete;
    ConflictBinPackingBound(ConflictBinPackingBound&&) noexcept;
    ConflictBinPackingBound& operator=(ConflictBinPackingBound&&) noexcept;

    [[nodiscard]] ConflictBinPackingResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline);

    [[nodiscard]] std::size_t bit_block_count() const noexcept;
    [[nodiscard]] std::uint64_t conflict_edge_count() const noexcept;
    [[nodiscard]] std::uint64_t memo_entry_count() const noexcept;
    [[nodiscard]] std::uint64_t memory_bytes() const noexcept;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}
