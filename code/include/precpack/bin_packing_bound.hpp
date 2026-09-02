#pragma once

#include "precpack/types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace precpack {

struct BinPackingBoundLimits {
    double call_time_limit_seconds = 1.0;
    std::uint64_t search_node_limit = 250'000U;
    std::uint64_t nondominated_load_limit_per_state = 50U;
    std::uint64_t memo_entry_limit = 200'000U;
    int maximum_item_count = 400;
    bool enable_conflicts = true;
    double conflict_call_time_limit_seconds = 1.0;
    std::uint64_t conflict_search_node_limit = 250'000U;
};

struct BinPackingBoundResult {
    bool attempted = false;
    bool completed = false;
    bool timed_out = false;
    bool node_limited = false;
    bool load_limited = false;
    bool item_limited = false;
    bool target_reached = false;
    bool useful_test_completed = false;
    bool conflict_aware = false;
    bool ordinary_phase_completed = false;
    bool conflict_phase_attempted = false;
    bool conflict_phase_completed = false;
    int lower_bound = 0;
    int optimum = 0;  // Valid only when completed is true.
    std::uint64_t search_nodes = 0;
    std::uint64_t nondominated_loads = 0;
    std::uint64_t memo_hits = 0;
    std::uint64_t conflict_edges = 0;
    std::uint64_t conflict_search_nodes = 0;
    std::uint64_t conflict_maximal_loads = 0;
    std::uint64_t conflict_memo_hits = 0;
    std::uint64_t ordinary_memo_hits = 0;
    double seconds = 0.0;
};

class BinPackingBound {
public:
    BinPackingBound(const Instance& instance, BinPackingBoundLimits limits);
    ~BinPackingBound();

    BinPackingBound(const BinPackingBound&) = delete;
    BinPackingBound& operator=(const BinPackingBound&) = delete;
    BinPackingBound(BinPackingBound&&) noexcept;
    BinPackingBound& operator=(BinPackingBound&&) noexcept;

    [[nodiscard]] BinPackingBoundResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline);
    [[nodiscard]] BinPackingBoundResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline,
        double maximum_call_seconds);
    [[nodiscard]] BinPackingBoundResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline,
        double maximum_call_seconds,
        int useful_lower_bound_target);

    [[nodiscard]] bool lookup_exact(const std::uint64_t* remaining_items,
                                    int* optimum);

    void disable_conflicts() noexcept;

    [[nodiscard]] std::size_t bit_block_count() const noexcept;
    [[nodiscard]] std::uint64_t conflict_edge_count() const noexcept;
    [[nodiscard]] std::uint64_t memo_entry_count() const noexcept;
    [[nodiscard]] std::uint64_t memory_bytes() const noexcept;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}
