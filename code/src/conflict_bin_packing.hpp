#pragma once

#include "precpack/types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace precpack::internal {

[[nodiscard]] bool items_have_same_bin_conflict(
    const Instance& instance,
    int lhs,
    int rhs) noexcept;

struct ConflictBinPackingResult {
    bool attempted = false;
    bool completed = false;
    bool target_test_completed = false;
    bool target_reached = false;
    bool timed_out = false;
    bool node_limited = false;
    bool load_limited = false;
    int lower_bound = 0;
    int optimum = 0;  // Valid only when completed is true.
    std::uint64_t search_nodes = 0;
    std::uint64_t maximal_loads = 0;
    std::uint64_t memo_hits = 0;
    std::uint64_t ordinary_memo_hits = 0;
    double seconds = 0.0;
};

struct ExactRelaxationLookup {
    void* context = nullptr;
    bool (*find)(void*, const std::uint64_t*, int*) = nullptr;
};

class ConflictBinPackingEngine {
public:
    ConflictBinPackingEngine(const Instance& instance,
                             std::uint64_t memo_entry_limit,
                             int maximum_item_count);
    ~ConflictBinPackingEngine();

    ConflictBinPackingEngine(const ConflictBinPackingEngine&) = delete;
    ConflictBinPackingEngine& operator=(
        const ConflictBinPackingEngine&) = delete;
    ConflictBinPackingEngine(ConflictBinPackingEngine&&) noexcept;
    ConflictBinPackingEngine& operator=(
        ConflictBinPackingEngine&&) noexcept;

    [[nodiscard]] ConflictBinPackingResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline,
        double maximum_call_seconds,
        std::uint64_t search_node_limit,
        std::uint64_t maximal_load_limit,
        int ordinary_lower_bound,
        int useful_lower_bound_target,
        ExactRelaxationLookup ordinary_lookup);

    [[nodiscard]] bool has_active_conflict(
        const std::uint64_t* remaining_items) const noexcept;
    [[nodiscard]] int quick_lower_bound(
        const std::uint64_t* remaining_items);
    [[nodiscard]] std::uint64_t conflict_edge_count() const noexcept;
    [[nodiscard]] std::uint64_t memo_entry_count() const noexcept;
    [[nodiscard]] std::uint64_t memory_bytes() const noexcept;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}
