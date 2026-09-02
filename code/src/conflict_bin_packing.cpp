#include "conflict_bin_packing.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace precpack::internal {
namespace {

using Clock = std::chrono::steady_clock;

constexpr std::size_t kMaximumConflictBlocks = 16U;
constexpr int kMaximumConflictItems =
    static_cast<int>(64U * kMaximumConflictBlocks);
using ConflictMask = std::array<std::uint64_t, kMaximumConflictBlocks>;

[[nodiscard]] std::uint64_t mix64(std::uint64_t value) noexcept {
    value ^= value >> 30U;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27U;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return value;
}

[[nodiscard]] std::uint64_t hash_words(const std::uint64_t* words,
                                       std::size_t count) noexcept {
    std::uint64_t hash = 0x6a09e667f3bcc909ULL ^
                         mix64(static_cast<std::uint64_t>(count));
    if (count == 1U) {
        return mix64(hash ^ mix64(words[0] + 0x9e3779b97f4a7c15ULL));
    }
    if (count == 2U) {
        hash = mix64(hash ^ mix64(words[0] + 0x9e3779b97f4a7c15ULL));
        return mix64(hash ^
                     mix64(words[1] + 2U * 0x9e3779b97f4a7c15ULL));
    }
    for (std::size_t index = 0; index < count; ++index) {
        hash = mix64(hash ^ mix64(
            words[index] + 0x9e3779b97f4a7c15ULL * (index + 1U)));
    }
    return hash;
}

[[nodiscard]] bool bit_is_set(const std::uint64_t* words,
                              int item) noexcept {
    return ((words[static_cast<std::size_t>(item) >> 6U] >>
             (static_cast<unsigned>(item) & 63U)) &
            std::uint64_t{1}) != 0U;
}

void set_bit(std::uint64_t* words, int item) noexcept {
    words[static_cast<std::size_t>(item) >> 6U] |=
        std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
}

void clear_bit(std::uint64_t* words, int item) noexcept {
    words[static_cast<std::size_t>(item) >> 6U] &=
        ~(std::uint64_t{1} << (static_cast<unsigned>(item) & 63U));
}

[[nodiscard]] int ceil_div(std::int64_t numerator,
                           int denominator) noexcept {
    return static_cast<int>((numerator + denominator - 1) / denominator);
}

template <class T>
[[nodiscard]] std::uint64_t vector_memory_bytes(
    const std::vector<T>& values) noexcept {
    const std::uint64_t count = values.capacity();
    if (count > std::numeric_limits<std::uint64_t>::max() / sizeof(T)) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    return count * sizeof(T);
}

[[nodiscard]] std::uint64_t saturated_add(std::uint64_t lhs,
                                          std::uint64_t rhs) noexcept {
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    return lhs > maximum - rhs ? maximum : lhs + rhs;
}

class ExactSubsetMemo {
public:
    ExactSubsetMemo(std::size_t blocks, std::uint64_t maximum_entries)
        : blocks_(blocks), maximum_entries_(maximum_entries) {
        if (maximum_entries_ == 0U) {
            return;
        }
        std::size_t slot_capacity = 1024U;
        while (maximum_entries_ * 10U >=
               static_cast<std::uint64_t>(slot_capacity) * 7U) {
            if (slot_capacity >
                std::numeric_limits<std::size_t>::max() / 2U) {
                throw std::length_error(
                    "conflict-aware BINLB hash table is too large");
            }
            slot_capacity *= 2U;
        }
        slots_.assign(slot_capacity, 0U);
        const std::size_t entry_capacity =
            static_cast<std::size_t>(maximum_entries_);
        if (blocks_ != 0U &&
            entry_capacity > keys_.max_size() / blocks_) {
            throw std::length_error(
                "conflict-aware BINLB memo is too large");
        }
        keys_.reserve(entry_capacity * blocks_);
        hashes_.reserve(entry_capacity);
        lower_bounds_.reserve(entry_capacity);
        upper_bounds_.reserve(entry_capacity);
    }

    [[nodiscard]] bool find(const std::uint64_t* key,
                            int* value) const noexcept {
        if (slots_.empty()) {
            return false;
        }
        const std::size_t entry = find_entry(
            key, hash_words(key, blocks_));
        if (entry == kMissing ||
            lower_bounds_[entry] != upper_bounds_[entry]) {
            return false;
        }
        *value = lower_bounds_[entry];
        return true;
    }

    void insert(const std::uint64_t* key, int value) {
        if (value < 0 || value >=
                static_cast<int>(kUnknownUpper)) {
            throw std::invalid_argument("invalid exact subset memo value");
        }
        const std::uint64_t hash = hash_words(key, blocks_);
        const std::size_t entry = ensure_entry(key, hash);
        if (entry == kMissing) {
            return;
        }
        lower_bounds_[entry] = static_cast<std::uint16_t>(value);
        upper_bounds_[entry] = static_cast<std::uint16_t>(value);
    }

    [[nodiscard]] int query_target(const std::uint64_t* key,
                                   int bin_limit) const noexcept {
        if (slots_.empty()) {
            return 0;
        }
        const std::size_t entry = find_entry(
            key, hash_words(key, blocks_));
        if (entry == kMissing) {
            return 0;
        }
        if (lower_bounds_[entry] > bin_limit) {
            return -1;
        }
        if (upper_bounds_[entry] <= bin_limit) {
            return 1;
        }
        return 0;
    }

    void record_feasible(const std::uint64_t* key, int bin_limit) {
        if (bin_limit < 0 || bin_limit >=
                static_cast<int>(kUnknownUpper)) {
            return;
        }
        const std::size_t entry = ensure_entry(
            key, hash_words(key, blocks_));
        if (entry == kMissing) {
            return;
        }
        upper_bounds_[entry] = std::min(
            upper_bounds_[entry], static_cast<std::uint16_t>(bin_limit));
        if (lower_bounds_[entry] > upper_bounds_[entry]) {
            throw std::logic_error(
                "inconsistent feasible conflict BINLB memo bound");
        }
    }

    void record_infeasible(const std::uint64_t* key, int bin_limit) {
        if (bin_limit < 0 || bin_limit + 1 >=
                static_cast<int>(kUnknownUpper)) {
            return;
        }
        const std::size_t entry = ensure_entry(
            key, hash_words(key, blocks_));
        if (entry == kMissing) {
            return;
        }
        lower_bounds_[entry] = std::max(
            lower_bounds_[entry],
            static_cast<std::uint16_t>(bin_limit + 1));
        if (lower_bounds_[entry] > upper_bounds_[entry]) {
            throw std::logic_error(
                "inconsistent infeasible conflict BINLB memo bound");
        }
    }

    [[nodiscard]] std::uint64_t size() const noexcept {
        return lower_bounds_.size();
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        std::uint64_t bytes = vector_memory_bytes(slots_);
        bytes = saturated_add(bytes, vector_memory_bytes(keys_));
        bytes = saturated_add(bytes, vector_memory_bytes(hashes_));
        bytes = saturated_add(bytes, vector_memory_bytes(lower_bounds_));
        return saturated_add(bytes, vector_memory_bytes(upper_bounds_));
    }

private:
    static constexpr std::size_t kMissing =
        std::numeric_limits<std::size_t>::max();
    static constexpr std::uint16_t kUnknownUpper =
        std::numeric_limits<std::uint16_t>::max();

    [[nodiscard]] std::size_t find_entry(
        const std::uint64_t* key,
        std::uint64_t hash) const noexcept {
        std::size_t slot = static_cast<std::size_t>(hash) &
                           (slots_.size() - 1U);
        while (slots_[slot] != 0U) {
            const std::size_t entry = slots_[slot] - 1U;
            if (hashes_[entry] == hash && keys_equal(entry, key)) {
                return entry;
            }
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
        return kMissing;
    }

    [[nodiscard]] std::size_t ensure_entry(
        const std::uint64_t* key,
        std::uint64_t hash) {
        if (maximum_entries_ == 0U) {
            return kMissing;
        }
        const std::size_t existing = find_entry(key, hash);
        if (existing != kMissing) {
            return existing;
        }
        if (lower_bounds_.size() >= maximum_entries_) {
            return kMissing;
        }
        const std::size_t entry = lower_bounds_.size();
        keys_.insert(keys_.end(), key, key + blocks_);
        hashes_.push_back(hash);
        lower_bounds_.push_back(0U);
        upper_bounds_.push_back(kUnknownUpper);
        std::size_t slot = static_cast<std::size_t>(hash) &
                           (slots_.size() - 1U);
        while (slots_[slot] != 0U) {
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
        slots_[slot] = static_cast<std::uint32_t>(entry + 1U);
        return entry;
    }

    [[nodiscard]] bool keys_equal(std::size_t entry,
                                  const std::uint64_t* key) const noexcept {
        const std::uint64_t* stored = keys_.data() + entry * blocks_;
        if (blocks_ == 1U) {
            return stored[0] == key[0];
        }
        if (blocks_ == 2U) {
            return stored[0] == key[0] && stored[1] == key[1];
        }
        return std::equal(stored, stored + blocks_, key);
    }

    std::size_t blocks_ = 0;
    std::uint64_t maximum_entries_ = 0;
    std::vector<std::uint32_t> slots_;
    std::vector<std::uint64_t> keys_;
    std::vector<std::uint64_t> hashes_;
    std::vector<std::uint16_t> lower_bounds_;
    std::vector<std::uint16_t> upper_bounds_;
};

struct NodeResult {
    bool completed = false;
    int lower_bound = 0;
    int optimum = 0;
};

enum class TargetDecision : unsigned char {
    kUnknown,
    kFeasible,
    kInfeasible,
};

}

bool items_have_same_bin_conflict(const Instance& instance,
                                  int lhs,
                                  int rhs) noexcept {
    if (instance.has_positive_separation_path(lhs, rhs) ||
        instance.has_positive_separation_path(rhs, lhs)) {
        return true;
    }

    int from = -1;
    int to = -1;
    if (instance.reaches(lhs, rhs)) {
        from = lhs;
        to = rhs;
    } else if (instance.reaches(rhs, lhs)) {
        from = rhs;
        to = lhs;
    } else {
        return false;
    }

    std::int64_t interval_weight =
        static_cast<std::int64_t>(
            instance.items[static_cast<std::size_t>(from)].weight) +
        instance.items[static_cast<std::size_t>(to)].weight;
    if (interval_weight > instance.capacity) {
        return true;
    }

    const std::uint64_t* descendants = instance.reachable_row(from);
    const std::uint64_t* ancestors = instance.reaching_row(to);
    for (std::size_t block = 0; block < instance.reachability_blocks; ++block) {
        std::uint64_t value = descendants[block] & ancestors[block];
        while (value != 0U) {
            const unsigned bit = std::countr_zero(value);
            const int item = static_cast<int>(block * 64U + bit);
            interval_weight +=
                instance.items[static_cast<std::size_t>(item)].weight;
            if (interval_weight > instance.capacity) {
                return true;
            }
            value &= value - 1U;
        }
    }
    return false;
}

class ConflictBinPackingEngine::Impl {
public:
    Impl(const Instance& instance,
         std::uint64_t memo_entry_limit,
         int maximum_item_count)
        : instance_(instance),
          n_(instance.size()),
          capacity_(instance.capacity),
          blocks_((static_cast<std::size_t>(n_) + 63U) / 64U),
          maximum_item_count_(maximum_item_count),
          memo_(blocks_, memo_entry_limit) {
        if (n_ <= 0 || n_ > kMaximumConflictItems || capacity_ <= 0 ||
            maximum_item_count_ <= 0 ||
            maximum_item_count_ > kMaximumConflictItems ||
            memo_entry_limit >=
                std::numeric_limits<std::uint32_t>::max()) {
            throw std::invalid_argument(
                "invalid conflict-aware BINLB limits");
        }
        weights_.resize(static_cast<std::size_t>(n_));
        lb2_units_.resize(static_cast<std::size_t>(n_));
        lb3_units_.resize(static_cast<std::size_t>(n_));
        for (int item = 0; item < n_; ++item) {
            const int weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            if (weight <= 0 || weight > capacity_) {
                throw std::invalid_argument(
                    "invalid conflict-aware BINLB item weight");
            }
            weights_[static_cast<std::size_t>(item)] = weight;
            if (2LL * weight > capacity_) {
                lb2_units_[static_cast<std::size_t>(item)] = 6;
            } else if (2LL * weight == capacity_) {
                lb2_units_[static_cast<std::size_t>(item)] = 3;
            }
            if (3LL * weight > 2LL * capacity_) {
                lb3_units_[static_cast<std::size_t>(item)] = 6;
            } else if (3LL * weight == 2LL * capacity_) {
                lb3_units_[static_cast<std::size_t>(item)] = 4;
            } else if (3LL * weight > capacity_) {
                lb3_units_[static_cast<std::size_t>(item)] = 3;
            } else if (3LL * weight == capacity_) {
                lb3_units_[static_cast<std::size_t>(item)] = 2;
            }
        }

        build_conflict_graph();
        build_orders();

        const std::size_t depth_count = static_cast<std::size_t>(
            std::min(n_, maximum_item_count_) + 1);
        state_masks_.assign(depth_count * blocks_, 0U);
        load_masks_.assign(depth_count * blocks_, 0U);
        candidates_.assign(
            depth_count * static_cast<std::size_t>(maximum_item_count_), -1);
        greedy_bin_loads_.assign(
            static_cast<std::size_t>(maximum_item_count_), 0);
        greedy_bin_forbidden_.assign(
            static_cast<std::size_t>(maximum_item_count_) * blocks_, 0U);
        clique_mask_.assign(blocks_, 0U);
    }

    [[nodiscard]] ConflictBinPackingResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline,
        double maximum_call_seconds,
        std::uint64_t search_node_limit,
        std::uint64_t maximal_load_limit,
        int ordinary_lower_bound,
        int useful_lower_bound_target,
        ExactRelaxationLookup ordinary_lookup) {
        if (!std::isfinite(maximum_call_seconds) ||
            maximum_call_seconds <= 0.0 || search_node_limit == 0U ||
            maximal_load_limit == 0U || ordinary_lower_bound < 0) {
            throw std::invalid_argument(
                "invalid conflict-aware BINLB call limits");
        }

        ConflictBinPackingResult result;
        result.lower_bound = ordinary_lower_bound;
        if (!has_active_conflict(remaining_items)) {
            result.completed = true;
            result.optimum = ordinary_lower_bound;
            result.target_reached =
                ordinary_lower_bound >= useful_lower_bound_target;
            return result;
        }

        int item_count = 0;
        std::int64_t total_weight = 0;
        int lb2_units = 0;
        int lb3_units = 0;
        std::uint64_t* root = state_row(0);
        std::copy(remaining_items, remaining_items + blocks_, root);
        const unsigned tail = static_cast<unsigned>(n_) & 63U;
        if (tail != 0U) {
            root[blocks_ - 1U] &= (std::uint64_t{1} << tail) - 1U;
        }
        for (std::size_t block = 0; block < blocks_; ++block) {
            std::uint64_t value = root[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                const int item = static_cast<int>(block * 64U + bit);
                ++item_count;
                total_weight += weights_[static_cast<std::size_t>(item)];
                lb2_units += lb2_units_[static_cast<std::size_t>(item)];
                lb3_units += lb3_units_[static_cast<std::size_t>(item)];
                value &= value - 1U;
            }
        }
        if (item_count > maximum_item_count_) {
            return result;
        }

        result.attempted = true;
        const auto start = Clock::now();
        call_end_ = std::min(
            global_deadline.end_time(),
            start + std::chrono::duration_cast<Clock::duration>(
                        std::chrono::duration<double>(maximum_call_seconds)));
        call_node_limit_ = search_node_limit;
        call_load_limit_ = maximal_load_limit;
        ordinary_lookup_ = ordinary_lookup;
        call_nodes_ = 0U;
        call_loads_ = 0U;
        call_memo_hits_ = 0U;
        call_ordinary_memo_hits_ = 0U;
        aborted_ = false;
        timed_out_ = false;
        node_limited_ = false;
        load_limited_ = false;

        int root_lower = cheap_lower_bound(
            root, total_weight, lb2_units, lb3_units, false);
        root_lower = std::max(root_lower, ordinary_lower_bound);
        result.lower_bound = root_lower;
        if (root_lower >= useful_lower_bound_target) {
            result.target_test_completed = true;
            result.target_reached = true;
        } else if (useful_lower_bound_target <
                   std::numeric_limits<int>::max()) {
            const TargetDecision decision = can_pack_with_bins(
                0, total_weight, item_count, lb2_units, lb3_units,
                useful_lower_bound_target - 1, ordinary_lower_bound);
            if (decision != TargetDecision::kUnknown) {
                result.target_test_completed = true;
                result.target_reached =
                    decision == TargetDecision::kInfeasible;
                if (result.target_reached) {
                    result.lower_bound = useful_lower_bound_target;
                } else if (useful_lower_bound_target ==
                           ordinary_lower_bound + 1) {
                    result.completed = true;
                    result.optimum = ordinary_lower_bound;
                    memo_.insert(root, ordinary_lower_bound);
                }
            }
        } else {
            const NodeResult root_result = solve_state(
                0, total_weight, item_count, lb2_units, lb3_units,
                ordinary_lower_bound);
            result.completed = root_result.completed;
            result.optimum = root_result.completed
                ? root_result.optimum
                : 0;
            result.lower_bound = root_result.completed
                ? root_result.optimum
                : std::max(root_lower, root_result.lower_bound);
            result.target_reached =
                result.lower_bound >= useful_lower_bound_target;
        }

        result.timed_out = timed_out_;
        result.node_limited = node_limited_;
        result.load_limited = load_limited_;
        result.search_nodes = call_nodes_;
        result.maximal_loads = call_loads_;
        result.memo_hits = call_memo_hits_;
        result.ordinary_memo_hits = call_ordinary_memo_hits_;
        result.seconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        if (result.completed && Clock::now() > call_end_) {
            result.completed = false;
            result.optimum = 0;
            result.timed_out = true;
        }
        return result;
    }

    [[nodiscard]] bool has_active_conflict(
        const std::uint64_t* remaining_items) const noexcept {
        if (conflict_edge_count_ == 0U) {
            return false;
        }
        for (std::size_t block = 0; block < blocks_; ++block) {
            std::uint64_t value = remaining_items[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                const int item = static_cast<int>(block * 64U + bit);
                if (item >= n_) {
                    break;
                }
                const std::uint64_t* conflicts = conflict_row(item);
                for (std::size_t other_block = 0; other_block < blocks_;
                     ++other_block) {
                    if ((conflicts[other_block] &
                         remaining_items[other_block]) != 0U) {
                        return true;
                    }
                }
                value &= value - 1U;
            }
        }
        return false;
    }

    [[nodiscard]] int quick_lower_bound(
        const std::uint64_t* remaining_items) {
        std::uint64_t* state = state_row(0);
        std::copy(remaining_items, remaining_items + blocks_, state);
        const unsigned tail = static_cast<unsigned>(n_) & 63U;
        if (tail != 0U) {
            state[blocks_ - 1U] &= (std::uint64_t{1} << tail) - 1U;
        }
        std::int64_t total_weight = 0;
        int lb2_units = 0;
        int lb3_units = 0;
        for (std::size_t block = 0; block < blocks_; ++block) {
            std::uint64_t value = state[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                const int item = static_cast<int>(block * 64U + bit);
                total_weight += weights_[static_cast<std::size_t>(item)];
                lb2_units += lb2_units_[static_cast<std::size_t>(item)];
                lb3_units += lb3_units_[static_cast<std::size_t>(item)];
                value &= value - 1U;
            }
        }
        return cheap_lower_bound(
            state, total_weight, lb2_units, lb3_units, false);
    }

    [[nodiscard]] std::uint64_t conflict_edges() const noexcept {
        return conflict_edge_count_;
    }

    [[nodiscard]] std::uint64_t memo_entries() const noexcept {
        return memo_.size();
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        std::uint64_t bytes = memo_.memory_bytes();
        bytes = saturated_add(bytes, vector_memory_bytes(conflict_masks_));
        bytes = saturated_add(bytes, vector_memory_bytes(degrees_));
        bytes = saturated_add(bytes, vector_memory_bytes(order_));
        bytes = saturated_add(bytes, vector_memory_bytes(clique_order_));
        bytes = saturated_add(bytes, vector_memory_bytes(weights_));
        bytes = saturated_add(bytes, vector_memory_bytes(lb2_units_));
        bytes = saturated_add(bytes, vector_memory_bytes(lb3_units_));
        bytes = saturated_add(bytes, vector_memory_bytes(state_masks_));
        bytes = saturated_add(bytes, vector_memory_bytes(load_masks_));
        bytes = saturated_add(bytes, vector_memory_bytes(candidates_));
        bytes = saturated_add(bytes, vector_memory_bytes(greedy_bin_loads_));
        bytes = saturated_add(
            bytes, vector_memory_bytes(greedy_bin_forbidden_));
        return saturated_add(bytes, vector_memory_bytes(clique_mask_));
    }

private:
    [[nodiscard]] std::uint64_t* conflict_row(int item) noexcept {
        return conflict_masks_.data() +
               static_cast<std::size_t>(item) * blocks_;
    }

    [[nodiscard]] const std::uint64_t* conflict_row(
        int item) const noexcept {
        return conflict_masks_.data() +
               static_cast<std::size_t>(item) * blocks_;
    }

    [[nodiscard]] std::uint64_t* state_row(int depth) noexcept {
        return state_masks_.data() +
               static_cast<std::size_t>(depth) * blocks_;
    }

    [[nodiscard]] std::uint64_t* load_row(int depth) noexcept {
        return load_masks_.data() +
               static_cast<std::size_t>(depth) * blocks_;
    }

    [[nodiscard]] int* candidate_row(int depth) noexcept {
        return candidates_.data() +
               static_cast<std::size_t>(depth) * maximum_item_count_;
    }

    void add_conflict(int lhs, int rhs) noexcept {
        if (bit_is_set(conflict_row(lhs), rhs)) {
            return;
        }
        set_bit(conflict_row(lhs), rhs);
        set_bit(conflict_row(rhs), lhs);
        ++degrees_[static_cast<std::size_t>(lhs)];
        ++degrees_[static_cast<std::size_t>(rhs)];
        ++conflict_edge_count_;
    }

    void build_conflict_graph() {
        conflict_masks_.assign(static_cast<std::size_t>(n_) * blocks_, 0U);
        degrees_.assign(static_cast<std::size_t>(n_), 0);

        for (int lhs = 0; lhs < n_; ++lhs) {
            for (int rhs = lhs + 1; rhs < n_; ++rhs) {
                if (items_have_same_bin_conflict(instance_, lhs, rhs)) {
                    add_conflict(lhs, rhs);
                }
            }
        }
    }

    void build_orders() {
        order_.resize(static_cast<std::size_t>(n_));
        std::iota(order_.begin(), order_.end(), 0);
        std::sort(order_.begin(), order_.end(), [&](int lhs, int rhs) {
            const int lhs_weight = weights_[static_cast<std::size_t>(lhs)];
            const int rhs_weight = weights_[static_cast<std::size_t>(rhs)];
            if (lhs_weight != rhs_weight) {
                return lhs_weight > rhs_weight;
            }
            const int lhs_degree = degrees_[static_cast<std::size_t>(lhs)];
            const int rhs_degree = degrees_[static_cast<std::size_t>(rhs)];
            return lhs_degree != rhs_degree ? lhs_degree > rhs_degree
                                             : lhs < rhs;
        });
        clique_order_ = order_;
        std::sort(clique_order_.begin(), clique_order_.end(),
                  [&](int lhs, int rhs) {
            const int lhs_degree = degrees_[static_cast<std::size_t>(lhs)];
            const int rhs_degree = degrees_[static_cast<std::size_t>(rhs)];
            if (lhs_degree != rhs_degree) {
                return lhs_degree > rhs_degree;
            }
            const int lhs_weight = weights_[static_cast<std::size_t>(lhs)];
            const int rhs_weight = weights_[static_cast<std::size_t>(rhs)];
            return lhs_weight != rhs_weight ? lhs_weight > rhs_weight
                                             : lhs < rhs;
        });
    }

    [[nodiscard]] bool consume_search_node() noexcept {
        if (aborted_) {
            return false;
        }
        if (call_nodes_ >= call_node_limit_) {
            node_limited_ = true;
            aborted_ = true;
            return false;
        }
        ++call_nodes_;
        if ((call_nodes_ & 255U) == 1U && Clock::now() >= call_end_) {
            timed_out_ = true;
            aborted_ = true;
            return false;
        }
        return true;
    }

    [[nodiscard]] int clique_lower_bound(const std::uint64_t* state) {
        std::fill(clique_mask_.begin(), clique_mask_.end(), 0U);
        int count = 0;
        for (const int item : clique_order_) {
            if (!bit_is_set(state, item)) {
                continue;
            }
            const std::uint64_t* conflicts = conflict_row(item);
            bool adjacent = true;
            for (std::size_t block = 0; block < blocks_; ++block) {
                if ((clique_mask_[block] & ~conflicts[block]) != 0U) {
                    adjacent = false;
                    break;
                }
            }
            if (adjacent) {
                set_bit(clique_mask_.data(), item);
                ++count;
            }
        }
        return count;
    }

    [[nodiscard]] int cheap_lower_bound(const std::uint64_t* state,
                                        std::int64_t total_weight,
                                        int lb2_units,
                                        int lb3_units,
                                        bool query_ordinary) {
        if (total_weight == 0) {
            return 0;
        }
        int bound = std::max({ceil_div(total_weight, capacity_),
                              ceil_div(lb2_units, 6),
                              ceil_div(lb3_units, 6),
                              clique_lower_bound(state)});
        int ordinary_optimum = 0;
        if (query_ordinary && ordinary_lookup_.find != nullptr &&
            ordinary_lookup_.find(
                ordinary_lookup_.context, state, &ordinary_optimum)) {
            ++call_ordinary_memo_hits_;
            bound = std::max(bound, ordinary_optimum);
        }
        return bound;
    }

    [[nodiscard]] int greedy_upper_bound(const std::uint64_t* state) {
        std::fill(greedy_bin_loads_.begin(), greedy_bin_loads_.end(), 0);
        std::fill(greedy_bin_forbidden_.begin(),
                  greedy_bin_forbidden_.end(), 0U);
        int bin_count = 0;
        for (const int item : order_) {
            if (!bit_is_set(state, item)) {
                continue;
            }
            const int weight = weights_[static_cast<std::size_t>(item)];
            int best_bin = -1;
            int best_load = -1;
            for (int bin = 0; bin < bin_count; ++bin) {
                const int load =
                    greedy_bin_loads_[static_cast<std::size_t>(bin)];
                if (load + weight > capacity_ || load <= best_load) {
                    continue;
                }
                const std::uint64_t* forbidden =
                    greedy_bin_forbidden_.data() +
                    static_cast<std::size_t>(bin) * blocks_;
                if (!bit_is_set(forbidden, item)) {
                    best_bin = bin;
                    best_load = load;
                }
            }
            if (best_bin < 0) {
                best_bin = bin_count++;
            }
            greedy_bin_loads_[static_cast<std::size_t>(best_bin)] += weight;
            std::uint64_t* forbidden = greedy_bin_forbidden_.data() +
                static_cast<std::size_t>(best_bin) * blocks_;
            const std::uint64_t* conflicts = conflict_row(item);
            for (std::size_t block = 0; block < blocks_; ++block) {
                forbidden[block] |= conflicts[block];
            }
        }
        return bin_count;
    }

    [[nodiscard]] int minimum_improving_load(
        std::int64_t state_weight, int best) const noexcept {
        const std::int64_t threshold = state_weight -
            static_cast<std::int64_t>(best - 2) * capacity_;
        return static_cast<int>(std::max<std::int64_t>(0, threshold));
    }

    [[nodiscard]] int minimum_target_load(
        std::int64_t state_weight, int bin_limit) const noexcept {
        const std::int64_t threshold = state_weight -
            static_cast<std::int64_t>(bin_limit - 1) * capacity_;
        return static_cast<int>(std::max<std::int64_t>(0, threshold));
    }

    [[nodiscard]] bool load_is_maximal(int depth,
                                       int candidate_count,
                                       const ConflictMask& blocked,
                                       int residual_capacity) const noexcept {
        const int* candidates = candidates_.data() +
            static_cast<std::size_t>(depth) * maximum_item_count_;
        const std::uint64_t* load = load_masks_.data() +
            static_cast<std::size_t>(depth) * blocks_;
        for (int position = 0; position < candidate_count; ++position) {
            const int item = candidates[position];
            if (!bit_is_set(load, item) &&
                weights_[static_cast<std::size_t>(item)] <=
                    residual_capacity &&
                !bit_is_set(blocked.data(), item)) {
                return false;
            }
        }
        return true;
    }

    bool enumerate_target_loads(
        int depth,
        int position,
        int candidate_count,
        std::int64_t remaining_candidate_weight,
        std::int64_t state_weight,
        int state_item_count,
        int load_weight,
        int load_item_count,
        int load_lb2_units,
        int load_lb3_units,
        int state_lb2_units,
        int state_lb3_units,
        int bin_limit,
        std::uint64_t& state_load_count,
        const ConflictMask& blocked,
        bool& feasible) {
        if (!consume_search_node()) {
            return false;
        }
        if (static_cast<std::int64_t>(load_weight) +
                remaining_candidate_weight <
            minimum_target_load(state_weight, bin_limit)) {
            return true;
        }
        if (position == candidate_count) {
            if (load_weight < minimum_target_load(state_weight, bin_limit) ||
                !load_is_maximal(
                    depth, candidate_count, blocked,
                    capacity_ - load_weight)) {
                return true;
            }
            if (state_load_count >= call_load_limit_) {
                load_limited_ = true;
                aborted_ = true;
                return false;
            }
            ++state_load_count;
            ++call_loads_;

            const std::uint64_t* state = state_row(depth);
            const std::uint64_t* load = load_row(depth);
            std::uint64_t* child = state_row(depth + 1);
            for (std::size_t block = 0; block < blocks_; ++block) {
                child[block] = state[block] & ~load[block];
            }
            const std::int64_t child_weight = state_weight - load_weight;
            const int child_lb2 = state_lb2_units - load_lb2_units;
            const int child_lb3 = state_lb3_units - load_lb3_units;
            const int child_lower = cheap_lower_bound(
                child, child_weight, child_lb2, child_lb3, false);
            if (child_lower > bin_limit - 1) {
                memo_.record_infeasible(child, bin_limit - 1);
                return true;
            }
            const TargetDecision child_result = can_pack_with_bins(
                depth + 1, child_weight,
                state_item_count - load_item_count,
                child_lb2, child_lb3, bin_limit - 1, 0);
            if (child_result == TargetDecision::kUnknown) {
                return false;
            }
            if (child_result == TargetDecision::kFeasible) {
                feasible = true;
            }
            return true;
        }

        const int item = candidate_row(depth)[position];
        const int weight = weights_[static_cast<std::size_t>(item)];
        const std::int64_t remaining_after =
            remaining_candidate_weight - weight;
        std::uint64_t* load = load_row(depth);
        if (load_weight + weight <= capacity_ &&
            !bit_is_set(blocked.data(), item)) {
            set_bit(load, item);
            ConflictMask next_blocked;
            const std::uint64_t* conflicts = conflict_row(item);
            for (std::size_t block = 0; block < blocks_; ++block) {
                next_blocked[block] = blocked[block] | conflicts[block];
            }
            if (!enumerate_target_loads(
                    depth, position + 1, candidate_count,
                    remaining_after, state_weight, state_item_count,
                    load_weight + weight, load_item_count + 1,
                    load_lb2_units +
                        lb2_units_[static_cast<std::size_t>(item)],
                    load_lb3_units +
                        lb3_units_[static_cast<std::size_t>(item)],
                    state_lb2_units, state_lb3_units, bin_limit,
                    state_load_count, next_blocked, feasible)) {
                clear_bit(load, item);
                return false;
            }
            clear_bit(load, item);
            if (feasible) {
                return true;
            }
        }
        return enumerate_target_loads(
            depth, position + 1, candidate_count, remaining_after,
            state_weight, state_item_count, load_weight, load_item_count,
            load_lb2_units, load_lb3_units, state_lb2_units,
            state_lb3_units, bin_limit, state_load_count, blocked,
            feasible);
    }

    [[nodiscard]] TargetDecision can_pack_with_bins(
        int depth,
        std::int64_t total_weight,
        int item_count,
        int lb2_units,
        int lb3_units,
        int bin_limit,
        int external_lower_bound) {
        if (!consume_search_node()) {
            return TargetDecision::kUnknown;
        }
        if (item_count == 0) {
            return TargetDecision::kFeasible;
        }
        if (bin_limit <= 0) {
            return TargetDecision::kInfeasible;
        }
        if (item_count <= bin_limit) {
            return TargetDecision::kFeasible;
        }

        std::uint64_t* state = state_row(depth);
        const int cached_decision = memo_.query_target(state, bin_limit);
        if (cached_decision != 0) {
            ++call_memo_hits_;
            return cached_decision > 0
                ? TargetDecision::kFeasible
                : TargetDecision::kInfeasible;
        }
        int lower_bound = std::max(
            external_lower_bound,
            cheap_lower_bound(state, total_weight, lb2_units,
                              lb3_units, false));
        if (lower_bound <= bin_limit && external_lower_bound == 0 &&
            ordinary_lookup_.find != nullptr) {
            int ordinary_optimum = 0;
            if (ordinary_lookup_.find(
                    ordinary_lookup_.context, state, &ordinary_optimum)) {
                ++call_ordinary_memo_hits_;
                lower_bound = std::max(lower_bound, ordinary_optimum);
            }
        }
        if (lower_bound > bin_limit) {
            memo_.record_infeasible(state, bin_limit);
            return TargetDecision::kInfeasible;
        }
        if (greedy_upper_bound(state) <= bin_limit) {
            memo_.record_feasible(state, bin_limit);
            return TargetDecision::kFeasible;
        }

        int anchor = -1;
        for (const int item : order_) {
            if (bit_is_set(state, item)) {
                anchor = item;
                break;
            }
        }
        if (anchor < 0) {
            throw std::logic_error(
                "nonempty conflict-aware BINLB decision state has no item");
        }
        std::uint64_t* load = load_row(depth);
        std::fill(load, load + blocks_, 0U);
        set_bit(load, anchor);
        const int anchor_weight = weights_[static_cast<std::size_t>(anchor)];
        const std::uint64_t* anchor_conflicts = conflict_row(anchor);
        ConflictMask blocked;
        std::copy(anchor_conflicts, anchor_conflicts + blocks_,
                  blocked.begin());

        int candidate_count = 0;
        std::int64_t candidate_weight = 0;
        int* candidates = candidate_row(depth);
        for (const int item : order_) {
            if (item == anchor || !bit_is_set(state, item) ||
                bit_is_set(anchor_conflicts, item)) {
                continue;
            }
            const int weight = weights_[static_cast<std::size_t>(item)];
            if (anchor_weight + weight <= capacity_) {
                candidates[candidate_count++] = item;
                candidate_weight += weight;
            }
        }

        std::uint64_t state_load_count = 0U;
        bool feasible = false;
        if (!enumerate_target_loads(
                depth, 0, candidate_count, candidate_weight, total_weight,
                item_count, anchor_weight, 1,
                lb2_units_[static_cast<std::size_t>(anchor)],
                lb3_units_[static_cast<std::size_t>(anchor)],
                lb2_units, lb3_units, bin_limit,
                state_load_count, blocked, feasible) || aborted_) {
            return TargetDecision::kUnknown;
        }
        if (feasible) {
            memo_.record_feasible(state, bin_limit);
            return TargetDecision::kFeasible;
        }
        memo_.record_infeasible(state, bin_limit);
        return TargetDecision::kInfeasible;
    }

    bool enumerate_loads(int depth,
                         int position,
                         int candidate_count,
                         std::int64_t remaining_candidate_weight,
                         std::int64_t state_weight,
                         int state_item_count,
                         int load_weight,
                         int load_item_count,
                         int load_lb2_units,
                         int load_lb3_units,
                         int state_lb2_units,
                         int state_lb3_units,
                         int state_lower_bound,
                         int& best,
                         std::uint64_t& state_load_count,
                         const ConflictMask& blocked) {
        if (!consume_search_node()) {
            return false;
        }
        if (best == state_lower_bound) {
            return true;
        }
        if (static_cast<std::int64_t>(load_weight) +
                remaining_candidate_weight <
            minimum_improving_load(state_weight, best)) {
            return true;
        }
        if (position == candidate_count) {
            if (load_weight < minimum_improving_load(state_weight, best) ||
                !load_is_maximal(
                    depth, candidate_count, blocked,
                    capacity_ - load_weight)) {
                return true;
            }
            if (state_load_count >= call_load_limit_) {
                load_limited_ = true;
                aborted_ = true;
                return false;
            }
            ++state_load_count;
            ++call_loads_;

            const std::uint64_t* state = state_row(depth);
            const std::uint64_t* load = load_row(depth);
            std::uint64_t* child = state_row(depth + 1);
            for (std::size_t block = 0; block < blocks_; ++block) {
                child[block] = state[block] & ~load[block];
            }
            const std::int64_t child_weight = state_weight - load_weight;
            const int child_lb2 = state_lb2_units - load_lb2_units;
            const int child_lb3 = state_lb3_units - load_lb3_units;
            const int child_lower = cheap_lower_bound(
                child, child_weight, child_lb2, child_lb3, false);
            if (1 + child_lower >= best) {
                return true;
            }
            const NodeResult child_result = solve_state(
                depth + 1, child_weight,
                state_item_count - load_item_count,
                child_lb2, child_lb3, 0);
            if (!child_result.completed) {
                return false;
            }
            best = std::min(best, 1 + child_result.optimum);
            return true;
        }

        const int item = candidate_row(depth)[position];
        const int weight = weights_[static_cast<std::size_t>(item)];
        const std::int64_t remaining_after =
            remaining_candidate_weight - weight;
        std::uint64_t* load = load_row(depth);
        if (load_weight + weight <= capacity_ &&
            !bit_is_set(blocked.data(), item)) {
            set_bit(load, item);
            ConflictMask next_blocked;
            const std::uint64_t* conflicts = conflict_row(item);
            for (std::size_t block = 0; block < blocks_; ++block) {
                next_blocked[block] = blocked[block] | conflicts[block];
            }
            if (!enumerate_loads(
                    depth, position + 1, candidate_count,
                    remaining_after, state_weight, state_item_count,
                    load_weight + weight, load_item_count + 1,
                    load_lb2_units +
                        lb2_units_[static_cast<std::size_t>(item)],
                    load_lb3_units +
                        lb3_units_[static_cast<std::size_t>(item)],
                    state_lb2_units, state_lb3_units, state_lower_bound,
                    best, state_load_count, next_blocked)) {
                clear_bit(load, item);
                return false;
            }
            clear_bit(load, item);
            if (best == state_lower_bound) {
                return true;
            }
        }
        return enumerate_loads(
            depth, position + 1, candidate_count, remaining_after,
            state_weight, state_item_count, load_weight, load_item_count,
            load_lb2_units, load_lb3_units, state_lb2_units,
            state_lb3_units, state_lower_bound, best, state_load_count,
            blocked);
    }

    [[nodiscard]] NodeResult solve_state(int depth,
                                         std::int64_t total_weight,
                                         int item_count,
                                         int lb2_units,
                                         int lb3_units,
                                         int external_lower_bound) {
        NodeResult result;
        if (!consume_search_node()) {
            return result;
        }
        if (item_count == 0) {
            result.completed = true;
            return result;
        }
        std::uint64_t* state = state_row(depth);
        int memo_value = 0;
        if (memo_.find(state, &memo_value)) {
            ++call_memo_hits_;
            result.completed = true;
            result.lower_bound = memo_value;
            result.optimum = memo_value;
            return result;
        }

        result.lower_bound = std::max(
            external_lower_bound,
            cheap_lower_bound(state, total_weight, lb2_units,
                              lb3_units, false));
        int best = greedy_upper_bound(state);
        if (external_lower_bound == 0 && result.lower_bound < best &&
            ordinary_lookup_.find != nullptr) {
            int ordinary_optimum = 0;
            if (ordinary_lookup_.find(
                    ordinary_lookup_.context, state, &ordinary_optimum)) {
                ++call_ordinary_memo_hits_;
                result.lower_bound = std::max(
                    result.lower_bound, ordinary_optimum);
            }
        }
        if (best == result.lower_bound) {
            result.completed = true;
            result.optimum = best;
            memo_.insert(state, best);
            return result;
        }

        int anchor = -1;
        for (const int item : order_) {
            if (bit_is_set(state, item)) {
                anchor = item;
                break;
            }
        }
        if (anchor < 0) {
            throw std::logic_error(
                "nonempty conflict-aware BINLB state has no item");
        }
        std::uint64_t* load = load_row(depth);
        std::fill(load, load + blocks_, 0U);
        set_bit(load, anchor);
        const int anchor_weight = weights_[static_cast<std::size_t>(anchor)];
        const std::uint64_t* anchor_conflicts = conflict_row(anchor);
        ConflictMask blocked;
        std::copy(anchor_conflicts, anchor_conflicts + blocks_,
                  blocked.begin());

        int candidate_count = 0;
        std::int64_t candidate_weight = 0;
        int* candidates = candidate_row(depth);
        for (const int item : order_) {
            if (item == anchor || !bit_is_set(state, item) ||
                bit_is_set(anchor_conflicts, item)) {
                continue;
            }
            const int weight = weights_[static_cast<std::size_t>(item)];
            if (anchor_weight + weight <= capacity_) {
                candidates[candidate_count++] = item;
                candidate_weight += weight;
            }
        }

        std::uint64_t state_load_count = 0U;
        const bool enumerated = enumerate_loads(
            depth, 0, candidate_count, candidate_weight, total_weight,
            item_count, anchor_weight, 1,
            lb2_units_[static_cast<std::size_t>(anchor)],
            lb3_units_[static_cast<std::size_t>(anchor)],
            lb2_units, lb3_units, result.lower_bound, best,
            state_load_count, blocked);
        if (!enumerated || aborted_) {
            return result;
        }
        result.completed = true;
        result.optimum = best;
        memo_.insert(state, best);
        return result;
    }

    const Instance& instance_;
    int n_ = 0;
    int capacity_ = 0;
    std::size_t blocks_ = 0;
    int maximum_item_count_ = 0;
    std::uint64_t conflict_edge_count_ = 0U;
    std::vector<std::uint64_t> conflict_masks_;
    std::vector<int> degrees_;
    std::vector<int> order_;
    std::vector<int> clique_order_;
    std::vector<int> weights_;
    std::vector<int> lb2_units_;
    std::vector<int> lb3_units_;
    ExactSubsetMemo memo_;

    std::vector<std::uint64_t> state_masks_;
    std::vector<std::uint64_t> load_masks_;
    std::vector<int> candidates_;
    std::vector<int> greedy_bin_loads_;
    std::vector<std::uint64_t> greedy_bin_forbidden_;
    std::vector<std::uint64_t> clique_mask_;

    Clock::time_point call_end_{};
    std::uint64_t call_node_limit_ = 0U;
    std::uint64_t call_load_limit_ = 0U;
    ExactRelaxationLookup ordinary_lookup_;
    std::uint64_t call_nodes_ = 0U;
    std::uint64_t call_loads_ = 0U;
    std::uint64_t call_memo_hits_ = 0U;
    std::uint64_t call_ordinary_memo_hits_ = 0U;
    bool aborted_ = false;
    bool timed_out_ = false;
    bool node_limited_ = false;
    bool load_limited_ = false;
};

ConflictBinPackingEngine::ConflictBinPackingEngine(
    const Instance& instance,
    std::uint64_t memo_entry_limit,
    int maximum_item_count)
    : impl_(std::make_unique<Impl>(
          instance, memo_entry_limit, maximum_item_count)) {}

ConflictBinPackingEngine::~ConflictBinPackingEngine() = default;
ConflictBinPackingEngine::ConflictBinPackingEngine(
    ConflictBinPackingEngine&&) noexcept = default;
ConflictBinPackingEngine& ConflictBinPackingEngine::operator=(
    ConflictBinPackingEngine&&) noexcept = default;

ConflictBinPackingResult ConflictBinPackingEngine::solve(
    const std::uint64_t* remaining_items,
    Deadline& global_deadline,
    double maximum_call_seconds,
    std::uint64_t search_node_limit,
    std::uint64_t maximal_load_limit,
    int ordinary_lower_bound,
    int useful_lower_bound_target,
    ExactRelaxationLookup ordinary_lookup) {
    if (remaining_items == nullptr) {
        throw std::invalid_argument(
            "null conflict-aware BINLB item set");
    }
    return impl_->solve(
        remaining_items, global_deadline, maximum_call_seconds,
        search_node_limit, maximal_load_limit, ordinary_lower_bound,
        useful_lower_bound_target, ordinary_lookup);
}

bool ConflictBinPackingEngine::has_active_conflict(
    const std::uint64_t* remaining_items) const noexcept {
    return impl_->has_active_conflict(remaining_items);
}

int ConflictBinPackingEngine::quick_lower_bound(
    const std::uint64_t* remaining_items) {
    if (remaining_items == nullptr) {
        throw std::invalid_argument(
            "null conflict-aware BINLB item set");
    }
    return impl_->quick_lower_bound(remaining_items);
}

std::uint64_t ConflictBinPackingEngine::conflict_edge_count() const noexcept {
    return impl_->conflict_edges();
}

std::uint64_t ConflictBinPackingEngine::memo_entry_count() const noexcept {
    return impl_->memo_entries();
}

std::uint64_t ConflictBinPackingEngine::memory_bytes() const noexcept {
    return impl_->memory_bytes();
}

}
