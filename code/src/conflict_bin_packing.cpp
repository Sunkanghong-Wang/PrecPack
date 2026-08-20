#include "precpack/conflict_bin_packing.hpp"

#include "precpack/bin_packing_bound.hpp"
#include "precpack/conflict_graph.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace precpack {
namespace {

using Clock = std::chrono::steady_clock;

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
    for (std::size_t index = 0; index < count; ++index) {
        hash = mix64(hash ^ mix64(words[index] +
                                  0x9e3779b97f4a7c15ULL * (index + 1U)));
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

[[nodiscard]] int ceil_div(std::int64_t numerator, int denominator) noexcept {
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
        if (maximum_entries_ > 0U) {
            slots_.assign(1024U, 0U);
        }
    }

    [[nodiscard]] bool find(const std::uint64_t* key,
                            int* value) const noexcept {
        if (slots_.empty()) {
            return false;
        }
        const std::uint64_t hash = hash_words(key, blocks_);
        std::size_t slot = static_cast<std::size_t>(hash) & (slots_.size() - 1U);
        while (slots_[slot] != 0U) {
            const std::size_t entry = slots_[slot] - 1U;
            if (hashes_[entry] == hash && keys_equal(entry, key)) {
                *value = values_[entry];
                return true;
            }
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
        return false;
    }

    void insert(const std::uint64_t* key, int value) {
        if (maximum_entries_ == 0U || values_.size() >= maximum_entries_) {
            return;
        }
        int existing = 0;
        if (find(key, &existing)) {
            return;
        }
        if ((values_.size() + 1U) * 10U >= slots_.size() * 7U) {
            rehash(slots_.size() * 2U);
        }
        const std::uint64_t hash = hash_words(key, blocks_);
        const std::size_t entry = values_.size();
        keys_.insert(keys_.end(), key, key + blocks_);
        hashes_.push_back(hash);
        values_.push_back(value);
        std::size_t slot = static_cast<std::size_t>(hash) & (slots_.size() - 1U);
        while (slots_[slot] != 0U) {
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
        slots_[slot] = static_cast<std::uint32_t>(entry + 1U);
    }

    [[nodiscard]] std::uint64_t size() const noexcept {
        return values_.size();
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        std::uint64_t bytes = vector_memory_bytes(slots_);
        bytes = saturated_add(bytes, vector_memory_bytes(keys_));
        bytes = saturated_add(bytes, vector_memory_bytes(hashes_));
        return saturated_add(bytes, vector_memory_bytes(values_));
    }

private:
    [[nodiscard]] bool keys_equal(std::size_t entry,
                                  const std::uint64_t* key) const noexcept {
        const std::uint64_t* stored = keys_.data() + entry * blocks_;
        for (std::size_t block = 0; block < blocks_; ++block) {
            if (stored[block] != key[block]) {
                return false;
            }
        }
        return true;
    }

    void rehash(std::size_t requested) {
        std::size_t capacity = 1024U;
        while (capacity < requested) {
            capacity *= 2U;
        }
        std::vector<std::uint32_t> replacement(capacity, 0U);
        for (std::size_t entry = 0; entry < values_.size(); ++entry) {
            std::size_t slot = static_cast<std::size_t>(hashes_[entry]) &
                               (capacity - 1U);
            while (replacement[slot] != 0U) {
                slot = (slot + 1U) & (capacity - 1U);
            }
            replacement[slot] = static_cast<std::uint32_t>(entry + 1U);
        }
        slots_.swap(replacement);
    }

    std::size_t blocks_ = 0;
    std::uint64_t maximum_entries_ = 0;
    std::vector<std::uint32_t> slots_;
    std::vector<std::uint64_t> keys_;
    std::vector<std::uint64_t> hashes_;
    std::vector<int> values_;
};

struct NodeResult {
    bool completed = false;
    int lower_bound = 0;
    int optimum = 0;
};

}

class ConflictBinPackingBound::Impl {
public:
    Impl(const Instance& instance,
         ConflictBinPackingLimits limits,
         BinPackingBound* ordinary_bound)
        : instance_(instance),
          limits_(limits),
          ordinary_bound_(ordinary_bound),
          n_(instance.size()),
          capacity_(instance.capacity),
          blocks_((static_cast<std::size_t>(n_) + 63U) / 64U),
          memo_(blocks_, limits.memo_entry_limit) {
        if (n_ <= 0 || capacity_ <= 0 || limits_.call_time_limit_seconds <= 0.0 ||
            limits_.total_time_limit_seconds <= 0.0 ||
            limits_.search_node_limit == 0U ||
            limits_.maximal_load_limit_per_state == 0U ||
            limits_.memo_entry_limit >=
                std::numeric_limits<std::uint32_t>::max() ||
            limits_.maximum_item_count <= 0) {
            throw std::invalid_argument("invalid conflict BINLB limits");
        }
        for (const Item& item : instance_.items) {
            if (item.weight <= 0 || item.weight > capacity_) {
                throw std::invalid_argument("invalid conflict BINLB item weight");
            }
        }

        const std::vector<std::pair<int, int>> edges =
            build_bppc_relaxation_conflict_edges(instance_);
        conflict_edge_count_ = edges.size();
        conflict_masks_.assign(static_cast<std::size_t>(n_) * blocks_, 0U);
        degrees_.assign(static_cast<std::size_t>(n_), 0);
        for (const auto [lhs, rhs] : edges) {
            set_bit(conflict_row(lhs), rhs);
            set_bit(conflict_row(rhs), lhs);
            ++degrees_[static_cast<std::size_t>(lhs)];
            ++degrees_[static_cast<std::size_t>(rhs)];
        }
        order_.resize(static_cast<std::size_t>(n_));
        std::iota(order_.begin(), order_.end(), 0);
        std::sort(order_.begin(), order_.end(), [&](int lhs, int rhs) {
            const int lhs_weight =
                instance_.items[static_cast<std::size_t>(lhs)].weight;
            const int rhs_weight =
                instance_.items[static_cast<std::size_t>(rhs)].weight;
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
            const int lhs_weight =
                instance_.items[static_cast<std::size_t>(lhs)].weight;
            const int rhs_weight =
                instance_.items[static_cast<std::size_t>(rhs)].weight;
            return lhs_weight != rhs_weight ? lhs_weight > rhs_weight
                                             : lhs < rhs;
        });

        const std::size_t depth_count = static_cast<std::size_t>(
            std::min(n_, limits_.maximum_item_count) + 1);
        state_masks_.assign(depth_count * blocks_, 0U);
        load_masks_.assign(depth_count * blocks_, 0U);
        candidates_.assign(
            depth_count * static_cast<std::size_t>(limits_.maximum_item_count),
            -1);
        greedy_bin_loads_.assign(static_cast<std::size_t>(n_), 0);
        greedy_bin_masks_.assign(static_cast<std::size_t>(n_) * blocks_, 0U);
        clique_mask_.assign(blocks_, 0U);
    }

    [[nodiscard]] ConflictBinPackingResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline) {
        ConflictBinPackingResult result;
        const double remaining_total_budget =
            limits_.total_time_limit_seconds - cumulative_seconds_;
        if (remaining_total_budget <= 0.0) {
            result.total_budget_exhausted = true;
            return result;
        }
        int item_count = 0;
        std::int64_t total_weight = 0;
        for (int item = 0; item < n_; ++item) {
            if (bit_is_set(remaining_items, item)) {
                ++item_count;
                total_weight +=
                    instance_.items[static_cast<std::size_t>(item)].weight;
            }
        }
        if (item_count > limits_.maximum_item_count) {
            result.item_limited = true;
            return result;
        }
        result.attempted = true;
        const auto start = Clock::now();
        const double local_seconds = std::min(
            limits_.call_time_limit_seconds, remaining_total_budget);
        call_end_ = std::min(
            global_deadline.end_time(),
            start + std::chrono::duration_cast<Clock::duration>(
                        std::chrono::duration<double>(local_seconds)));
        call_nodes_ = 0U;
        call_loads_ = 0U;
        call_memo_hits_ = 0U;
        call_ordinary_memo_hits_ = 0U;
        aborted_ = false;
        timed_out_ = false;
        node_limited_ = false;
        load_limited_ = false;

        std::copy(remaining_items, remaining_items + blocks_,
                  state_masks_.begin());
        const unsigned tail = static_cast<unsigned>(n_) & 63U;
        if (tail != 0U) {
            state_masks_[blocks_ - 1U] &=
                (std::uint64_t{1} << tail) - 1U;
        }
        const NodeResult root = solve_state(0, total_weight, item_count);
        result.completed = root.completed;
        result.lower_bound = root.lower_bound;
        result.optimum = root.completed ? root.optimum : 0;
        result.timed_out = timed_out_;
        result.node_limited = node_limited_;
        result.load_limited = load_limited_;
        result.search_nodes = call_nodes_;
        result.maximal_loads = call_loads_;
        result.memo_hits = call_memo_hits_;
        result.ordinary_memo_hits = call_ordinary_memo_hits_;
        result.seconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        cumulative_seconds_ += result.seconds;
        return result;
    }

    [[nodiscard]] std::size_t blocks() const noexcept { return blocks_; }
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
        bytes = saturated_add(bytes, vector_memory_bytes(state_masks_));
        bytes = saturated_add(bytes, vector_memory_bytes(load_masks_));
        bytes = saturated_add(bytes, vector_memory_bytes(candidates_));
        bytes = saturated_add(bytes, vector_memory_bytes(greedy_bin_loads_));
        bytes = saturated_add(bytes, vector_memory_bytes(greedy_bin_masks_));
        return saturated_add(bytes, vector_memory_bytes(clique_mask_));
    }

private:
    [[nodiscard]] std::uint64_t* conflict_row(int item) noexcept {
        return conflict_masks_.data() + static_cast<std::size_t>(item) * blocks_;
    }
    [[nodiscard]] const std::uint64_t* conflict_row(int item) const noexcept {
        return conflict_masks_.data() + static_cast<std::size_t>(item) * blocks_;
    }
    [[nodiscard]] std::uint64_t* state_row(int depth) noexcept {
        return state_masks_.data() + static_cast<std::size_t>(depth) * blocks_;
    }
    [[nodiscard]] std::uint64_t* load_row(int depth) noexcept {
        return load_masks_.data() + static_cast<std::size_t>(depth) * blocks_;
    }
    [[nodiscard]] int* candidate_row(int depth) noexcept {
        return candidates_.data() + static_cast<std::size_t>(depth) *
                                        limits_.maximum_item_count;
    }

    [[nodiscard]] bool consume_search_node() noexcept {
        if (aborted_) {
            return false;
        }
        if (call_nodes_ >= limits_.search_node_limit) {
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

    [[nodiscard]] bool conflicts_with_load(
        int item, const std::uint64_t* load) const noexcept {
        const std::uint64_t* conflicts = conflict_row(item);
        for (std::size_t block = 0; block < blocks_; ++block) {
            if ((conflicts[block] & load[block]) != 0U) {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] int clique_lower_bound(const std::uint64_t* state) {
        std::fill(clique_mask_.begin(), clique_mask_.end(), 0U);
        int count = 0;
        for (const int item : clique_order_) {
            if (!bit_is_set(state, item)) {
                continue;
            }
            const std::uint64_t* conflicts = conflict_row(item);
            bool adjacent_to_clique = true;
            for (std::size_t block = 0; block < blocks_; ++block) {
                if ((clique_mask_[block] & ~conflicts[block]) != 0U) {
                    adjacent_to_clique = false;
                    break;
                }
            }
            if (adjacent_to_clique) {
                set_bit(clique_mask_.data(), item);
                ++count;
            }
        }
        return count;
    }

    [[nodiscard]] int cheap_lower_bound(const std::uint64_t* state,
                                        std::int64_t total_weight,
                                        bool query_ordinary = true) {
        if (total_weight == 0) {
            return 0;
        }
        int bound = std::max(ceil_div(total_weight, capacity_),
                             clique_lower_bound(state));
        int ordinary_optimum = 0;
        if (query_ordinary && ordinary_bound_ != nullptr &&
            ordinary_bound_->lookup_exact(state, &ordinary_optimum)) {
            ++call_ordinary_memo_hits_;
            bound = std::max(bound, ordinary_optimum);
        }
        return bound;
    }

    [[nodiscard]] int greedy_upper_bound(const std::uint64_t* state) {
        std::fill(greedy_bin_loads_.begin(), greedy_bin_loads_.end(), 0);
        std::fill(greedy_bin_masks_.begin(), greedy_bin_masks_.end(), 0U);
        int bin_count = 0;
        for (const int item : order_) {
            if (!bit_is_set(state, item)) {
                continue;
            }
            const int weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            int best_bin = -1;
            int best_load = -1;
            for (int bin = 0; bin < bin_count; ++bin) {
                if (greedy_bin_loads_[static_cast<std::size_t>(bin)] + weight >
                    capacity_) {
                    continue;
                }
                const std::uint64_t* bin_mask =
                    greedy_bin_masks_.data() +
                    static_cast<std::size_t>(bin) * blocks_;
                if (!conflicts_with_load(item, bin_mask) &&
                    greedy_bin_loads_[static_cast<std::size_t>(bin)] >
                        best_load) {
                    best_bin = bin;
                    best_load =
                        greedy_bin_loads_[static_cast<std::size_t>(bin)];
                }
            }
            if (best_bin < 0) {
                best_bin = bin_count++;
            }
            greedy_bin_loads_[static_cast<std::size_t>(best_bin)] += weight;
            set_bit(greedy_bin_masks_.data() +
                        static_cast<std::size_t>(best_bin) * blocks_,
                    item);
        }
        return bin_count;
    }

    [[nodiscard]] bool load_is_maximal(const std::uint64_t* state,
                                       const std::uint64_t* load,
                                       int load_weight) const noexcept {
        const int residual_capacity = capacity_ - load_weight;
        for (const int item : order_) {
            if (!bit_is_set(state, item) || bit_is_set(load, item) ||
                instance_.items[static_cast<std::size_t>(item)].weight >
                    residual_capacity) {
                continue;
            }
            if (!conflicts_with_load(item, load)) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] int minimum_improving_load(
        std::int64_t state_weight, int best) const noexcept {
        const std::int64_t threshold =
            state_weight - static_cast<std::int64_t>(best - 2) * capacity_;
        return static_cast<int>(std::max<std::int64_t>(0, threshold));
    }

    bool enumerate_loads(int depth,
                         int position,
                         int candidate_count,
                         std::int64_t remaining_candidate_weight,
                         std::int64_t state_weight,
                         int state_item_count,
                         int load_weight,
                         int load_item_count,
                         int state_lower_bound,
                         int& best,
                         std::uint64_t& state_load_count) {
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
            std::uint64_t* state = state_row(depth);
            std::uint64_t* load = load_row(depth);
            if (load_weight < minimum_improving_load(state_weight, best) ||
                !load_is_maximal(state, load, load_weight)) {
                return true;
            }
            if (state_load_count >=
                limits_.maximal_load_limit_per_state) {
                load_limited_ = true;
                aborted_ = true;
                return false;
            }
            ++state_load_count;
            ++call_loads_;

            std::uint64_t* child = state_row(depth + 1);
            for (std::size_t block = 0; block < blocks_; ++block) {
                child[block] = state[block] & ~load[block];
            }
            const std::int64_t child_weight = state_weight - load_weight;
            const int child_lower =
                cheap_lower_bound(child, child_weight, false);
            if (1 + child_lower >= best) {
                return true;
            }
            const NodeResult child_result = solve_state(
                depth + 1, child_weight,
                state_item_count - load_item_count);
            if (!child_result.completed) {
                return false;
            }
            best = std::min(best, 1 + child_result.optimum);
            return true;
        }

        const int item = candidate_row(depth)[position];
        const int weight =
            instance_.items[static_cast<std::size_t>(item)].weight;
        std::uint64_t* load = load_row(depth);
        if (load_weight + weight <= capacity_ &&
            !conflicts_with_load(item, load)) {
            set_bit(load, item);
            if (!enumerate_loads(
                    depth, position + 1, candidate_count,
                    remaining_candidate_weight - weight, state_weight,
                    state_item_count, load_weight + weight,
                    load_item_count + 1, state_lower_bound, best,
                    state_load_count)) {
                clear_bit(load, item);
                return false;
            }
            clear_bit(load, item);
            if (best == state_lower_bound) {
                return true;
            }
        }
        return enumerate_loads(
            depth, position + 1, candidate_count,
            remaining_candidate_weight - weight, state_weight,
            state_item_count, load_weight, load_item_count,
            state_lower_bound, best, state_load_count);
    }

    [[nodiscard]] NodeResult solve_state(int depth,
                                         std::int64_t total_weight,
                                         int item_count) {
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

        result.lower_bound = cheap_lower_bound(state, total_weight);
        int best = greedy_upper_bound(state);
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
            throw std::logic_error("nonempty conflict BINLB state has no item");
        }
        std::uint64_t* load = load_row(depth);
        std::fill(load, load + blocks_, 0U);
        set_bit(load, anchor);
        const int anchor_weight =
            instance_.items[static_cast<std::size_t>(anchor)].weight;
        int candidate_count = 0;
        std::int64_t candidate_weight = 0;
        int* candidates = candidate_row(depth);
        for (const int item : order_) {
            if (item == anchor || !bit_is_set(state, item) ||
                bit_is_set(conflict_row(anchor), item)) {
                continue;
            }
            const int weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            if (anchor_weight + weight <= capacity_) {
                candidates[candidate_count++] = item;
                candidate_weight += weight;
            }
        }

        std::uint64_t state_load_count = 0U;
        const bool enumerated = enumerate_loads(
            depth, 0, candidate_count, candidate_weight, total_weight,
            item_count, anchor_weight, 1, result.lower_bound, best,
            state_load_count);
        if (!enumerated || aborted_) {
            return result;
        }
        result.completed = true;
        result.optimum = best;
        memo_.insert(state, best);
        return result;
    }

    const Instance& instance_;
    ConflictBinPackingLimits limits_;
    BinPackingBound* ordinary_bound_ = nullptr;
    int n_ = 0;
    int capacity_ = 0;
    std::size_t blocks_ = 0;
    std::uint64_t conflict_edge_count_ = 0;
    std::vector<std::uint64_t> conflict_masks_;
    std::vector<int> degrees_;
    std::vector<int> order_;
    std::vector<int> clique_order_;
    ExactSubsetMemo memo_;

    std::vector<std::uint64_t> state_masks_;
    std::vector<std::uint64_t> load_masks_;
    std::vector<int> candidates_;
    std::vector<int> greedy_bin_loads_;
    std::vector<std::uint64_t> greedy_bin_masks_;
    std::vector<std::uint64_t> clique_mask_;

    Clock::time_point call_end_{};
    double cumulative_seconds_ = 0.0;
    std::uint64_t call_nodes_ = 0U;
    std::uint64_t call_loads_ = 0U;
    std::uint64_t call_memo_hits_ = 0U;
    std::uint64_t call_ordinary_memo_hits_ = 0U;
    bool aborted_ = false;
    bool timed_out_ = false;
    bool node_limited_ = false;
    bool load_limited_ = false;
};

ConflictBinPackingBound::ConflictBinPackingBound(
    const Instance& instance,
    ConflictBinPackingLimits limits,
    BinPackingBound* ordinary_bound)
    : impl_(std::make_unique<Impl>(instance, limits, ordinary_bound)) {}

ConflictBinPackingBound::~ConflictBinPackingBound() = default;
ConflictBinPackingBound::ConflictBinPackingBound(
    ConflictBinPackingBound&&) noexcept = default;
ConflictBinPackingBound& ConflictBinPackingBound::operator=(
    ConflictBinPackingBound&&) noexcept = default;

ConflictBinPackingResult ConflictBinPackingBound::solve(
    const std::uint64_t* remaining_items,
    Deadline& global_deadline) {
    if (remaining_items == nullptr) {
        throw std::invalid_argument("null conflict BINLB item set");
    }
    return impl_->solve(remaining_items, global_deadline);
}

std::size_t ConflictBinPackingBound::bit_block_count() const noexcept {
    return impl_->blocks();
}

std::uint64_t ConflictBinPackingBound::conflict_edge_count() const noexcept {
    return impl_->conflict_edges();
}

std::uint64_t ConflictBinPackingBound::memo_entry_count() const noexcept {
    return impl_->memo_entries();
}

std::uint64_t ConflictBinPackingBound::memory_bytes() const noexcept {
    return impl_->memory_bytes();
}

}
