#include "precpack/bin_packing_bound.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
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

[[nodiscard]] std::uint64_t hash_counts(const std::uint16_t* counts,
                                        std::size_t count) noexcept {
    std::uint64_t hash = 0x510e527fade682d1ULL ^
                         mix64(static_cast<std::uint64_t>(count));
    for (std::size_t index = 0; index < count; ++index) {
        hash = mix64(hash ^ mix64(
            static_cast<std::uint64_t>(counts[index]) +
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

[[nodiscard]] int ceil_div(std::int64_t numerator,
                           std::int64_t denominator) noexcept {
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

class ExactMultiplicityMemo {
public:
    ExactMultiplicityMemo(std::size_t class_count,
                          std::uint64_t maximum_entries)
        : class_count_(class_count), maximum_entries_(maximum_entries) {
        if (maximum_entries_ > 0U) {
            slots_.assign(1024U, 0U);
        }
    }

    [[nodiscard]] bool find(const std::uint16_t* key,
                            int* value) const noexcept {
        if (slots_.empty()) {
            return false;
        }
        const std::uint64_t hash = hash_counts(key, class_count_);
        std::size_t slot = static_cast<std::size_t>(hash) &
                           (slots_.size() - 1U);
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

    void insert(const std::uint16_t* key, int value) {
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
        const std::uint64_t hash = hash_counts(key, class_count_);
        const std::size_t entry = values_.size();
        keys_.insert(keys_.end(), key, key + class_count_);
        hashes_.push_back(hash);
        values_.push_back(value);
        std::size_t slot = static_cast<std::size_t>(hash) &
                           (slots_.size() - 1U);
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
                                  const std::uint16_t* key) const noexcept {
        const std::uint16_t* stored =
            keys_.data() + entry * class_count_;
        for (std::size_t index = 0; index < class_count_; ++index) {
            if (stored[index] != key[index]) {
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

    std::size_t class_count_ = 0;
    std::uint64_t maximum_entries_ = 0;
    std::vector<std::uint32_t> slots_;
    std::vector<std::uint16_t> keys_;
    std::vector<std::uint64_t> hashes_;
    std::vector<int> values_;
};

struct NodeResult {
    bool completed = false;
    int lower_bound = 0;
    int optimum = 0;
};

}

class BinPackingBound::Impl {
public:
    Impl(const Instance& instance, BinPackingBoundLimits limits)
        : instance_(instance),
          limits_(limits),
          n_(instance.size()),
          capacity_(instance.capacity),
          blocks_((static_cast<std::size_t>(n_) + 63U) / 64U) {
        if (n_ <= 0 || capacity_ <= 0 ||
            limits_.call_time_limit_seconds <= 0.0 ||
            limits_.search_node_limit == 0U ||
            limits_.nondominated_load_limit_per_state == 0U ||
            limits_.memo_entry_limit >=
                std::numeric_limits<std::uint32_t>::max() ||
            limits_.maximum_item_count <= 0 ||
            limits_.maximum_item_count >
                std::numeric_limits<std::uint16_t>::max()) {
            throw std::invalid_argument("invalid ordinary BINLB limits");
        }
        weights_.reserve(static_cast<std::size_t>(n_));
        for (const Item& item : instance_.items) {
            if (item.weight <= 0 || item.weight > capacity_) {
                throw std::invalid_argument("invalid ordinary BINLB item weight");
            }
            weights_.push_back(item.weight);
        }
        std::sort(weights_.begin(), weights_.end(), std::greater<>());
        weights_.erase(std::unique(weights_.begin(), weights_.end()),
                       weights_.end());
        class_count_ = weights_.size();
        if (class_count_ >
            static_cast<std::size_t>(
                std::numeric_limits<std::uint16_t>::max())) {
            throw std::invalid_argument(
                "too many ordinary BINLB weight classes");
        }
        item_class_.assign(static_cast<std::size_t>(n_), 0U);
        for (int item = 0; item < n_; ++item) {
            const int weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            const auto position = std::lower_bound(
                weights_.begin(), weights_.end(), weight, std::greater<>());
            item_class_[static_cast<std::size_t>(item)] =
                static_cast<std::uint16_t>(position - weights_.begin());
        }

        lb2_units_.assign(class_count_, 0);
        lb3_units_.assign(class_count_, 0);
        for (std::size_t item_class = 0; item_class < class_count_;
             ++item_class) {
            const std::int64_t weight = weights_[item_class];
            const std::int64_t capacity = capacity_;
            if (2 * weight > capacity) {
                lb2_units_[item_class] = 6;
            } else if (2 * weight == capacity) {
                lb2_units_[item_class] = 3;
            }
            if (3 * weight > 2 * capacity) {
                lb3_units_[item_class] = 6;
            } else if (3 * weight == 2 * capacity) {
                lb3_units_[item_class] = 4;
            } else if (3 * weight > capacity) {
                lb3_units_[item_class] = 3;
            } else if (3 * weight == capacity) {
                lb3_units_[item_class] = 2;
            }
        }

        memo_ = std::make_unique<ExactMultiplicityMemo>(
            class_count_, limits_.memo_entry_limit);
        const std::size_t depth_count = static_cast<std::size_t>(
            std::min(n_, limits_.maximum_item_count) + 1);
        state_counts_.assign(depth_count * class_count_, 0U);
        load_counts_.assign(depth_count * class_count_, 0U);
        candidate_classes_.assign(depth_count * class_count_, 0U);
        lookup_counts_.assign(class_count_, 0U);
        greedy_bin_loads_.assign(
            static_cast<std::size_t>(limits_.maximum_item_count), 0);

        constexpr int kMaximumDominanceWorkspace = 1'000'000;
        dominance_sum_limit_ = std::min(
            weights_.front(), kMaximumDominanceWorkspace);
        subset_reachable_.assign(
            static_cast<std::size_t>(dominance_sum_limit_ + 1), 0U);
    }

    [[nodiscard]] BinPackingBoundResult solve(
        const std::uint64_t* remaining_items,
        Deadline& global_deadline,
        double maximum_call_seconds) {
        if (maximum_call_seconds <= 0.0) {
            throw std::invalid_argument(
                "ordinary BINLB call budget must be positive");
        }
        BinPackingBoundResult result;
        int item_count = 0;
        std::int64_t total_weight = 0;
        int lb2_units = 0;
        int lb3_units = 0;
        std::uint16_t* root = state_row(0);
        std::fill(root, root + class_count_, 0U);
        for (int item = 0; item < n_; ++item) {
            if (!bit_is_set(remaining_items, item)) {
                continue;
            }
            ++item_count;
            const std::size_t item_class =
                item_class_[static_cast<std::size_t>(item)];
            ++root[item_class];
            total_weight += weights_[item_class];
            lb2_units += lb2_units_[item_class];
            lb3_units += lb3_units_[item_class];
        }
        result.lower_bound = lower_bound(total_weight, lb2_units, lb3_units);
        if (item_count > limits_.maximum_item_count) {
            result.item_limited = true;
            return result;
        }
        result.attempted = true;
        const auto start = Clock::now();
        call_end_ = std::min(
            global_deadline.end_time(),
            start + std::chrono::duration_cast<Clock::duration>(
                        std::chrono::duration<double>(
                            std::min(limits_.call_time_limit_seconds,
                                     maximum_call_seconds))));
        call_nodes_ = 0U;
        call_loads_ = 0U;
        call_memo_hits_ = 0U;
        aborted_ = false;
        timed_out_ = false;
        node_limited_ = false;
        load_limited_ = false;

        const NodeResult root_result = solve_state(
            0, total_weight, item_count, lb2_units, lb3_units);
        result.completed = root_result.completed;
        result.lower_bound = root_result.lower_bound;
        result.optimum = root_result.completed ? root_result.optimum : 0;
        result.timed_out = timed_out_;
        result.node_limited = node_limited_;
        result.load_limited = load_limited_;
        result.search_nodes = call_nodes_;
        result.nondominated_loads = call_loads_;
        result.memo_hits = call_memo_hits_;
        result.seconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }

    [[nodiscard]] bool lookup_exact(const std::uint64_t* remaining_items,
                                    int* optimum) {
        std::fill(lookup_counts_.begin(), lookup_counts_.end(), 0U);
        for (std::size_t block = 0; block < blocks_; ++block) {
            std::uint64_t value = remaining_items[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                const std::size_t item = block * 64U + bit;
                if (item < static_cast<std::size_t>(n_)) {
                    ++lookup_counts_[item_class_[item]];
                }
                value &= value - 1U;
            }
        }
        return memo_->find(lookup_counts_.data(), optimum);
    }

    [[nodiscard]] std::size_t blocks() const noexcept { return blocks_; }
    [[nodiscard]] std::uint64_t memo_entries() const noexcept {
        return memo_->size();
    }
    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        std::uint64_t bytes = memo_->memory_bytes();
        bytes = saturated_add(bytes, vector_memory_bytes(weights_));
        bytes = saturated_add(bytes, vector_memory_bytes(item_class_));
        bytes = saturated_add(bytes, vector_memory_bytes(lb2_units_));
        bytes = saturated_add(bytes, vector_memory_bytes(lb3_units_));
        bytes = saturated_add(bytes, vector_memory_bytes(state_counts_));
        bytes = saturated_add(bytes, vector_memory_bytes(load_counts_));
        bytes = saturated_add(bytes, vector_memory_bytes(candidate_classes_));
        bytes = saturated_add(bytes, vector_memory_bytes(lookup_counts_));
        bytes = saturated_add(bytes, vector_memory_bytes(greedy_bin_loads_));
        return saturated_add(bytes, vector_memory_bytes(subset_reachable_));
    }

private:
    [[nodiscard]] std::uint16_t* state_row(int depth) noexcept {
        return state_counts_.data() +
               static_cast<std::size_t>(depth) * class_count_;
    }
    [[nodiscard]] std::uint16_t* load_row(int depth) noexcept {
        return load_counts_.data() +
               static_cast<std::size_t>(depth) * class_count_;
    }
    [[nodiscard]] std::uint16_t* candidate_row(int depth) noexcept {
        return candidate_classes_.data() +
               static_cast<std::size_t>(depth) * class_count_;
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

    [[nodiscard]] int lower_bound(std::int64_t total_weight,
                                  int lb2_units,
                                  int lb3_units) const noexcept {
        if (total_weight == 0) {
            return 0;
        }
        return std::max({ceil_div(total_weight, capacity_),
                         ceil_div(lb2_units, 6),
                         ceil_div(lb3_units, 6)});
    }

    [[nodiscard]] int greedy_upper_bound(const std::uint16_t* state,
                                         int item_count) noexcept {
        std::fill(greedy_bin_loads_.begin(),
                  greedy_bin_loads_.begin() + item_count, 0);
        int bin_count = 0;
        for (std::size_t item_class = 0; item_class < class_count_;
             ++item_class) {
            const int weight = weights_[item_class];
            for (int copy = 0; copy < state[item_class]; ++copy) {
                int best_bin = -1;
                int best_load = -1;
                for (int bin = 0; bin < bin_count; ++bin) {
                    const int load =
                        greedy_bin_loads_[static_cast<std::size_t>(bin)];
                    if (load + weight <= capacity_ && load > best_load) {
                        best_bin = bin;
                        best_load = load;
                    }
                }
                if (best_bin < 0) {
                    best_bin = bin_count++;
                }
                greedy_bin_loads_[static_cast<std::size_t>(best_bin)] +=
                    weight;
            }
        }
        return bin_count;
    }

    [[nodiscard]] int minimum_improving_load(
        std::int64_t state_weight, int best) const noexcept {
        const std::int64_t threshold =
            state_weight - static_cast<std::int64_t>(best - 2) * capacity_;
        return static_cast<int>(std::max<std::int64_t>(0, threshold));
    }

    [[nodiscard]] bool load_is_nondominated(
        const std::uint16_t* state,
        const std::uint16_t* load,
        std::size_t anchor_class,
        int load_weight) {
        const int waste = capacity_ - load_weight;
        int largest_excluded = 0;
        for (std::size_t item_class = 0; item_class < class_count_;
             ++item_class) {
            if (state[item_class] <= load[item_class]) {
                continue;
            }
            const int weight = weights_[item_class];
            if (weight <= waste) {
                return false;
            }
            largest_excluded = std::max(largest_excluded, weight);
        }
        if (largest_excluded == 0 ||
            largest_excluded > dominance_sum_limit_) {
            return true;
        }

        std::fill(subset_reachable_.begin(),
                  subset_reachable_.begin() + largest_excluded + 1, 0U);
        subset_reachable_[0] = 1U;
        int reachable_limit = 0;
        for (std::size_t item_class = 0; item_class < class_count_;
             ++item_class) {
            int copies = load[item_class];
            if (item_class == anchor_class) {
                --copies;
            }
            const int weight = weights_[item_class];
            for (int copy = 0; copy < copies; ++copy) {
                const int next_limit = std::min(
                    largest_excluded, reachable_limit + weight);
                for (int sum = next_limit; sum >= weight; --sum) {
                    const unsigned char previous =
                        subset_reachable_[static_cast<std::size_t>(
                            sum - weight)];
                    if ((previous & 1U) != 0U) {
                        unsigned char flags = 1U;
                        if (sum > weight || (previous & 2U) != 0U) {
                            flags |= 2U;
                        }
                        subset_reachable_[static_cast<std::size_t>(sum)] |=
                            flags;
                    }
                }
                reachable_limit = next_limit;
            }
        }
        for (std::size_t item_class = 0; item_class < class_count_;
             ++item_class) {
            if (state[item_class] <= load[item_class]) {
                continue;
            }
            const int excluded_weight = weights_[item_class];
            const int first_sum = std::max(1, excluded_weight - waste);
            const int last_sum = std::min(excluded_weight, reachable_limit);
            for (int sum = first_sum; sum <= last_sum; ++sum) {
                const unsigned char flags =
                    subset_reachable_[static_cast<std::size_t>(sum)];
                if ((flags & 1U) != 0U &&
                    (sum < excluded_weight || (flags & 2U) != 0U)) {
                    return false;
                }
            }
        }
        return true;
    }

    bool enumerate_loads(int depth,
                         int position,
                         int candidate_count,
                         std::int64_t remaining_candidate_weight,
                         std::int64_t state_weight,
                         int state_item_count,
                         std::size_t anchor_class,
                         int load_weight,
                         int load_item_count,
                         int load_lb2_units,
                         int load_lb3_units,
                         int state_lb2_units,
                         int state_lb3_units,
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
            std::uint16_t* state = state_row(depth);
            std::uint16_t* load = load_row(depth);
            if (load_weight < minimum_improving_load(state_weight, best)) {
                return true;
            }
            const std::int64_t child_weight = state_weight - load_weight;
            const int child_lower = lower_bound(
                child_weight, state_lb2_units - load_lb2_units,
                state_lb3_units - load_lb3_units);
            if (1 + child_lower >= best) {
                return true;
            }
            if (!load_is_nondominated(
                    state, load, anchor_class, load_weight)) {
                return true;
            }
            if (state_load_count >=
                limits_.nondominated_load_limit_per_state) {
                load_limited_ = true;
                aborted_ = true;
                return false;
            }
            ++state_load_count;
            ++call_loads_;

            std::uint16_t* child = state_row(depth + 1);
            for (std::size_t item_class = 0; item_class < class_count_;
                 ++item_class) {
                child[item_class] = static_cast<std::uint16_t>(
                    state[item_class] - load[item_class]);
            }
            const NodeResult child_result = solve_state(
                depth + 1, child_weight,
                state_item_count - load_item_count,
                state_lb2_units - load_lb2_units,
                state_lb3_units - load_lb3_units);
            if (!child_result.completed) {
                return false;
            }
            best = std::min(best, 1 + child_result.optimum);
            return true;
        }

        const std::size_t item_class =
            candidate_row(depth)[static_cast<std::size_t>(position)];
        std::uint16_t* state = state_row(depth);
        std::uint16_t* load = load_row(depth);
        const int weight = weights_[item_class];
        const int base_count = load[item_class];
        const int available = state[item_class] - base_count;
        const int maximum_copies = std::min(
            available, (capacity_ - load_weight) / weight);
        const std::int64_t remaining_after =
            remaining_candidate_weight -
            static_cast<std::int64_t>(available) * weight;
        for (int copies = maximum_copies; copies >= 0; --copies) {
            load[item_class] = static_cast<std::uint16_t>(
                base_count + copies);
            if (!enumerate_loads(
                    depth, position + 1, candidate_count, remaining_after,
                    state_weight, state_item_count, anchor_class,
                    load_weight + copies * weight,
                    load_item_count + copies,
                    load_lb2_units + copies * lb2_units_[item_class],
                    load_lb3_units + copies * lb3_units_[item_class],
                    state_lb2_units, state_lb3_units, state_lower_bound,
                    best, state_load_count)) {
                load[item_class] = static_cast<std::uint16_t>(base_count);
                return false;
            }
            if (best == state_lower_bound) {
                load[item_class] = static_cast<std::uint16_t>(base_count);
                return true;
            }
        }
        load[item_class] = static_cast<std::uint16_t>(base_count);
        return true;
    }

    [[nodiscard]] NodeResult solve_state(int depth,
                                         std::int64_t total_weight,
                                         int item_count,
                                         int lb2_units,
                                         int lb3_units) {
        NodeResult result;
        result.lower_bound = lower_bound(total_weight, lb2_units, lb3_units);
        if (!consume_search_node()) {
            return result;
        }
        if (item_count == 0) {
            result.completed = true;
            return result;
        }
        std::uint16_t* state = state_row(depth);
        int memo_value = 0;
        if (memo_->find(state, &memo_value)) {
            ++call_memo_hits_;
            result.completed = true;
            result.lower_bound = memo_value;
            result.optimum = memo_value;
            return result;
        }

        int best = greedy_upper_bound(state, item_count);
        if (best == result.lower_bound) {
            result.completed = true;
            result.optimum = best;
            memo_->insert(state, best);
            return result;
        }

        std::size_t anchor_class = 0;
        while (anchor_class < class_count_ && state[anchor_class] == 0U) {
            ++anchor_class;
        }
        if (anchor_class == class_count_) {
            throw std::logic_error("nonempty ordinary BINLB state has no item");
        }
        std::uint16_t* load = load_row(depth);
        std::fill(load, load + class_count_, 0U);
        load[anchor_class] = 1U;
        const int anchor_weight = weights_[anchor_class];
        int candidate_count = 0;
        std::int64_t candidate_weight = 0;
        std::uint16_t* candidates = candidate_row(depth);
        for (std::size_t item_class = 0; item_class < class_count_;
             ++item_class) {
            const int available =
                state[item_class] - (item_class == anchor_class ? 1 : 0);
            if (available <= 0 ||
                anchor_weight + weights_[item_class] > capacity_) {
                continue;
            }
            candidates[static_cast<std::size_t>(candidate_count++)] =
                static_cast<std::uint16_t>(item_class);
            candidate_weight +=
                static_cast<std::int64_t>(available) * weights_[item_class];
        }

        std::uint64_t state_load_count = 0U;
        const bool enumerated = enumerate_loads(
            depth, 0, candidate_count, candidate_weight, total_weight,
            item_count, anchor_class, anchor_weight, 1,
            lb2_units_[anchor_class], lb3_units_[anchor_class],
            lb2_units, lb3_units, result.lower_bound, best,
            state_load_count);
        if (!enumerated || aborted_) {
            return result;
        }
        result.completed = true;
        result.optimum = best;
        memo_->insert(state, best);
        return result;
    }

    const Instance& instance_;
    BinPackingBoundLimits limits_;
    int n_ = 0;
    int capacity_ = 0;
    std::size_t blocks_ = 0;
    std::size_t class_count_ = 0;
    std::vector<int> weights_;
    std::vector<std::uint16_t> item_class_;
    std::vector<int> lb2_units_;
    std::vector<int> lb3_units_;
    std::unique_ptr<ExactMultiplicityMemo> memo_;

    std::vector<std::uint16_t> state_counts_;
    std::vector<std::uint16_t> load_counts_;
    std::vector<std::uint16_t> candidate_classes_;
    std::vector<std::uint16_t> lookup_counts_;
    std::vector<int> greedy_bin_loads_;
    int dominance_sum_limit_ = 0;
    std::vector<unsigned char> subset_reachable_;

    Clock::time_point call_end_{};
    std::uint64_t call_nodes_ = 0U;
    std::uint64_t call_loads_ = 0U;
    std::uint64_t call_memo_hits_ = 0U;
    bool aborted_ = false;
    bool timed_out_ = false;
    bool node_limited_ = false;
    bool load_limited_ = false;
};

BinPackingBound::BinPackingBound(const Instance& instance,
                                 BinPackingBoundLimits limits)
    : impl_(std::make_unique<Impl>(instance, limits)) {}

BinPackingBound::~BinPackingBound() = default;
BinPackingBound::BinPackingBound(BinPackingBound&&) noexcept = default;
BinPackingBound& BinPackingBound::operator=(BinPackingBound&&) noexcept =
    default;

BinPackingBoundResult BinPackingBound::solve(
    const std::uint64_t* remaining_items,
    Deadline& global_deadline) {
    return solve(remaining_items, global_deadline,
                 std::numeric_limits<double>::max());
}

BinPackingBoundResult BinPackingBound::solve(
    const std::uint64_t* remaining_items,
    Deadline& global_deadline,
    double maximum_call_seconds) {
    if (remaining_items == nullptr) {
        throw std::invalid_argument("null ordinary BINLB item set");
    }
    return impl_->solve(
        remaining_items, global_deadline, maximum_call_seconds);
}

bool BinPackingBound::lookup_exact(const std::uint64_t* remaining_items,
                                   int* optimum) {
    if (remaining_items == nullptr || optimum == nullptr) {
        throw std::invalid_argument("null ordinary BINLB memo query");
    }
    return impl_->lookup_exact(remaining_items, optimum);
}

std::size_t BinPackingBound::bit_block_count() const noexcept {
    return impl_->blocks();
}

std::uint64_t BinPackingBound::memo_entry_count() const noexcept {
    return impl_->memo_entries();
}

std::uint64_t BinPackingBound::memory_bytes() const noexcept {
    return impl_->memory_bytes();
}

}
