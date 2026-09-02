#include "precpack/bbr.hpp"

#include "precpack/algorithms.hpp"
#include "precpack/bin_packing_bound.hpp"
#include "precpack/dff.hpp"
#include "precpack/solver_profile.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace precpack {
namespace {

using Clock = std::chrono::steady_clock;

constexpr std::uint32_t kInvalidState =
    std::numeric_limits<std::uint32_t>::max();
constexpr std::uint64_t kMaximumStateCount = kInvalidState;

[[nodiscard]] int ceil_div_i64(std::int64_t numerator,
                               std::int64_t denominator,
                               const std::source_location location =
                                   std::source_location::current()) {
    if (numerator < 0 || denominator <= 0) {
        throw std::logic_error(
            "invalid BBR integer division at " +
            std::string(location.file_name()) + ":" +
            std::to_string(location.line()) +
            " (numerator=" + std::to_string(numerator) +
            ", denominator=" + std::to_string(denominator) + ")");
    }
    return static_cast<int>((numerator + denominator - 1) / denominator);
}

[[nodiscard]] std::uint64_t mix64(std::uint64_t value) noexcept {
    value ^= value >> 30U;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27U;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return value;
}

[[nodiscard]] std::uint64_t saturated_add(std::uint64_t lhs,
                                          std::uint64_t rhs) noexcept {
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    return lhs > maximum - rhs ? maximum : lhs + rhs;
}

[[nodiscard]] std::uint64_t saturated_multiply(std::size_t count,
                                               std::size_t width) noexcept {
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    if (width != 0U && count > maximum / width) {
        return maximum;
    }
    return static_cast<std::uint64_t>(count) * width;
}

template <class T>
[[nodiscard]] std::uint64_t vector_memory_bytes(
    const std::vector<T>& values) noexcept {
    return saturated_multiply(values.capacity(), sizeof(T));
}

[[nodiscard]] std::uint64_t hash_words(const std::uint64_t* words,
                                       std::size_t count) noexcept {
    if (count == 2U) {
        std::uint64_t hash = 0x6a09e667f3bcc909ULL ^ mix64(2U);
        hash = mix64(hash ^ mix64(words[0] + 0x9e3779b97f4a7c15ULL));
        return mix64(hash ^
                     mix64(words[1] + 2U * 0x9e3779b97f4a7c15ULL));
    }
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
            1U) != 0U;
}

void set_bit(std::uint64_t* words, int item) noexcept {
    words[static_cast<std::size_t>(item) >> 6U] |=
        std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
}

void clear_bit(std::uint64_t* words, int item) noexcept {
    words[static_cast<std::size_t>(item) >> 6U] &=
        ~(std::uint64_t{1} << (static_cast<unsigned>(item) & 63U));
}

[[nodiscard]] bool mask_subset(const std::uint64_t* subset,
                               const std::uint64_t* superset,
                               std::size_t blocks) noexcept {
    if (blocks == 2U) {
        return (subset[0] & ~superset[0]) == 0U &&
               (subset[1] & ~superset[1]) == 0U;
    }
    for (std::size_t block = 0; block < blocks; ++block) {
        if ((subset[block] & ~superset[block]) != 0U) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool masks_equal(const std::uint64_t* lhs,
                               const std::uint64_t* rhs,
                               std::size_t blocks) noexcept {
    if (blocks == 2U) {
        return lhs[0] == rhs[0] && lhs[1] == rhs[1];
    }
    for (std::size_t block = 0; block < blocks; ++block) {
        if (lhs[block] != rhs[block]) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] int mask_popcount(const std::uint64_t* words,
                                std::size_t blocks) noexcept {
    if (blocks == 2U) {
        return std::popcount(words[0]) + std::popcount(words[1]);
    }
    int count = 0;
    for (std::size_t block = 0; block < blocks; ++block) {
        count += std::popcount(words[block]);
    }
    return count;
}

struct DffCandidate {
    std::int64_t capacity = 1;
    std::vector<std::int64_t> contribution;
    std::int64_t total = 0;
    int root_bound = 0;
};

struct BbrPrecomputed {
    explicit BbrPrecomputed(const Instance& input,
                            bool enable_item_dominance,
                            bool enable_generalized_item_dominance,
                            bool enable_complete_dff,
                            int dff_transform_limit,
                            const Deadline& deadline)
        : instance(input),
          n(input.size()),
          capacity(input.capacity),
          blocks((static_cast<std::size_t>(n) + 63U) / 64U),
          all_mask(blocks, std::numeric_limits<std::uint64_t>::max()),
          item_hash(static_cast<std::size_t>(n), 0U) {
        const unsigned tail = static_cast<unsigned>(n) & 63U;
        if (tail != 0U) {
            all_mask.back() = (std::uint64_t{1} << tail) - 1U;
        }
        for (int item = 0; item < n; ++item) {
            item_hash[static_cast<std::size_t>(item)] = mix64(
                0x243f6a8885a308d3ULL +
                0x9e3779b97f4a7c15ULL * static_cast<std::uint64_t>(item + 1));
        }
        for (const Arc& arc : instance.arcs) {
            maximum_separation = std::max(maximum_separation, arc.separation);
            salbp_semantics = salbp_semantics && arc.separation == 0;
            bppp_semantics = bppp_semantics && arc.separation == 1;
        }
        cooldown_levels = std::max(0, maximum_separation - 1);
        key_words = blocks * static_cast<std::size_t>(1 + cooldown_levels);
        build_weight_byte_lookup();
        build_direct_masks();
        build_closure_masks();
        build_order_and_jackson_pairs(enable_item_dominance,
                                      enable_generalized_item_dominance,
                                      deadline);
        build_generalized_dominator_masks(
            enable_generalized_item_dominance);
        complete_dff_used = enable_complete_dff &&
                            (salbp_semantics || bppp_semantics);
        build_dff_tables(complete_dff_used, dff_transform_limit);
    }

    [[nodiscard]] const std::uint64_t* pred_zero_row(int item) const noexcept {
        return pred_zero.data() + static_cast<std::size_t>(item) * blocks;
    }

    [[nodiscard]] const std::uint64_t* pred_positive_row(int item) const noexcept {
        return pred_positive.data() + static_cast<std::size_t>(item) * blocks;
    }

    [[nodiscard]] const std::uint64_t* predecessor_closure_row(
        int item) const noexcept {
        return instance.reaching_row(item);
    }

    [[nodiscard]] const std::uint64_t* successor_closure_row(
        int item) const noexcept {
        return instance.reachable_row(item);
    }

    [[nodiscard]] const std::uint64_t* zero_successor_row(
        int item) const noexcept {
        return zero_successor_mask.data() +
               static_cast<std::size_t>(item) * blocks;
    }

    [[nodiscard]] const std::uint64_t* successors_above_row(
        int level,
        int item) const noexcept {
        return successors_above.data() +
               (static_cast<std::size_t>(level) * n + item) * blocks;
    }

    [[nodiscard]] const std::uint64_t* generalized_dominator_row(
        int item) const noexcept {
        return generalized_dominator_mask.data() +
               static_cast<std::size_t>(item) * blocks;
    }

    [[nodiscard]] const std::int64_t* transformed_weights_for_item(
        int item) const noexcept {
        return dff_contribution.data() +
               static_cast<std::size_t>(item) * dff_capacity.size();
    }

    [[nodiscard]] std::uint64_t assigned_set_hash(
        const std::uint64_t* assigned) const noexcept {
        std::uint64_t hash = 0x13198a2e03707344ULL;
        for_each_set_bit(assigned, [&](int item) {
            hash ^= item_hash[static_cast<std::size_t>(item)];
        });
        return hash;
    }

    template <class Callback>
    void for_each_set_bit(const std::uint64_t* mask, Callback&& callback) const {
        for (std::size_t block = 0; block < blocks; ++block) {
            std::uint64_t value = mask[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                const int item = static_cast<int>(block * 64U + bit);
                callback(item);
                value &= value - 1U;
            }
        }
    }

    [[nodiscard]] std::int64_t masked_weight_sum(
        const std::uint64_t* mask,
        const std::uint64_t* assigned) const noexcept {
        std::int64_t sum = 0;
        for (std::size_t segment = 0; segment < weight_byte_segments;
             ++segment) {
            const std::size_t block = segment >> 3U;
            const unsigned shift = static_cast<unsigned>(segment & 7U) * 8U;
            const unsigned value = static_cast<unsigned>(
                (mask[block] & ~assigned[block]) >> shift) & 0xffU;
            sum += weight_byte_lookup[segment * 256U + value];
        }
        return sum;
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        std::uint64_t bytes = 0;
        const auto add = [&](std::uint64_t amount) {
            bytes = saturated_add(bytes, amount);
        };
        add(vector_memory_bytes(all_mask));
        add(vector_memory_bytes(item_hash));
        add(vector_memory_bytes(weight_byte_lookup));
        add(vector_memory_bytes(pred_zero));
        add(vector_memory_bytes(pred_positive));
        add(vector_memory_bytes(zero_successor_mask));
        add(vector_memory_bytes(successors_above));
        add(vector_memory_bytes(predecessor_closure_capacity_bins));
        add(vector_memory_bytes(successor_closure_capacity_bins));
        add(vector_memory_bytes(closure_bound_order));
        add(vector_memory_bytes(closure_bound_static_potential));
        add(vector_memory_bytes(zero_successor_offset));
        add(vector_memory_bytes(zero_successors));
        add(vector_memory_bytes(branch_order));
        add(vector_memory_bytes(jackson_dominators));
        for (const auto& values : jackson_dominators) {
            add(vector_memory_bytes(values));
        }
        add(vector_memory_bytes(generalized_dominators));
        for (const auto& values : generalized_dominators) {
            add(vector_memory_bytes(values));
        }
        add(vector_memory_bytes(generalized_dominator_mask));
        add(vector_memory_bytes(has_successor));
        add(vector_memory_bytes(has_successor_mask));
        add(vector_memory_bytes(dff_capacity));
        add(vector_memory_bytes(dff_contribution));
        add(vector_memory_bytes(dff_total));
        return bytes;
    }

    const Instance& instance;
    int n = 0;
    int capacity = 0;
    std::size_t blocks = 0;
    int maximum_separation = 0;
    int cooldown_levels = 0;
    std::size_t key_words = 0;
    bool salbp_semantics = true;
    bool bppp_semantics = true;
    std::vector<std::uint64_t> all_mask;
    std::vector<std::uint64_t> item_hash;
    std::size_t weight_byte_segments = 0;
    std::vector<std::int64_t> weight_byte_lookup;
    std::vector<std::uint64_t> pred_zero;
    std::vector<std::uint64_t> pred_positive;
    std::vector<std::uint64_t> zero_successor_mask;
    std::vector<std::uint64_t> successors_above;
    std::vector<int> predecessor_closure_capacity_bins;
    std::vector<int> successor_closure_capacity_bins;
    std::vector<int> closure_bound_order;
    std::vector<std::int64_t> closure_bound_static_potential;
    std::vector<int> zero_successor_offset;
    std::vector<int> zero_successors;
    std::vector<int> branch_order;
    std::vector<std::vector<int>> jackson_dominators;
    std::vector<std::vector<int>> generalized_dominators;
    std::vector<std::uint64_t> generalized_dominator_mask;
    std::uint64_t generalized_dominance_search_nodes = 0;
    double generalized_dominance_seconds = 0.0;
    std::vector<unsigned char> has_successor;
    std::vector<std::uint64_t> has_successor_mask;
    std::vector<std::int64_t> dff_capacity;
    std::vector<std::int64_t> dff_contribution;
    std::vector<std::int64_t> dff_total;
    int root_dff_bound = 0;
    bool complete_dff_used = false;

private:
    void build_weight_byte_lookup() {
        weight_byte_segments = (static_cast<std::size_t>(n) + 7U) / 8U;
        weight_byte_lookup.assign(weight_byte_segments * 256U, 0);
        for (std::size_t segment = 0; segment < weight_byte_segments;
             ++segment) {
            std::int64_t* table =
                weight_byte_lookup.data() + segment * 256U;
            for (unsigned value = 1; value < 256U; ++value) {
                const unsigned previous = value & (value - 1U);
                const unsigned bit = std::countr_zero(value);
                const std::size_t item = segment * 8U + bit;
                table[value] = table[previous];
                if (item < static_cast<std::size_t>(n)) {
                    table[value] += instance.items[item].weight;
                }
            }
        }
    }

    void build_generalized_dominator_masks(bool enabled) {
        if (!enabled || !bppp_semantics) {
            return;
        }
        generalized_dominator_mask.assign(
            static_cast<std::size_t>(n) * blocks, 0U);
        for (int dominated = 0; dominated < n; ++dominated) {
            std::uint64_t* row = generalized_dominator_mask.data() +
                static_cast<std::size_t>(dominated) * blocks;
            for (const int dominator :
                 generalized_dominators[static_cast<std::size_t>(dominated)]) {
                set_bit(row, dominator);
            }
        }
    }

    void build_direct_masks() {
        pred_zero.assign(static_cast<std::size_t>(n) * blocks, 0U);
        pred_positive.assign(static_cast<std::size_t>(n) * blocks, 0U);
        zero_successor_mask.assign(static_cast<std::size_t>(n) * blocks, 0U);
        successors_above.assign(
            static_cast<std::size_t>(cooldown_levels) * n * blocks, 0U);
        std::vector<int> zero_degree(static_cast<std::size_t>(n), 0);
        for (const Arc& arc : instance.arcs) {
            std::uint64_t* predecessor_row =
                (arc.separation == 0 ? pred_zero.data() : pred_positive.data()) +
                static_cast<std::size_t>(arc.to) * blocks;
            set_bit(predecessor_row, arc.from);
            if (arc.separation == 0) {
                ++zero_degree[static_cast<std::size_t>(arc.from)];
                set_bit(zero_successor_mask.data() +
                            static_cast<std::size_t>(arc.from) * blocks,
                        arc.to);
            }
            for (int level = 1; level <= cooldown_levels; ++level) {
                if (arc.separation > level) {
                    std::uint64_t* row = successors_above.data() +
                        (static_cast<std::size_t>(level - 1) * n + arc.from) *
                            blocks;
                    set_bit(row, arc.to);
                }
            }
        }
        zero_successor_offset.assign(static_cast<std::size_t>(n + 1), 0);
        for (int item = 0; item < n; ++item) {
            zero_successor_offset[static_cast<std::size_t>(item + 1)] =
                zero_successor_offset[static_cast<std::size_t>(item)] +
                zero_degree[static_cast<std::size_t>(item)];
        }
        zero_successors.resize(
            static_cast<std::size_t>(zero_successor_offset.back()));
        std::vector<int> cursor = zero_successor_offset;
        for (const Arc& arc : instance.arcs) {
            if (arc.separation == 0) {
                zero_successors[static_cast<std::size_t>(
                    cursor[static_cast<std::size_t>(arc.from)]++)] = arc.to;
            }
        }
    }

    void build_closure_masks() {
        predecessor_closure_capacity_bins.resize(static_cast<std::size_t>(n));
        successor_closure_capacity_bins.resize(static_cast<std::size_t>(n));
        for (int item = 0; item < n; ++item) {
            std::int64_t predecessor_weight =
                instance.items[static_cast<std::size_t>(item)].weight;
            std::int64_t successor_weight = predecessor_weight;
            for_each_set_bit(predecessor_closure_row(item), [&](int predecessor) {
                predecessor_weight += instance.items[
                    static_cast<std::size_t>(predecessor)].weight;
            });
            for_each_set_bit(successor_closure_row(item), [&](int successor) {
                successor_weight += instance.items[
                    static_cast<std::size_t>(successor)].weight;
            });
            predecessor_closure_capacity_bins[static_cast<std::size_t>(item)] =
                ceil_div_i64(predecessor_weight, capacity);
            successor_closure_capacity_bins[static_cast<std::size_t>(item)] =
                ceil_div_i64(successor_weight, capacity);
        }

        std::vector<std::int64_t> global_earliest(static_cast<std::size_t>(n),
                                                  0);
        std::vector<std::int64_t> global_tail(static_cast<std::size_t>(n), 0);
        for (const int item : instance.topological_order) {
            for (const auto& [predecessor, separation] :
                 instance.predecessor_arcs[static_cast<std::size_t>(item)]) {
                global_earliest[static_cast<std::size_t>(item)] = std::max(
                    global_earliest[static_cast<std::size_t>(item)],
                    global_earliest[static_cast<std::size_t>(predecessor)] +
                        separation);
            }
        }
        for (auto order = instance.topological_order.rbegin();
             order != instance.topological_order.rend(); ++order) {
            const int item = *order;
            for (const auto& [successor, separation] :
                 instance.successor_arcs[static_cast<std::size_t>(item)]) {
                global_tail[static_cast<std::size_t>(item)] = std::max(
                    global_tail[static_cast<std::size_t>(item)],
                    static_cast<std::int64_t>(separation) +
                        global_tail[static_cast<std::size_t>(successor)]);
            }
        }
        closure_bound_static_potential.resize(static_cast<std::size_t>(n));
        closure_bound_order.resize(static_cast<std::size_t>(n));
        std::iota(closure_bound_order.begin(), closure_bound_order.end(), 0);
        for (int item = 0; item < n; ++item) {
            const std::size_t index = static_cast<std::size_t>(item);
            const std::int64_t prefix_upper = std::max(
                static_cast<std::int64_t>(
                    predecessor_closure_capacity_bins[index]),
                global_earliest[index] + 1);
            const std::int64_t suffix_upper = std::max(
                static_cast<std::int64_t>(
                    successor_closure_capacity_bins[index]),
                global_tail[index] + 1);
            closure_bound_static_potential[index] =
                prefix_upper + suffix_upper - 1;
        }
        std::stable_sort(
            closure_bound_order.begin(), closure_bound_order.end(),
            [&](int lhs, int rhs) {
                return std::pair{
                           closure_bound_static_potential[
                               static_cast<std::size_t>(lhs)],
                           -lhs} >
                       std::pair{
                           closure_bound_static_potential[
                               static_cast<std::size_t>(rhs)],
                           -rhs};
            });
    }

    void build_order_and_jackson_pairs(
        bool enable_item_dominance,
        bool enable_generalized_item_dominance,
        const Deadline& deadline) {
        jackson_dominators.assign(static_cast<std::size_t>(n), {});
        generalized_dominators.assign(static_cast<std::size_t>(n), {});
        has_successor.assign(static_cast<std::size_t>(n), 0U);
        has_successor_mask.assign(blocks, 0U);
        std::vector<int> dominance_count(static_cast<std::size_t>(n), 0);
        std::vector<std::int64_t> positional_weight(static_cast<std::size_t>(n), 0);
        std::vector<int> successor_count(static_cast<std::size_t>(n), 0);
        for (int item = 0; item < n; ++item) {
            const std::uint64_t* successors = successor_closure_row(item);
            successor_count[static_cast<std::size_t>(item)] =
                mask_popcount(successors, blocks);
            if (successor_count[static_cast<std::size_t>(item)] > 0) {
                has_successor[static_cast<std::size_t>(item)] = 1U;
                set_bit(has_successor_mask.data(), item);
            }
            positional_weight[static_cast<std::size_t>(item)] =
                instance.items[static_cast<std::size_t>(item)].weight;
            for_each_set_bit(successors, [&](int successor) {
                positional_weight[static_cast<std::size_t>(item)] +=
                    instance.items[static_cast<std::size_t>(successor)].weight;
            });
        }
        const bool generalized_separations =
            !salbp_semantics && !bppp_semantics;
        std::vector<std::uint64_t> direct_successors;
        std::vector<std::uint64_t> positive_successors;
        if (generalized_separations && enable_item_dominance) {
            direct_successors.assign(static_cast<std::size_t>(n) * blocks, 0U);
            positive_successors.assign(static_cast<std::size_t>(n) * blocks,
                                       0U);
            for (const Arc& arc : instance.arcs) {
                set_bit(direct_successors.data() +
                            static_cast<std::size_t>(arc.from) * blocks,
                        arc.to);
                if (arc.separation > 0) {
                    set_bit(positive_successors.data() +
                                static_cast<std::size_t>(arc.from) * blocks,
                            arc.to);
                }
            }
        }
        const auto separation_profile_covers =
            [&](int dominator, int dominated) noexcept {
                const std::size_t dominator_offset =
                    static_cast<std::size_t>(dominator) * blocks;
                const std::size_t dominated_offset =
                    static_cast<std::size_t>(dominated) * blocks;
                if (!mask_subset(direct_successors.data() + dominated_offset,
                                 direct_successors.data() + dominator_offset,
                                 blocks) ||
                    !mask_subset(positive_successors.data() + dominated_offset,
                                 positive_successors.data() + dominator_offset,
                                 blocks)) {
                    return false;
                }
                for (int level = 0; level < cooldown_levels; ++level) {
                    if (!mask_subset(successors_above_row(level, dominated),
                                     successors_above_row(level, dominator),
                                     blocks)) {
                        return false;
                    }
                }
                return true;
            };
        if (enable_item_dominance) {
            for (int dominated = 0; dominated < n; ++dominated) {
                const std::uint64_t* dominated_successors =
                    successor_closure_row(dominated);
                for (int dominator = 0; dominator < n; ++dominator) {
                    if (dominator == dominated ||
                        instance.items[static_cast<std::size_t>(dominator)].weight <
                            instance.items[static_cast<std::size_t>(dominated)].weight) {
                        continue;
                    }
                    if (bit_is_set(successor_closure_row(dominator), dominated) ||
                        bit_is_set(dominated_successors, dominator)) {
                        continue;
                    }
                    const bool profile_covers = generalized_separations
                        ? separation_profile_covers(dominator, dominated)
                        : mask_subset(dominated_successors,
                                      successor_closure_row(dominator), blocks);
                    if (!profile_covers) {
                        continue;
                    }
                    const bool equal_weight =
                        instance.items[static_cast<std::size_t>(dominator)].weight ==
                        instance.items[static_cast<std::size_t>(dominated)].weight;
                    const bool equal_successors = generalized_separations
                        ? separation_profile_covers(dominated, dominator)
                        : masks_equal(dominated_successors,
                                      successor_closure_row(dominator), blocks);
                    if (equal_weight && equal_successors && dominator > dominated) {
                        continue;
                    }
                    jackson_dominators[static_cast<std::size_t>(dominated)]
                        .push_back(dominator);
                    if (salbp_semantics) {
                        ++dominance_count[
                            static_cast<std::size_t>(dominator)];
                    }
                }
            }
        }
        if (enable_generalized_item_dominance &&
            (salbp_semantics || bppp_semantics)) {
            const auto started = Clock::now();
            build_generalized_item_dominance(dominance_count, deadline);
            generalized_dominance_seconds =
                std::chrono::duration<double>(Clock::now() - started).count();
        }
        branch_order.resize(static_cast<std::size_t>(n));
        std::iota(branch_order.begin(), branch_order.end(), 0);
        std::stable_sort(branch_order.begin(), branch_order.end(),
                         [&](int lhs, int rhs) {
            return std::tuple{
                       dominance_count[static_cast<std::size_t>(lhs)],
                       positional_weight[static_cast<std::size_t>(lhs)],
                       instance.items[static_cast<std::size_t>(lhs)].weight,
                       successor_count[static_cast<std::size_t>(lhs)], -lhs} >
                   std::tuple{
                       dominance_count[static_cast<std::size_t>(rhs)],
                       positional_weight[static_cast<std::size_t>(rhs)],
                       instance.items[static_cast<std::size_t>(rhs)].weight,
                       successor_count[static_cast<std::size_t>(rhs)], -rhs};
                         });
    }

    void build_generalized_item_dominance(
        std::vector<int>& dominance_count,
        const Deadline& deadline) {
        if (deadline.expired()) {
            return;
        }

        constexpr std::size_t kArcPollMask = 1023U;
        constexpr std::size_t kStructuralPollMask = 63U;
        constexpr std::size_t kCandidatePollMask = 31U;

        std::vector<std::uint64_t> direct_out(
            static_cast<std::size_t>(n) * blocks, 0U);
        for (std::size_t arc_index = 0; arc_index < instance.arcs.size();
             ++arc_index) {
            if ((arc_index & kArcPollMask) == kArcPollMask &&
                deadline.expired()) {
                return;
            }
            const Arc& arc = instance.arcs[arc_index];
            set_bit(direct_out.data() +
                        static_cast<std::size_t>(arc.from) * blocks,
                    arc.to);
        }
        const auto edge = [&](int from, int to) noexcept {
            return bit_is_set(direct_out.data() +
                                  static_cast<std::size_t>(from) * blocks,
                              to);
        };
        const auto already_basic = [&](int dominator, int dominated) {
            const auto& values = jackson_dominators[
                static_cast<std::size_t>(dominated)];
            return std::find(values.begin(), values.end(), dominator) !=
                   values.end();
        };

        constexpr std::uint64_t kTotalSearchNodeLimit = 2'000'000U;
        constexpr std::uint64_t kPairSearchNodeLimit = 100'000U;

        for (int dominated = 0; dominated < n; ++dominated) {
            if (deadline.expired()) {
                return;
            }
            std::vector<int> source_vertices{dominated};
            for_each_set_bit(successor_closure_row(dominated),
                             [&](int item) { source_vertices.push_back(item); });
            for (int dominator = 0; dominator < n; ++dominator) {
                if (generalized_dominance_search_nodes >=
                        kTotalSearchNodeLimit ||
                    ((static_cast<std::size_t>(dominator) &
                      kStructuralPollMask) == kStructuralPollMask &&
                     deadline.expired())) {
                    return;
                }
                if (dominator == dominated ||
                    instance.items[static_cast<std::size_t>(dominator)].weight <
                        instance.items[static_cast<std::size_t>(dominated)].weight ||
                    bit_is_set(successor_closure_row(dominator), dominated) ||
                    bit_is_set(successor_closure_row(dominated), dominator)) {
                    continue;
                }
                if (already_basic(dominator, dominated)) {
                    generalized_dominators[
                        static_cast<std::size_t>(dominated)]
                        .push_back(dominator);
                    continue;
                }
                const bool equal_root_weight =
                    instance.items[static_cast<std::size_t>(dominator)].weight ==
                    instance.items[static_cast<std::size_t>(dominated)].weight;
                if (equal_root_weight && dominator > dominated) {
                    continue;
                }

                std::vector<int> target_vertices{dominator};
                for_each_set_bit(successor_closure_row(dominator),
                                 [&](int item) {
                                     target_vertices.push_back(item);
                                 });
                if (source_vertices.size() > target_vertices.size()) {
                    continue;
                }

                std::vector<int> source_in(static_cast<std::size_t>(n), 0);
                std::vector<int> source_out(static_cast<std::size_t>(n), 0);
                std::vector<int> target_in(static_cast<std::size_t>(n), 0);
                std::vector<int> target_out(static_cast<std::size_t>(n), 0);
                for (std::size_t lhs_index = 0;
                     lhs_index < source_vertices.size(); ++lhs_index) {
                    if ((lhs_index & kStructuralPollMask) ==
                            kStructuralPollMask &&
                        deadline.expired()) {
                        return;
                    }
                    const int lhs = source_vertices[lhs_index];
                    for (const int rhs : source_vertices) {
                        if (edge(lhs, rhs)) {
                            ++source_out[static_cast<std::size_t>(lhs)];
                            ++source_in[static_cast<std::size_t>(rhs)];
                        }
                    }
                }
                for (std::size_t lhs_index = 0;
                     lhs_index < target_vertices.size(); ++lhs_index) {
                    if ((lhs_index & kStructuralPollMask) ==
                            kStructuralPollMask &&
                        deadline.expired()) {
                        return;
                    }
                    const int lhs = target_vertices[lhs_index];
                    for (const int rhs : target_vertices) {
                        if (edge(lhs, rhs)) {
                            ++target_out[static_cast<std::size_t>(lhs)];
                            ++target_in[static_cast<std::size_t>(rhs)];
                        }
                    }
                }
                if (source_out[static_cast<std::size_t>(dominated)] >
                        target_out[static_cast<std::size_t>(dominator)] ||
                    source_in[static_cast<std::size_t>(dominated)] >
                        target_in[static_cast<std::size_t>(dominator)]) {
                    continue;
                }

                std::vector<std::vector<int>> candidates(
                    static_cast<std::size_t>(n));
                bool candidate_failure = false;
                for (std::size_t source_index = 0;
                     source_index < source_vertices.size(); ++source_index) {
                    if ((source_index & kStructuralPollMask) ==
                            kStructuralPollMask &&
                        deadline.expired()) {
                        return;
                    }
                    const int source = source_vertices[source_index];
                    if (source == dominated) {
                        candidates[static_cast<std::size_t>(source)] = {
                            dominator};
                        continue;
                    }
                    auto& values = candidates[static_cast<std::size_t>(source)];
                    for (std::size_t target_index = 0;
                         target_index < target_vertices.size(); ++target_index) {
                        if ((target_index & kCandidatePollMask) ==
                                kCandidatePollMask &&
                            deadline.expired()) {
                            return;
                        }
                        const int target = target_vertices[target_index];
                        if (target == dominator ||
                            instance.items[static_cast<std::size_t>(source)].weight !=
                                instance.items[static_cast<std::size_t>(target)].weight ||
                            source_in[static_cast<std::size_t>(source)] >
                                target_in[static_cast<std::size_t>(target)] ||
                            source_out[static_cast<std::size_t>(source)] >
                                target_out[static_cast<std::size_t>(target)] ||
                            (edge(dominated, source) &&
                             !edge(dominator, target))) {
                            continue;
                        }
                        values.push_back(target);
                    }
                    if (values.empty()) {
                        candidate_failure = true;
                        break;
                    }
                }
                if (candidate_failure) {
                    continue;
                }

                std::vector<int> mapping(static_cast<std::size_t>(n), -1);
                std::vector<unsigned char> target_used(
                    static_cast<std::size_t>(n), 0U);
                mapping[static_cast<std::size_t>(dominated)] = dominator;
                target_used[static_cast<std::size_t>(dominator)] = 1U;
                std::uint64_t pair_nodes = 0;
                bool aborted = false;
                bool deadline_exhausted = false;
                const auto search = [&](auto&& self, std::size_t mapped_count)
                    -> bool {
                    ++pair_nodes;
                    ++generalized_dominance_search_nodes;
                    if (pair_nodes > kPairSearchNodeLimit ||
                        generalized_dominance_search_nodes >
                            kTotalSearchNodeLimit) {
                        aborted = true;
                        return false;
                    }
                    if ((pair_nodes & 255U) == 0U && deadline.expired()) {
                        aborted = true;
                        deadline_exhausted = true;
                        return false;
                    }
                    if (mapped_count == source_vertices.size()) {
                        return true;
                    }
                    int selected_source = -1;
                    std::size_t best_count =
                        std::numeric_limits<std::size_t>::max();
                    for (const int source : source_vertices) {
                        if (mapping[static_cast<std::size_t>(source)] >= 0) {
                            continue;
                        }
                        std::size_t viable = 0;
                        const auto& source_candidates = candidates[
                            static_cast<std::size_t>(source)];
                        for (std::size_t target_index = 0;
                             target_index < source_candidates.size();
                             ++target_index) {
                            if ((target_index & kCandidatePollMask) ==
                                    kCandidatePollMask &&
                                deadline.expired()) {
                                aborted = true;
                                deadline_exhausted = true;
                                return false;
                            }
                            const int target = source_candidates[target_index];
                            if (target_used[static_cast<std::size_t>(target)] != 0U) {
                                continue;
                            }
                            bool compatible = true;
                            for (const int other : source_vertices) {
                                const int mapped = mapping[
                                    static_cast<std::size_t>(other)];
                                if (mapped < 0) {
                                    continue;
                                }
                                if ((edge(source, other) &&
                                     !edge(target, mapped)) ||
                                    (edge(other, source) &&
                                     !edge(mapped, target))) {
                                    compatible = false;
                                    break;
                                }
                            }
                            viable += compatible ? 1U : 0U;
                        }
                        if (viable < best_count) {
                            best_count = viable;
                            selected_source = source;
                            if (viable == 0U) {
                                break;
                            }
                        }
                    }
                    if (selected_source < 0 || best_count == 0U) {
                        return false;
                    }
                    const auto& selected_candidates = candidates[
                        static_cast<std::size_t>(selected_source)];
                    for (std::size_t target_index = 0;
                         target_index < selected_candidates.size();
                         ++target_index) {
                        if ((target_index & kCandidatePollMask) ==
                                kCandidatePollMask &&
                            deadline.expired()) {
                            aborted = true;
                            deadline_exhausted = true;
                            return false;
                        }
                        const int target = selected_candidates[target_index];
                        if (target_used[static_cast<std::size_t>(target)] != 0U) {
                            continue;
                        }
                        bool compatible = true;
                        for (const int other : source_vertices) {
                            const int mapped = mapping[
                                static_cast<std::size_t>(other)];
                            if (mapped < 0) {
                                continue;
                            }
                            if ((edge(selected_source, other) &&
                                 !edge(target, mapped)) ||
                                (edge(other, selected_source) &&
                                 !edge(mapped, target))) {
                                compatible = false;
                                break;
                            }
                        }
                        if (!compatible) {
                            continue;
                        }
                        mapping[static_cast<std::size_t>(selected_source)] =
                            target;
                        target_used[static_cast<std::size_t>(target)] = 1U;
                        if (self(self, mapped_count + 1U)) {
                            return true;
                        }
                        target_used[static_cast<std::size_t>(target)] = 0U;
                        mapping[static_cast<std::size_t>(selected_source)] = -1;
                        if (aborted) {
                            return false;
                        }
                    }
                    return false;
                };
                if (search(search, 1U)) {
                    generalized_dominators[
                        static_cast<std::size_t>(dominated)]
                        .push_back(dominator);
                    ++dominance_count[static_cast<std::size_t>(dominator)];
                }
                if (deadline_exhausted) {
                    return;
                }
                if (deadline.expired()) {
                    return;
                }
            }
        }
    }

    void build_dff_tables(bool complete_family, int transform_limit) {
        std::vector<DffCandidate> candidates;
        const auto append_candidate = [&](std::int64_t transformed_capacity,
                                          std::vector<std::int64_t> values) {
            DffCandidate candidate;
            candidate.capacity = transformed_capacity;
            candidate.contribution = std::move(values);
            candidate.total = std::accumulate(candidate.contribution.begin(),
                                              candidate.contribution.end(),
                                              std::int64_t{0});
            candidate.root_bound = ceil_div_i64(candidate.total,
                                                candidate.capacity);
            for (const DffCandidate& previous : candidates) {
                if (previous.capacity == candidate.capacity &&
                    previous.contribution == candidate.contribution) {
                    return;
                }
            }
            candidates.push_back(std::move(candidate));
        };

        for (int parameter = 1; parameter <= 20; ++parameter) {
            std::vector<std::int64_t> values(static_cast<std::size_t>(n), 0);
            const std::int64_t scale =
                static_cast<std::int64_t>(parameter) * capacity;
            for (int item = 0; item < n; ++item) {
                const std::int64_t weight =
                    instance.items[static_cast<std::size_t>(item)].weight;
                const std::int64_t numerator =
                    static_cast<std::int64_t>(parameter + 1) * weight;
                values[static_cast<std::size_t>(item)] =
                    numerator % capacity == 0
                        ? static_cast<std::int64_t>(parameter) * weight
                        : (numerator / capacity) * capacity;
            }
            append_candidate(scale, std::move(values));
        }

        std::vector<int> thresholds{capacity / 2};
        for (const Item& item : instance.items) {
            if (item.weight > 0 && 2LL * item.weight <= capacity) {
                thresholds.push_back(item.weight);
            }
        }
        std::sort(thresholds.begin(), thresholds.end());
        thresholds.erase(std::unique(thresholds.begin(), thresholds.end()),
                         thresholds.end());
        for (const int threshold : thresholds) {
            if (threshold <= 0 || 2LL * threshold > capacity) {
                continue;
            }
            std::vector<std::int64_t> values(static_cast<std::size_t>(n), 0);
            for (int item = 0; item < n; ++item) {
                const int weight =
                    instance.items[static_cast<std::size_t>(item)].weight;
                if (weight > capacity - threshold) {
                    values[static_cast<std::size_t>(item)] = capacity;
                } else if (weight >= threshold) {
                    values[static_cast<std::size_t>(item)] = weight;
                }
            }
            append_candidate(capacity, std::move(values));
        }
        const auto candidate_better = [](const DffCandidate& lhs,
                                         const DffCandidate& rhs) {
            if (lhs.root_bound != rhs.root_bound) {
                return lhs.root_bound > rhs.root_bound;
            }
            const long double lhs_scaled =
                static_cast<long double>(lhs.total) * rhs.capacity;
            const long double rhs_scaled =
                static_cast<long double>(rhs.total) * lhs.capacity;
            if (lhs_scaled != rhs_scaled) {
                return lhs_scaled > rhs_scaled;
            }
            return lhs.capacity < rhs.capacity;
        };
        std::stable_sort(candidates.begin(), candidates.end(), candidate_better);
        constexpr std::size_t kMaximumBaseTransforms = 7U;
        if (candidates.size() > kMaximumBaseTransforms) {
            candidates.resize(kMaximumBaseTransforms);
        }

        if (complete_family) {
            std::vector<int> weights;
            weights.reserve(static_cast<std::size_t>(n));
            for (const Item& item : instance.items) {
                weights.push_back(item.weight);
            }
            constexpr std::size_t kMaximumAdditionalTransforms = 8U;
            const DffTransformSet transforms =
                select_ranked_complete_dff_transforms(
                    weights, capacity, true,
                    kMaximumAdditionalTransforms + candidates.size());
            std::vector<DffCandidate> additional;
            additional.reserve(transforms.size());
            for (std::size_t transform = 0; transform < transforms.size();
                 ++transform) {
                DffCandidate candidate;
                candidate.capacity = transforms.capacities[transform];
                const std::int64_t* row = transforms.row(transform);
                candidate.contribution.assign(row, row + n);
                candidate.total = std::accumulate(
                    candidate.contribution.begin(),
                    candidate.contribution.end(), std::int64_t{0});
                candidate.root_bound =
                    ceil_div_i64(candidate.total, candidate.capacity);
                additional.push_back(std::move(candidate));
            }
            std::stable_sort(additional.begin(), additional.end(),
                             candidate_better);
            std::size_t appended = 0;
            for (DffCandidate& candidate : additional) {
                const bool duplicate = std::any_of(
                    candidates.begin(), candidates.end(),
                    [&](const DffCandidate& selected) {
                        return selected.capacity == candidate.capacity &&
                               selected.contribution == candidate.contribution;
                    });
                if (!duplicate) {
                    candidates.push_back(std::move(candidate));
                    if (++appended == kMaximumAdditionalTransforms) {
                        break;
                    }
                }
            }
        }
        for (const DffCandidate& candidate : candidates) {
            root_dff_bound = std::max(root_dff_bound,
                                      candidate.root_bound);
        }
        if (transform_limit > 0 &&
            candidates.size() > static_cast<std::size_t>(transform_limit)) {
            candidates.resize(static_cast<std::size_t>(transform_limit));
        }
        dff_capacity.reserve(candidates.size());
        dff_total.reserve(candidates.size());
        dff_contribution.resize(candidates.size() * static_cast<std::size_t>(n));
        for (const DffCandidate& candidate : candidates) {
            dff_capacity.push_back(candidate.capacity);
            dff_total.push_back(candidate.total);
        }
        for (int item = 0; item < n; ++item) {
            for (std::size_t transform = 0; transform < candidates.size();
                 ++transform) {
                dff_contribution[static_cast<std::size_t>(item) *
                                     candidates.size() +
                                 transform] =
                    candidates[transform]
                        .contribution[static_cast<std::size_t>(item)];
            }
        }
    }
};

struct StateStore {
    static constexpr std::size_t kChunkShift = 14U;
    static constexpr std::size_t kChunkStates = std::size_t{1} << kChunkShift;
    static constexpr std::size_t kChunkMask = kChunkStates - 1U;

private:
    struct Chunk {
        Chunk(std::size_t chunk_capacity,
              std::size_t key_word_count,
              std::size_t transform_count,
              bool compact_dff,
              bool store_assigned_hash,
              bool store_profile_link)
            : capacity(chunk_capacity),
              keys(new std::uint64_t[chunk_capacity * key_word_count]),
              hashes(new std::uint64_t[chunk_capacity]),
              assigned_hashes(store_assigned_hash
                                  ? new std::uint64_t[chunk_capacity]
                                  : nullptr),
              parents(new std::uint32_t[chunk_capacity]),
              depths(new int[chunk_capacity]),
              bounds(new int[chunk_capacity]),
              versions(new std::uint32_t[chunk_capacity]),
              queued(new unsigned char[chunk_capacity]),
              assigned_weights(new std::int64_t[chunk_capacity]),
              assigned_counts(new int[chunk_capacity]),
              compact_dff_sums(transform_count == 0U || !compact_dff
                                   ? nullptr
                                   : new std::uint32_t[
                                         chunk_capacity * transform_count]),
              wide_dff_sums(transform_count == 0U || compact_dff
                                ? nullptr
                                : new std::int64_t[
                                      chunk_capacity * transform_count]),
              profile_next(store_profile_link
                               ? new std::uint32_t[chunk_capacity]
                               : nullptr) {}

        std::size_t capacity = 0U;
        std::unique_ptr<std::uint64_t[]> keys;
        std::unique_ptr<std::uint64_t[]> hashes;
        std::unique_ptr<std::uint64_t[]> assigned_hashes;
        std::unique_ptr<std::uint32_t[]> parents;
        std::unique_ptr<int[]> depths;
        std::unique_ptr<int[]> bounds;
        std::unique_ptr<std::uint32_t[]> versions;
        std::unique_ptr<unsigned char[]> queued;
        std::unique_ptr<std::int64_t[]> assigned_weights;
        std::unique_ptr<int[]> assigned_counts;
        std::unique_ptr<std::uint32_t[]> compact_dff_sums;
        std::unique_ptr<std::int64_t[]> wide_dff_sums;
        std::unique_ptr<std::uint32_t[]> profile_next;
    };

    template <class T, std::unique_ptr<T[]> Chunk::*Member>
    class Field {
    public:
        explicit Field(StateStore* store) noexcept : store_(store) {}

        [[nodiscard]] T& operator[](std::size_t state) noexcept {
            Chunk& chunk = *store_->chunks_[state >> kChunkShift];
            return (chunk.*Member)[state & kChunkMask];
        }

        [[nodiscard]] const T& operator[](std::size_t state) const noexcept {
            const Chunk& chunk = *store_->chunks_[state >> kChunkShift];
            return (chunk.*Member)[state & kChunkMask];
        }

    private:
        StateStore* store_ = nullptr;
    };

public:
    StateStore(std::size_t key_word_count,
               const std::vector<std::int64_t>& dff_totals,
               bool assigned_hash_equals_key_hash,
               bool profile_links_enabled)
        : hashes(this),
          parents(this),
          depths(this),
          bounds(this),
          versions(this),
          queued(this),
          assigned_weights(this),
          assigned_counts(this),
          key_words(key_word_count),
          dff_count(static_cast<int>(dff_totals.size())),
          // Every residual is obtained by subtracting nonnegative item
          // contributions from its root total. Hence the root totals certify
          // whether all per-state residuals fit losslessly in 32 bits.
          compact_dff(std::all_of(
              dff_totals.begin(), dff_totals.end(),
              [](std::int64_t value) {
                  return value >= 0 &&
                      static_cast<std::uint64_t>(value) <=
                          std::numeric_limits<std::uint32_t>::max();
              })),
          store_assigned_hash(!assigned_hash_equals_key_hash),
          // Cooldown-profile chains do not exist for SALBP-I or BPP-P.
          store_profile_link(profile_links_enabled) {
        chunks_.reserve(64U);
        refresh_memory_bytes();
    }

    StateStore(const StateStore&) = delete;
    StateStore& operator=(const StateStore&) = delete;

    [[nodiscard]] bool prepare_append(std::uint64_t other_memory_bytes,
                                      std::uint64_t memory_limit) {
        const std::size_t required = size() + 1U;
        if (required <= capacity_states_) {
            return saturated_add(other_memory_bytes, memory_bytes()) <=
                   memory_limit;
        }
        if (required > kMaximumStateCount) {
            return false;
        }
        const std::uint64_t remaining =
            kMaximumStateCount - capacity_states_;
        const std::size_t chunk_capacity = static_cast<std::size_t>(
            std::min<std::uint64_t>(kChunkStates, remaining));
        if (chunk_capacity == 0U) {
            return false;
        }
        if ((key_words != 0U &&
             chunk_capacity >
                 std::numeric_limits<std::size_t>::max() / key_words) ||
            (dff_count > 0 &&
             chunk_capacity >
                 std::numeric_limits<std::size_t>::max() /
                     static_cast<std::size_t>(dff_count))) {
            return false;
        }
        const std::uint64_t chunk_bytes =
            chunk_memory_bytes(chunk_capacity);
        const std::uint64_t required_memory = saturated_add(
            other_memory_bytes,
            saturated_add(memory_bytes(), chunk_bytes));
        if (required_memory > memory_limit) {
            return false;
        }
        try {
            chunks_.push_back(std::make_unique<Chunk>(
                chunk_capacity, key_words,
                static_cast<std::size_t>(dff_count), compact_dff,
                store_assigned_hash, store_profile_link));
        } catch (const std::bad_alloc&) {
            return false;
        } catch (const std::length_error&) {
            return false;
        }
        capacity_states_ += chunk_capacity;
        refresh_memory_bytes();
        return saturated_add(other_memory_bytes, memory_bytes()) <=
               memory_limit;
    }

    [[nodiscard]] std::size_t size() const noexcept { return size_; }

    [[nodiscard]] const std::uint64_t* key(std::uint32_t state) const noexcept {
        const std::size_t index = static_cast<std::size_t>(state);
        const Chunk& chunk = *chunks_[index >> kChunkShift];
        return chunk.keys.get() + (index & kChunkMask) * key_words;
    }

    [[nodiscard]] std::uint64_t* mutable_key(std::uint32_t state) noexcept {
        const std::size_t index = static_cast<std::size_t>(state);
        Chunk& chunk = *chunks_[index >> kChunkShift];
        return chunk.keys.get() + (index & kChunkMask) * key_words;
    }

    void load_dff(std::uint32_t state,
                  std::int64_t* destination) const noexcept {
        if (dff_count == 0) {
            return;
        }
        const std::size_t index = static_cast<std::size_t>(state);
        const Chunk& chunk = *chunks_[index >> kChunkShift];
        const std::size_t offset =
            (index & kChunkMask) * static_cast<std::size_t>(dff_count);
        if (compact_dff) {
            const std::uint32_t* source =
                chunk.compact_dff_sums.get() + offset;
            std::copy(source, source + dff_count, destination);
        } else {
            const std::int64_t* source = chunk.wide_dff_sums.get() + offset;
            std::copy(source, source + dff_count, destination);
        }
    }

    [[nodiscard]] std::uint64_t assigned_hash(
        std::size_t state) const noexcept {
        const Chunk& chunk = *chunks_[state >> kChunkShift];
        const std::size_t offset = state & kChunkMask;
        return store_assigned_hash ? chunk.assigned_hashes[offset]
                                   : chunk.hashes[offset];
    }

    [[nodiscard]] std::uint32_t profile_link(
        std::size_t state) const noexcept {
        assert(store_profile_link);
        const Chunk& chunk = *chunks_[state >> kChunkShift];
        return chunk.profile_next[state & kChunkMask];
    }

    void set_profile_link(std::size_t state, std::uint32_t next) noexcept {
        assert(store_profile_link);
        Chunk& chunk = *chunks_[state >> kChunkShift];
        chunk.profile_next[state & kChunkMask] = next;
    }

    [[nodiscard]] std::uint32_t append(const std::uint64_t* state_key,
                                       std::uint64_t hash,
                                       std::uint64_t assigned_hash,
                                       std::uint32_t parent,
                                       int depth,
                                       int bound,
                                       std::int64_t assigned_weight,
                                       int assigned_count,
                                       const std::int64_t* transformed_sums) {
        if (size() >= static_cast<std::size_t>(kInvalidState)) {
            throw std::overflow_error("BBR state identifier overflow");
        }
        const std::uint32_t state = static_cast<std::uint32_t>(size());
        const std::size_t index = static_cast<std::size_t>(state);
        Chunk& chunk = *chunks_[index >> kChunkShift];
        const std::size_t offset = index & kChunkMask;
        std::copy(state_key, state_key + key_words,
                  chunk.keys.get() + offset * key_words);
        chunk.hashes[offset] = hash;
        if (store_assigned_hash) {
            chunk.assigned_hashes[offset] = assigned_hash;
        } else {
            assert(hash == assigned_hash);
        }
        chunk.parents[offset] = parent;
        chunk.depths[offset] = depth;
        chunk.bounds[offset] = bound;
        chunk.versions[offset] = 1U;
        chunk.queued[offset] = 1U;
        chunk.assigned_weights[offset] = assigned_weight;
        chunk.assigned_counts[offset] = assigned_count;
        if (dff_count > 0) {
            const std::size_t dff_offset =
                offset * static_cast<std::size_t>(dff_count);
            if (compact_dff) {
                std::uint32_t* destination =
                    chunk.compact_dff_sums.get() + dff_offset;
                for (int transform = 0; transform < dff_count; ++transform) {
                    assert(transformed_sums[transform] >= 0);
                    assert(static_cast<std::uint64_t>(
                               transformed_sums[transform]) <=
                           std::numeric_limits<std::uint32_t>::max());
                    destination[transform] = static_cast<std::uint32_t>(
                        transformed_sums[transform]);
                }
            } else {
                std::copy(
                    transformed_sums, transformed_sums + dff_count,
                    chunk.wide_dff_sums.get() + dff_offset);
            }
        }
        if (store_profile_link) {
            chunk.profile_next[offset] = kInvalidState;
        }
        ++size_;
        return state;
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        return memory_bytes_cached;
    }

private:
    [[nodiscard]] std::uint64_t chunk_memory_bytes(
        std::size_t count) const noexcept {
        std::size_t scalar_bytes =
            sizeof(std::uint64_t) +
            2U * sizeof(std::uint32_t) +
            3U * sizeof(int) + sizeof(unsigned char) +
            sizeof(std::int64_t);
        if (store_assigned_hash) {
            scalar_bytes += sizeof(std::uint64_t);
        }
        if (store_profile_link) {
            scalar_bytes += sizeof(std::uint32_t);
        }
        std::uint64_t per_state = saturated_multiply(
            key_words, sizeof(std::uint64_t));
        per_state = saturated_add(
            per_state,
            saturated_multiply(static_cast<std::size_t>(dff_count),
                               compact_dff ? sizeof(std::uint32_t)
                                           : sizeof(std::int64_t)));
        per_state = saturated_add(per_state, scalar_bytes);
        return saturated_multiply(count, per_state);
    }

    void refresh_memory_bytes() noexcept {
        memory_bytes_cached = vector_memory_bytes(chunks_);
        for (const auto& chunk : chunks_) {
            memory_bytes_cached = saturated_add(
                memory_bytes_cached,
                saturated_add(sizeof(Chunk),
                              chunk_memory_bytes(chunk->capacity)));
        }
    }

public:
    Field<std::uint64_t, &Chunk::hashes> hashes;
    Field<std::uint32_t, &Chunk::parents> parents;
    Field<int, &Chunk::depths> depths;
    Field<int, &Chunk::bounds> bounds;
    Field<std::uint32_t, &Chunk::versions> versions;
    Field<unsigned char, &Chunk::queued> queued;
    Field<std::int64_t, &Chunk::assigned_weights> assigned_weights;
    Field<int, &Chunk::assigned_counts> assigned_counts;
    std::size_t key_words = 0;
    int dff_count = 0;
    bool compact_dff = false;
    bool store_assigned_hash = true;
    bool store_profile_link = true;
    std::size_t size_ = 0U;
    std::size_t capacity_states_ = 0U;
    std::vector<std::unique_ptr<Chunk>> chunks_;
    std::uint64_t memory_bytes_cached = 0;
};

class ExactStateTable {
public:
    ExactStateTable(StateStore& state_store, BbrStatistics& counters,
                    std::uint64_t memory_limit)
        : store_(state_store), statistics_(counters),
          memory_limit_(memory_limit), slots_(1024U, kInvalidState),
          fingerprints_(1024U, 0U) {}

    [[nodiscard]] std::uint32_t find(const std::uint64_t* key,
                                     std::uint64_t hash) {
        ++statistics_.hash_lookups;
        std::size_t slot = static_cast<std::size_t>(hash) & (slots_.size() - 1U);
        const unsigned char fingerprint =
            static_cast<unsigned char>(hash >> 56U);
        while (true) {
            ++statistics_.hash_probes;
            const std::uint32_t state = slots_[slot];
            if (state == kInvalidState) {
                return kInvalidState;
            }
            if (fingerprints_[slot] == fingerprint &&
                store_.hashes[static_cast<std::size_t>(state)] == hash &&
                masks_equal(store_.key(state), key, store_.key_words)) {
                return state;
            }
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
    }

    [[nodiscard]] bool prepare_insert(std::uint64_t other_memory_bytes) {
        if ((store_.size() + 1U) * 10U < slots_.size() * 7U) {
            return true;
        }
        const std::size_t new_size = slots_.size() * 2U;
        const std::uint64_t transient =
            store_.memory_bytes() + memory_bytes() + other_memory_bytes +
            static_cast<std::uint64_t>(new_size) *
                (sizeof(std::uint32_t) + sizeof(unsigned char));
        if (transient > memory_limit_) {
            return false;
        }
        std::vector<std::uint32_t> replacement(new_size, kInvalidState);
        std::vector<unsigned char> replacement_fingerprints(new_size, 0U);
        for (std::uint32_t state = 0;
             state < static_cast<std::uint32_t>(store_.size()); ++state) {
            std::size_t slot =
                static_cast<std::size_t>(store_.hashes[state]) & (new_size - 1U);
            while (replacement[slot] != kInvalidState) {
                slot = (slot + 1U) & (new_size - 1U);
            }
            replacement[slot] = state;
            replacement_fingerprints[slot] = static_cast<unsigned char>(
                store_.hashes[state] >> 56U);
        }
        slots_.swap(replacement);
        fingerprints_.swap(replacement_fingerprints);
        return true;
    }

    void insert(std::uint32_t state) {
        const std::uint64_t hash = store_.hashes[static_cast<std::size_t>(state)];
        std::size_t slot = static_cast<std::size_t>(hash) & (slots_.size() - 1U);
        while (slots_[slot] != kInvalidState) {
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
        slots_[slot] = state;
        fingerprints_[slot] = static_cast<unsigned char>(hash >> 56U);
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        return static_cast<std::uint64_t>(slots_.capacity()) *
                   sizeof(std::uint32_t) +
               static_cast<std::uint64_t>(fingerprints_.capacity()) *
                   sizeof(unsigned char);
    }

    void set_memory_limit(std::uint64_t memory_limit) noexcept {
        memory_limit_ = memory_limit;
    }

private:
    StateStore& store_;
    BbrStatistics& statistics_;
    std::uint64_t memory_limit_ = 0;
    std::vector<std::uint32_t> slots_;
    std::vector<unsigned char> fingerprints_;
};

class AssignedProfileTable {
public:
    AssignedProfileTable(StateStore& state_store,
                         std::size_t assigned_blocks,
                         int cooldown_levels,
                         std::uint64_t memory_limit)
        : store_(state_store),
          assigned_blocks_(assigned_blocks),
          cooldown_levels_(cooldown_levels),
          memory_limit_(memory_limit),
          slots_(1024U, kInvalidState) {}

    [[nodiscard]] bool dominated(const std::uint64_t* key,
                                 std::uint64_t assigned_hash,
                                 int depth) const noexcept {
        const std::uint32_t head = find_head(key, assigned_hash);
        if (head == kInvalidState) {
            return false;
        }
        const std::size_t cooldown_words =
            assigned_blocks_ * static_cast<std::size_t>(cooldown_levels_);
        for (std::uint32_t state = head; state != kInvalidState;
             state = store_.profile_link(static_cast<std::size_t>(state))) {
            if (store_.depths[static_cast<std::size_t>(state)] <= depth &&
                mask_subset(store_.key(state) + assigned_blocks_,
                            key + assigned_blocks_, cooldown_words)) {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] bool prepare_insert(const std::uint64_t* key,
                                      std::uint64_t assigned_hash,
                                      std::uint64_t other_memory_bytes) {
        if (find_head(key, assigned_hash) != kInvalidState) {
            return true;
        }
        if ((group_count_ + 1U) * 10U < slots_.size() * 7U) {
            return true;
        }
        const std::size_t new_size = slots_.size() * 2U;
        const std::uint64_t transient =
            store_.memory_bytes() + memory_bytes() + other_memory_bytes +
            static_cast<std::uint64_t>(new_size) * sizeof(std::uint32_t);
        if (transient > memory_limit_) {
            return false;
        }
        std::vector<std::uint32_t> replacement(new_size, kInvalidState);
        for (const std::uint32_t head : slots_) {
            if (head == kInvalidState) {
                continue;
            }
            std::size_t slot = static_cast<std::size_t>(
                                   store_.assigned_hash(
                                       static_cast<std::size_t>(head))) &
                               (new_size - 1U);
            while (replacement[slot] != kInvalidState) {
                slot = (slot + 1U) & (new_size - 1U);
            }
            replacement[slot] = head;
        }
        slots_.swap(replacement);
        return true;
    }

    void insert(std::uint32_t state) noexcept {
        const std::uint64_t assigned_hash =
            store_.assigned_hash(static_cast<std::size_t>(state));
        std::size_t slot = static_cast<std::size_t>(assigned_hash) &
                           (slots_.size() - 1U);
        while (true) {
            const std::uint32_t head = slots_[slot];
            if (head == kInvalidState) {
                slots_[slot] = state;
                ++group_count_;
                return;
            }
            if (store_.assigned_hash(static_cast<std::size_t>(head)) ==
                    assigned_hash &&
                masks_equal(store_.key(head), store_.key(state),
                            assigned_blocks_)) {
                store_.set_profile_link(static_cast<std::size_t>(state), head);
                slots_[slot] = state;
                return;
            }
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
    }

    [[nodiscard]] std::uint64_t memory_bytes() const noexcept {
        return static_cast<std::uint64_t>(slots_.capacity()) *
               sizeof(std::uint32_t);
    }

    void set_memory_limit(std::uint64_t memory_limit) noexcept {
        memory_limit_ = memory_limit;
    }

private:
    [[nodiscard]] std::uint32_t find_head(
        const std::uint64_t* key,
        std::uint64_t assigned_hash) const noexcept {
        std::size_t slot = static_cast<std::size_t>(assigned_hash) &
                           (slots_.size() - 1U);
        while (true) {
            const std::uint32_t head = slots_[slot];
            if (head == kInvalidState) {
                return kInvalidState;
            }
            if (store_.assigned_hash(static_cast<std::size_t>(head)) ==
                    assigned_hash &&
                masks_equal(store_.key(head), key, assigned_blocks_)) {
                return head;
            }
            slot = (slot + 1U) & (slots_.size() - 1U);
        }
    }

    StateStore& store_;
    std::size_t assigned_blocks_ = 0;
    int cooldown_levels_ = 0;
    std::uint64_t memory_limit_ = 0;
    std::size_t group_count_ = 0;
    std::vector<std::uint32_t> slots_;
};

struct QueueEntry {
    std::int64_t best_measure = 0;
    int longest_path_tie = 0;
    int lower_bound = 0;
    std::uint32_t state = kInvalidState;
    std::uint32_t version = 0;
};

[[nodiscard]] bool queue_entry_worse(const QueueEntry& lhs,
                                     const QueueEntry& rhs) noexcept {
    return std::tie(lhs.best_measure, lhs.longest_path_tie, lhs.lower_bound,
                    lhs.state) >
           std::tie(rhs.best_measure, rhs.longest_path_tie, rhs.lower_bound,
                    rhs.state);
}

class BbrEngine {
public:
    BbrEngine(const PreparedInstance& prepared,
              int initial_lower_bound,
              const Assignment& initial_incumbent,
              const Config& config,
              Deadline& deadline)
        : prepared_(prepared),
          instance_(prepared.search_instance),
          precomputed_(instance_, config.bbr_enable_jackson,
                       config.bbr_enable_generalized_item_dominance,
                       config.bbr_enable_complete_dff,
                       config.bbr_dff_transform_limit, deadline),
          initial_lower_bound_(initial_lower_bound),
          config_(config),
          deadline_(deadline),
          maximum_memory_limit_bytes_(
              checked_memory_limit(config.bbr_memory_limit_mb)),
          store_(precomputed_.key_words, precomputed_.dff_total,
                 precomputed_.salbp_semantics,
                 precomputed_.cooldown_levels > 0),
          memory_limit_bytes_(maximum_memory_limit_bytes_),
          exact_table_(store_, statistics_, memory_limit_bytes_),
          profile_table_(store_, precomputed_.blocks,
                         precomputed_.cooldown_levels, memory_limit_bytes_),
          incumbent_(initial_incumbent),
          upper_bound_(initial_incumbent.bin_count),
          current_dff_sums_(precomputed_.dff_capacity.size(), 0),
          current_key_(precomputed_.key_words, 0U),
          child_key_(precomputed_.key_words, 0U),
          probe_key_(precomputed_.key_words, 0U),
          load_mask_(precomputed_.blocks, 0U),
          replacement_union_(precomputed_.blocks, 0U),
          item_status_(static_cast<std::size_t>(precomputed_.n), 3U),
          remaining_zero_predecessors_(static_cast<std::size_t>(precomputed_.n), 0),
          earliest_offset_(static_cast<std::size_t>(precomputed_.n), 0),
          tail_offset_(static_cast<std::size_t>(precomputed_.n), 0),
          machine_order_(static_cast<std::size_t>(precomputed_.n), 0),
          machine_bucket_weights_(
              static_cast<std::size_t>(
                  std::max(precomputed_.n, initial_incumbent.bin_count) + 1),
              0),
          child_dff_sums_(precomputed_.dff_capacity.size(), 0),
          load_dff_sums_(precomputed_.dff_capacity.size(), 0) {
        initial_lower_bound_ = std::max(initial_lower_bound_,
                                        precomputed_.root_dff_bound);
        if (initial_lower_bound_ <= 0 || upper_bound_ <= 0 ||
            initial_lower_bound_ > upper_bound_) {
            throw std::invalid_argument("invalid initial bounds for BBR");
        }
        std::string diagnostic;
        if (!check_assignment(instance_, incumbent_, &diagnostic)) {
            throw std::logic_error(
                "prepared BBR incumbent is invalid: " + diagnostic);
        }
        statistics_.attempted = true;
        statistics_.exact_phase_attempted = true;
        statistics_.item_dominance_enabled = config_.bbr_enable_jackson;
        statistics_.generalized_item_dominance_enabled =
            config_.bbr_enable_generalized_item_dominance;
        statistics_.paper_queue_order_enabled =
            config_.bbr_enable_paper_queue_order;
        statistics_.complete_dff_enabled = config_.bbr_enable_complete_dff;
        statistics_.binlb_enabled = config_.bbr_enable_binlb;
        statistics_.conflict_binlb_enabled =
            config_.bbr_enable_binlb && config_.bbr_enable_conflict_binlb;
        statistics_.binlb_call_time_limit_seconds =
            config_.bbr_binlb_call_time_limit_seconds;
        statistics_.binlb_total_time_limit_seconds =
            config_.bbr_binlb_total_time_limit_seconds;
        statistics_.binlb_node_limit = config_.bbr_binlb_node_limit;
        statistics_.binlb_load_limit = config_.bbr_binlb_load_limit;
        statistics_.binlb_memo_limit = config_.bbr_binlb_memo_limit;
        statistics_.binlb_max_items = config_.bbr_binlb_max_items;
        statistics_.conflict_binlb_call_time_limit_seconds =
            config_.bbr_conflict_binlb_call_time_limit_seconds;
        statistics_.conflict_binlb_node_limit =
            config_.bbr_conflict_binlb_node_limit;
        statistics_.dff_transform_count = precomputed_.dff_capacity.size();
        statistics_.structured_preprocessing_enabled =
            prepared_.structured_preprocessing_enabled;
        statistics_.structured_search_items =
            static_cast<std::uint64_t>(precomputed_.n);
        statistics_.structured_items_removed =
            prepared_.original_item_count > precomputed_.n
                ? static_cast<std::uint64_t>(prepared_.original_item_count -
                                             precomputed_.n)
                : 0U;
        statistics_.structured_full_capacity_items =
            prepared_.full_capacity_removals.size();
        statistics_.structured_prefix_bins =
            prepared_.fixed_prefix_bins.size();
        statistics_.structured_suffix_bins =
            prepared_.fixed_suffix_bins.size();
        statistics_.structured_fixed_bins =
            static_cast<std::uint64_t>(prepared_.fixed_bin_offset);
        statistics_.generalized_item_dominance_search_nodes =
            precomputed_.generalized_dominance_search_nodes;
        statistics_.generalized_item_dominance_seconds =
            precomputed_.generalized_dominance_seconds;
        for (const auto& values : precomputed_.generalized_dominators) {
            statistics_.generalized_item_dominance_pairs += values.size();
        }
        statistics_.reverse_direction = prepared_.reversed;
        statistics_.memory_limit_bytes = maximum_memory_limit_bytes_;
        statistics_.seed = config_.seed;
        statistics_.time_limit_seconds = config_.time_limit_seconds;
        statistics_.initialization_time_limit_seconds =
            config_.initialization_time_limit_seconds > 0.0
                ? config_.initialization_time_limit_seconds
                : std::min(3.0, 0.1 * config_.time_limit_seconds);
        bbr12_item_order_ = precomputed_.branch_order;
        bbr12_rank_.resize(static_cast<std::size_t>(precomputed_.n));
        for (int rank = 0; rank < precomputed_.n; ++rank) {
            bbr12_rank_[static_cast<std::size_t>(bbr12_item_order_[
                static_cast<std::size_t>(rank)])] = rank;
        }
        bbr12_successor_offsets_.resize(
            static_cast<std::size_t>(precomputed_.n + 1), 0);
        for (int item = 0; item < precomputed_.n; ++item) {
            bbr12_successor_offsets_[static_cast<std::size_t>(item + 1)] =
                bbr12_successor_offsets_[static_cast<std::size_t>(item)] +
                precomputed_.zero_successor_offset[
                    static_cast<std::size_t>(item + 1)] -
                precomputed_.zero_successor_offset[
                    static_cast<std::size_t>(item)];
        }
        bbr12_successors_.resize(static_cast<std::size_t>(
            bbr12_successor_offsets_.back()));
        for (int item = 0; item < precomputed_.n; ++item) {
            const int source_begin = precomputed_.zero_successor_offset[
                static_cast<std::size_t>(item)];
            const int source_end = precomputed_.zero_successor_offset[
                static_cast<std::size_t>(item + 1)];
            const int target_begin = bbr12_successor_offsets_[
                static_cast<std::size_t>(item)];
            std::copy(
                precomputed_.zero_successors.begin() + source_begin,
                precomputed_.zero_successors.begin() + source_end,
                bbr12_successors_.begin() + target_begin);
            std::stable_sort(
                bbr12_successors_.begin() + target_begin,
                bbr12_successors_.begin() +
                    bbr12_successor_offsets_[
                        static_cast<std::size_t>(item + 1)],
                [&](int lhs, int rhs) {
                    return bbr12_rank_[static_cast<std::size_t>(lhs)] <
                           bbr12_rank_[static_cast<std::size_t>(rhs)];
                });
        }
        bbr12_eligible_.resize(static_cast<std::size_t>(precomputed_.n));
        constexpr int kMaximumFitMaskCapacity = 65'536;
        fit_item_masks_available_ =
            precomputed_.capacity <= kMaximumFitMaskCapacity;
        if (fit_item_masks_available_) {
            fit_item_masks_.assign(
                static_cast<std::size_t>(precomputed_.capacity + 1) *
                    precomputed_.blocks,
                0U);
            for (int item = 0; item < precomputed_.n; ++item) {
                const int weight =
                    instance_.items[static_cast<std::size_t>(item)].weight;
                std::uint64_t* row = fit_item_masks_.data() +
                    static_cast<std::size_t>(weight) * precomputed_.blocks;
                row[static_cast<std::size_t>(item) >> 6U] |=
                    std::uint64_t{1} <<
                    (static_cast<unsigned>(item) & 63U);
            }
            for (int capacity = 1; capacity <= precomputed_.capacity;
                 ++capacity) {
                std::uint64_t* row = fit_item_masks_.data() +
                    static_cast<std::size_t>(capacity) * precomputed_.blocks;
                const std::uint64_t* previous = row - precomputed_.blocks;
                for (std::size_t block = 0; block < precomputed_.blocks;
                     ++block) {
                    row[block] |= previous[block];
                }
            }
        }
        fast_bppp_dominance_masks_enabled_ =
            fit_item_masks_available_ && precomputed_.bppp_semantics;
        if (fast_bppp_dominance_masks_enabled_) {
            boundary_eligible_item_mask_.assign(precomputed_.blocks, 0U);
        }
        if (statistics_.binlb_enabled) {
            BinPackingBoundLimits limits;
            limits.call_time_limit_seconds =
                config_.bbr_binlb_call_time_limit_seconds;
            limits.search_node_limit = config_.bbr_binlb_node_limit;
            limits.nondominated_load_limit_per_state =
                config_.bbr_binlb_load_limit;
            const std::uint64_t estimated_memo_entry_bytes =
                static_cast<std::uint64_t>(precomputed_.n) *
                    sizeof(std::uint16_t) +
                32U;
            const std::uint64_t auxiliary_memo_budget =
                maximum_memory_limit_bytes_ / 8U;
            limits.memo_entry_limit = std::min(
                config_.bbr_binlb_memo_limit,
                auxiliary_memo_budget /
                    std::max<std::uint64_t>(1U,
                                            estimated_memo_entry_bytes));
            limits.maximum_item_count = config_.bbr_binlb_max_items;
            limits.enable_conflicts = config_.bbr_enable_conflict_binlb;
            limits.conflict_call_time_limit_seconds =
                config_.bbr_conflict_binlb_call_time_limit_seconds;
            limits.conflict_search_node_limit =
                config_.bbr_conflict_binlb_node_limit;
            bin_packing_ = std::make_unique<BinPackingBound>(
                instance_, limits);
            statistics_.binlb_conflict_edges =
                bin_packing_->conflict_edge_count();
            conflict_remaining_.assign(precomputed_.blocks, 0U);
        }
        heaps_.resize(static_cast<std::size_t>(upper_bound_ + 1));
        refresh_fixed_memory_bytes();
        heap_memory_bytes_cached_ = vector_memory_bytes(heaps_);
        if (!ensure_memory_capacity(current_memory_bytes())) {
            throw std::bad_alloc();
        }
        update_peak_memory();
    }

    [[nodiscard]] BbrResult solve() {
        const auto search_start = Clock::now();
        try {
            if (deadline_.expired()) {
                timed_out_ = true;
                stop_ = true;
            } else {
                initialize_root();
                search();
            }
        } catch (const std::bad_alloc&) {
            memory_limited_ = true;
            stop_ = true;
        }
        return finish_result(search_start);
    }

private:
    [[nodiscard]] BbrResult finish_result(Clock::time_point search_start) {
        BbrResult result;
        result.attempted = true;
        result.optimal = optimal_ || upper_bound_ <= initial_lower_bound_;
        result.timed_out = timed_out_;
        result.memory_limited = memory_limited_;
        result.incumbent = incumbent_;
        result.certified_lower_bound = result.optimal
            ? upper_bound_
            : certified_open_lower_bound();
        result.certified_lower_bound = std::min(result.certified_lower_bound,
                                                upper_bound_);

        statistics_.timed_out = result.timed_out;
        statistics_.memory_limited = result.memory_limited;
        statistics_.peak_memory_bytes =
            std::max(statistics_.peak_memory_bytes, current_memory_bytes());
        statistics_.search_seconds =
            std::chrono::duration<double>(Clock::now() - search_start).count();
        statistics_.exact_search_seconds = statistics_.search_seconds;
        result.statistics = statistics_;
        return result;
    }

    [[nodiscard]] bool acquire_state_slot() const {
        if (store_.size() >= kMaximumStateCount) {
            throw std::overflow_error(
                "BBR state identifier range exhausted");
        }
        return true;
    }

    [[nodiscard]] static std::uint64_t checked_memory_limit(
        std::uint64_t megabytes) {
        constexpr std::uint64_t kMegabyte = 1024U * 1024U;
        if (megabytes == 0U ||
            megabytes > std::numeric_limits<std::uint64_t>::max() / kMegabyte) {
            throw std::invalid_argument("invalid BBR memory limit");
        }
        return megabytes * kMegabyte;
    }

    void initialize_root() {
        std::fill(child_key_.begin(), child_key_.end(), 0U);
        const int root_bound = compute_lower_bound(
            child_key_.data(), 0, 0, 0, precomputed_.dff_total.data());
        const std::uint64_t root_assigned_hash =
            precomputed_.assigned_set_hash(child_key_.data());
        const std::uint64_t root_hash = precomputed_.salbp_semantics
            ? root_assigned_hash
            : hash_words(child_key_.data(), precomputed_.key_words);
        if (!acquire_state_slot()) {
            return;
        }
        if (!prepare_store_append(
                exact_table_.memory_bytes() + profile_table_.memory_bytes() +
                    heap_memory_bytes() + fixed_memory_bytes()) ||
            !prepare_heap_push(0)) {
            memory_limited_ = true;
            stop_ = true;
            return;
        }
        const std::uint32_t root = store_.append(
            child_key_.data(), root_hash, root_assigned_hash, kInvalidState,
            0, root_bound, 0, 0, precomputed_.dff_total.data());
        exact_table_.insert(root);
        if (profile_dominance_enabled()) {
            profile_table_.insert(root);
        }
        ++statistics_.states_created;
        open_states_ = 1U;
        statistics_.peak_open_states = 1U;
        push_state_unchecked(root);
        update_peak_memory();
    }

    void search() {
        while (!stop_) {
            if (upper_bound_ <= initial_lower_bound_) {
                optimal_ = true;
                break;
            }
            if (deadline_.expired()) {
                timed_out_ = true;
                stop_ = true;
                break;
            }
            const std::uint32_t state = pop_next_state();
            if (state == kInvalidState) {
                optimal_ = true;
                break;
            }
            expand_state(state);
        }
    }

    [[nodiscard]] std::int64_t best_measure(
        std::uint32_t state) const noexcept {
        const int depth = store_.depths[static_cast<std::size_t>(state)];
        const std::int64_t idle =
            static_cast<std::int64_t>(depth) * precomputed_.capacity -
            store_.assigned_weights[static_cast<std::size_t>(state)];
        const int unassigned = precomputed_.n -
            store_.assigned_counts[static_cast<std::size_t>(state)];
        return 50 * idle - static_cast<std::int64_t>(depth) * unassigned;
    }

    [[nodiscard]] bool prepare_heap_push(int depth) {
        if (depth < 0 || static_cast<std::size_t>(depth) >= heaps_.size()) {
            return false;
        }
        auto& heap = heaps_[static_cast<std::size_t>(depth)];
        if (heap.size() < heap.capacity()) {
            return ensure_memory_capacity(current_memory_bytes());
        }
        std::size_t desired = heap.capacity() == 0U
                                  ? 256U
                                  : heap.capacity() * 2U;
        if (desired <= heap.capacity()) {
            return false;
        }
        std::uint64_t new_buffer =
            saturated_multiply(desired, sizeof(QueueEntry));
        if (!ensure_memory_capacity(
                saturated_add(current_memory_bytes(), new_buffer))) {
            desired = heap.size() + 1U;
            new_buffer = saturated_multiply(desired, sizeof(QueueEntry));
            if (!ensure_memory_capacity(
                    saturated_add(current_memory_bytes(), new_buffer))) {
                return false;
            }
        }
        try {
            const std::uint64_t old_buffer =
                vector_memory_bytes(heap);
            heap.reserve(desired);
            heap_memory_bytes_cached_ = saturated_add(
                heap_memory_bytes_cached_ - old_buffer,
                vector_memory_bytes(heap));
        } catch (const std::bad_alloc&) {
            return false;
        } catch (const std::length_error&) {
            return false;
        }
        if (!ensure_memory_capacity(current_memory_bytes())) {
            return false;
        }
        update_peak_memory();
        return true;
    }

    void push_state_unchecked(std::uint32_t state) {
        const int depth = store_.depths[static_cast<std::size_t>(state)];
        if (depth < 0) {
            throw std::logic_error("negative BBR state depth");
        }
        auto& heap = heaps_[static_cast<std::size_t>(depth)];
        const bool paper_order = statistics_.paper_queue_order_enabled;
        heap.push_back(QueueEntry{
            paper_order ? next_queue_machine_numerator_ : best_measure(state),
            paper_order ? next_queue_longest_path_tie_ : 0,
            store_.bounds[static_cast<std::size_t>(state)], state,
            store_.versions[static_cast<std::size_t>(state)]});
        std::push_heap(heap.begin(), heap.end(), queue_entry_worse);
    }

    [[nodiscard]] std::uint32_t pop_next_state() {
        if (heaps_.empty()) {
            return kInvalidState;
        }
        for (std::size_t checked = 0; checked < heaps_.size(); ++checked) {
            const std::size_t depth = next_depth_;
            next_depth_ = (next_depth_ + 1U) % heaps_.size();
            const std::uint32_t state = pop_valid_state_at_depth(depth);
            if (state != kInvalidState) {
                return state;
            }
        }
        return kInvalidState;
    }

    [[nodiscard]] std::uint32_t pop_valid_state_at_depth(
        std::size_t depth) {
        auto& heap = heaps_[depth];
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end(), queue_entry_worse);
            const QueueEntry entry = heap.back();
            heap.pop_back();
            const std::size_t state = entry.state;
            if (entry.state == kInvalidState || state >= store_.size() ||
                store_.queued[state] == 0U ||
                store_.versions[state] != entry.version ||
                store_.depths[state] != static_cast<int>(depth)) {
                continue;
            }
            if (store_.bounds[state] >= upper_bound_) {
                store_.queued[state] = 0U;
                --open_states_;
                ++statistics_.bound_prunes;
                continue;
            }
            return entry.state;
        }
        return kInvalidState;
    }

    void expand_state(std::uint32_t state) {
        current_state_ = state;
        current_version_ = store_.versions[static_cast<std::size_t>(state)];
        current_depth_ = store_.depths[static_cast<std::size_t>(state)];
        current_bound_ = store_.bounds[static_cast<std::size_t>(state)];
        current_assigned_weight_ =
            store_.assigned_weights[static_cast<std::size_t>(state)];
        current_assigned_count_ =
            store_.assigned_counts[static_cast<std::size_t>(state)];
        current_assigned_hash_ =
            store_.assigned_hash(static_cast<std::size_t>(state));
        std::copy(store_.key(state),
                  store_.key(state) + precomputed_.key_words,
                  current_key_.begin());
        if (!current_dff_sums_.empty()) {
            store_.load_dff(state, current_dff_sums_.data());
        }

        ++statistics_.states_expanded;
        prepare_load_enumeration();
        enumerate_loads();
        if (!stop_ &&
            store_.versions[static_cast<std::size_t>(state)] == current_version_ &&
            store_.queued[static_cast<std::size_t>(state)] != 0U) {
            store_.queued[static_cast<std::size_t>(state)] = 0U;
            --open_states_;
        }
    }

    [[nodiscard]] int cooldown_release(const std::uint64_t* key,
                                       int item) const noexcept {
        for (int level = precomputed_.cooldown_levels; level >= 1; --level) {
            const std::uint64_t* row = key +
                static_cast<std::size_t>(level) * precomputed_.blocks;
            if (bit_is_set(row, item)) {
                return level;
            }
        }
        return 0;
    }

    [[nodiscard]] int compute_dff_lower_bound(
        const std::uint64_t* key,
        std::int64_t numerator,
        std::size_t transform) const {
        const std::int64_t denominator =
            precomputed_.dff_capacity[transform];
        if (numerator < 0 || denominator <= 0) {
            std::int64_t expected_child =
                precomputed_.dff_total[transform];
            precomputed_.for_each_set_bit(key, [&](int item) {
                expected_child -=
                    precomputed_.transformed_weights_for_item(item)[transform];
            });
            std::int64_t expected_current =
                precomputed_.dff_total[transform];
            std::int64_t expected_load = 0;
            precomputed_.for_each_set_bit(
                current_key_.data(), [&](int item) {
                    expected_current -=
                        precomputed_.transformed_weights_for_item(
                            item)[transform];
                });
            precomputed_.for_each_set_bit(
                load_mask_.data(), [&](int item) {
                    expected_load +=
                        precomputed_.transformed_weights_for_item(
                            item)[transform];
                });
            throw std::logic_error(
                "invalid BBR DFF residual"
                " (transform=" + std::to_string(transform) +
                ", capacity=" + std::to_string(denominator) +
                ", actual=" + std::to_string(numerator) +
                ", expected_child=" + std::to_string(expected_child) +
                ", current=" +
                std::to_string(current_dff_sums_[transform]) +
                ", expected_current=" +
                std::to_string(expected_current) +
                ", load=" + std::to_string(load_dff_sums_[transform]) +
                ", expected_load=" + std::to_string(expected_load) +
                ", depth=" + std::to_string(current_depth_ + 1) +
                ", assigned_count=" +
                std::to_string(current_assigned_count_ +
                               current_load_count_) + ")");
        }
        return static_cast<int>((numerator + denominator - 1) / denominator);
    }

    [[nodiscard]] int compute_cheap_lower_bound(
        const std::uint64_t* key,
        int depth,
        std::int64_t assigned_weight,
        const std::int64_t* dff_sums) {
        const std::int64_t remaining_weight =
            instance_.total_weight - assigned_weight;
        next_queue_machine_numerator_ = 0;
        next_queue_longest_path_tie_ = 0;
        int bin_packing_bound =
            ceil_div_i64(remaining_weight, precomputed_.capacity);
        for (std::size_t transform = 0;
             transform < precomputed_.dff_capacity.size(); ++transform) {
            bin_packing_bound = std::max(
                bin_packing_bound,
                compute_dff_lower_bound(
                    key, dff_sums[transform], transform));
        }
        int lower_bound = std::max(initial_lower_bound_,
                                   depth + bin_packing_bound);
        if (lower_bound >= upper_bound_ || remaining_weight == 0) {
            return lower_bound;
        }

        return lower_bound;
    }

    [[nodiscard]] int compute_precedence_lower_bound(
        const std::uint64_t* key,
        int depth,
        int lower_bound) {
        const std::uint64_t* assigned = key;
        if (precomputed_.salbp_semantics) {
            return lower_bound;
        }

        remaining_machine_count_ = 0;
        remaining_maximum_tail_ = 0;
        remaining_maximum_earliest_ = 0;
        remaining_minimum_tail_ = precomputed_.n + 1;
        for (const int item : instance_.topological_order) {
            if (bit_is_set(assigned, item)) {
                continue;
            }
            int earliest = cooldown_release(key, item);
            for (const auto& [predecessor, separation] :
                 instance_.predecessor_arcs[static_cast<std::size_t>(item)]) {
                if (!bit_is_set(assigned, predecessor)) {
                    earliest = std::max(
                        earliest,
                        earliest_offset_[static_cast<std::size_t>(predecessor)] +
                            separation);
                }
            }
            earliest_offset_[static_cast<std::size_t>(item)] = earliest;
        }
        for (auto order = instance_.topological_order.rbegin();
             order != instance_.topological_order.rend(); ++order) {
            const int item = *order;
            if (bit_is_set(assigned, item)) {
                continue;
            }
            int tail = 0;
            for (const auto& [successor, separation] :
                 instance_.successor_arcs[static_cast<std::size_t>(item)]) {
                if (!bit_is_set(assigned, successor)) {
                    tail = std::max(
                        tail,
                        separation +
                            tail_offset_[static_cast<std::size_t>(successor)]);
                }
            }
            tail_offset_[static_cast<std::size_t>(item)] = tail;
            machine_order_[static_cast<std::size_t>(remaining_machine_count_++)] =
                item;
            remaining_maximum_tail_ = std::max(remaining_maximum_tail_, tail);
            remaining_maximum_earliest_ = std::max(
                remaining_maximum_earliest_,
                earliest_offset_[static_cast<std::size_t>(item)]);
            if (statistics_.paper_queue_order_enabled) {
                remaining_minimum_tail_ = std::min(
                    remaining_minimum_tail_, tail);
            }
            lower_bound = std::max(
                lower_bound,
                depth + earliest_offset_[static_cast<std::size_t>(item)] +
                    tail + 1);
            if (lower_bound >= upper_bound_) {
                return lower_bound;
            }
        }
        return lower_bound;
    }

    [[nodiscard]] int strengthen_post_memory_lower_bound(
        const std::uint64_t* key,
        int depth,
        std::int64_t assigned_weight,
        int assigned_count,
        int lower_bound) {
        const std::uint64_t* assigned = key;
        const std::int64_t remaining_weight =
            instance_.total_weight - assigned_weight;

        if (!precomputed_.salbp_semantics) {
            ++statistics_.machine_bound_calls;
            const auto apply_machine_bound = [&](const std::vector<int>& offsets,
                                                 int maximum_offset,
                                                 std::int64_t* raw_numerator) {
            if (raw_numerator == nullptr && depth + maximum_offset +
                    ceil_div_i64(remaining_weight, precomputed_.capacity) <=
                lower_bound) {
                return lower_bound;
            }

            if (maximum_offset < 0 ||
                static_cast<std::size_t>(maximum_offset) >=
                    machine_bucket_weights_.size()) {
                throw std::logic_error(
                    "precedence machine offset out of incumbent range");
            }
            std::fill(machine_bucket_weights_.begin(),
                      machine_bucket_weights_.begin() + maximum_offset + 1,
                      0);
            for (int position = 0; position < remaining_machine_count_;
                 ++position) {
                const int item =
                    machine_order_[static_cast<std::size_t>(position)];
                const int offset = offsets[static_cast<std::size_t>(item)];
                if (offset < 0 ||
                    static_cast<std::size_t>(offset) >=
                        machine_bucket_weights_.size()) {
                    throw std::logic_error(
                        "precedence machine offset out of incumbent range");
                }
                machine_bucket_weights_[static_cast<std::size_t>(offset)] +=
                    instance_.items[static_cast<std::size_t>(item)].weight;
            }
            std::int64_t prefix_weight = 0;
            std::int64_t numerator = 0;
            int bound = lower_bound;
            for (int offset = maximum_offset; offset >= 0; --offset) {
                const std::int64_t bucket_weight =
                    machine_bucket_weights_[static_cast<std::size_t>(offset)];
                if (bucket_weight == 0) {
                    continue;
                }
                prefix_weight += bucket_weight;
                numerator = std::max(
                    numerator,
                    prefix_weight + static_cast<std::int64_t>(offset) *
                                        precomputed_.capacity);
                bound = std::max(
                    bound,
                    depth + offset +
                        ceil_div_i64(prefix_weight, precomputed_.capacity));
            }
            if (raw_numerator != nullptr) {
                *raw_numerator = numerator;
            }
            return bound;
            };
            std::int64_t queue_machine_numerator = 0;
            lower_bound = apply_machine_bound(
                tail_offset_, remaining_maximum_tail_,
                statistics_.paper_queue_order_enabled
                    ? &queue_machine_numerator
                    : nullptr);
            if (statistics_.paper_queue_order_enabled) {
                next_queue_machine_numerator_ = queue_machine_numerator;
                next_queue_longest_path_tie_ =
                    remaining_minimum_tail_ <= precomputed_.n
                        ? remaining_minimum_tail_
                        : 0;
            }
            if (lower_bound < upper_bound_) {
                lower_bound = apply_machine_bound(
                    earliest_offset_, remaining_maximum_earliest_, nullptr);
            }
            if (lower_bound >= upper_bound_) {
                ++statistics_.machine_bound_prunes;
                return lower_bound;
            }
        } else if (statistics_.paper_queue_order_enabled) {
            next_queue_machine_numerator_ = remaining_weight;
            next_queue_longest_path_tie_ = 0;
        }

        if (config_.bbr_enable_closure_bound) {
            ++statistics_.closure_bound_calls;
            const int bound_before_closure = lower_bound;
            for (const int item : precomputed_.closure_bound_order) {
                if (static_cast<std::int64_t>(depth) +
                        precomputed_.closure_bound_static_potential[
                            static_cast<std::size_t>(item)] <=
                    lower_bound) {
                    break;
                }
                if (bit_is_set(assigned, item)) {
                    continue;
                }
                const int prefix_bin_upper_bound = std::max(
                    precomputed_.predecessor_closure_capacity_bins[
                        static_cast<std::size_t>(item)],
                    earliest_offset_[static_cast<std::size_t>(item)] + 1);
                const int suffix_bin_upper_bound = std::max(
                    precomputed_.successor_closure_capacity_bins[
                        static_cast<std::size_t>(item)],
                    tail_offset_[static_cast<std::size_t>(item)] + 1);
                if (static_cast<std::int64_t>(depth) +
                        prefix_bin_upper_bound + suffix_bin_upper_bound - 1 <=
                    lower_bound) {
                    continue;
                }
                const std::int64_t own_weight =
                    instance_.items[static_cast<std::size_t>(item)].weight;
                const std::int64_t prefix_weight = own_weight +
                    precomputed_.masked_weight_sum(
                        precomputed_.predecessor_closure_row(item), assigned);
                const std::int64_t suffix_weight = own_weight +
                    precomputed_.masked_weight_sum(
                        precomputed_.successor_closure_row(item), assigned);
                const int prefix_bins = std::max(
                    ceil_div_i64(prefix_weight, precomputed_.capacity),
                    earliest_offset_[static_cast<std::size_t>(item)] + 1);
                const int suffix_bins = std::max(
                    ceil_div_i64(suffix_weight, precomputed_.capacity),
                    tail_offset_[static_cast<std::size_t>(item)] + 1);
                lower_bound = std::max(
                    lower_bound,
                    depth + prefix_bins + suffix_bins - 1);
                if (lower_bound >= upper_bound_) {
                    break;
                }
            }
            if (lower_bound > bound_before_closure) {
                ++statistics_.closure_bound_improvements;
            }
            if (lower_bound >= upper_bound_) {
                ++statistics_.closure_bound_prunes;
            }
        }
        return strengthen_bin_packing_bounds(
            key, depth, assigned_count, lower_bound);
    }

    [[nodiscard]] int strengthen_bin_packing_bounds(
        const std::uint64_t* key,
        int depth,
        int assigned_count,
        int lower_bound) {
        if (bin_packing_ == nullptr ||
            lower_bound >= upper_bound_) {
            return lower_bound;
        }
        const int remaining_count = precomputed_.n - assigned_count;
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            conflict_remaining_[block] =
                precomputed_.all_mask[block] & ~key[block];
        }

        if (bin_packing_ != nullptr) {
            if (binlb_disabled_after_abort_) {
                ++statistics_.binlb_disabled_skips;
            } else if (binlb_budget_exhausted_) {
                ++statistics_.binlb_budget_skips;
            } else if (remaining_count > config_.bbr_binlb_max_items) {
                ++statistics_.binlb_item_skips;
            } else {
                const double remaining_budget =
                    config_.bbr_binlb_total_time_limit_seconds -
                    statistics_.binlb_seconds;
                if (remaining_budget <= 0.0) {
                    binlb_budget_exhausted_ = true;
                    ++statistics_.binlb_budget_skips;
                    return lower_bound;
                }
                ++statistics_.binlb_calls;
                const BinPackingBoundResult result = bin_packing_->solve(
                    conflict_remaining_.data(), deadline_, remaining_budget,
                    upper_bound_ - depth);
                if (!synchronize_memory_usage()) {
                    return upper_bound_;
                }
                statistics_.binlb_search_nodes += result.search_nodes;
                statistics_.binlb_loads += result.nondominated_loads;
                statistics_.binlb_memo_hits += result.memo_hits;
                statistics_.binlb_target_hits +=
                    result.target_reached ? 1U : 0U;
                statistics_.binlb_conflict_calls +=
                    result.conflict_phase_attempted ? 1U : 0U;
                statistics_.binlb_conflict_completed +=
                    result.conflict_phase_completed ? 1U : 0U;
                statistics_.binlb_conflict_search_nodes +=
                    result.conflict_search_nodes;
                statistics_.binlb_conflict_loads +=
                    result.conflict_maximal_loads;
                statistics_.binlb_conflict_memo_hits +=
                    result.conflict_memo_hits;
                statistics_.binlb_ordinary_memo_hits +=
                    result.ordinary_memo_hits;
                statistics_.binlb_seconds += result.seconds;
                if (statistics_.binlb_seconds >=
                    config_.bbr_binlb_total_time_limit_seconds) {
                    binlb_budget_exhausted_ = true;
                }
                statistics_.binlb_memo_entries =
                    bin_packing_->memo_entry_count();
                const int relaxation_bound = depth + result.lower_bound;
                if (relaxation_bound > lower_bound) {
                    ++statistics_.binlb_bound_improvements;
                    lower_bound = relaxation_bound;
                }
                if (lower_bound >= upper_bound_) {
                    ++statistics_.binlb_prunes;
                    return lower_bound;
                }
                if (result.item_limited) {
                    ++statistics_.binlb_item_skips;
                } else if (!result.completed &&
                           !result.useful_test_completed) {
                    if (result.attempted) {
                        ++statistics_.binlb_aborted;
                        statistics_.binlb_timeouts +=
                            result.timed_out ? 1U : 0U;
                        statistics_.binlb_node_limits +=
                            result.node_limited ? 1U : 0U;
                        statistics_.binlb_load_limits +=
                            result.load_limited ? 1U : 0U;
                        if (result.ordinary_phase_completed &&
                            result.conflict_aware) {
                            bin_packing_->disable_conflicts();
                        } else {
                            binlb_disabled_after_abort_ = true;
                        }
                    }
                } else {
                    ++statistics_.binlb_completed;
                }
            }
        }

        return lower_bound;
    }
    [[nodiscard]] int compute_lower_bound(
        const std::uint64_t* key,
        int depth,
        std::int64_t assigned_weight,
        int assigned_count,
        const std::int64_t* dff_sums) {
        int preliminary = compute_cheap_lower_bound(
            key, depth, assigned_weight, dff_sums);
        if (preliminary >= upper_bound_ ||
            assigned_weight == instance_.total_weight) {
            return preliminary;
        }
        preliminary = compute_precedence_lower_bound(
            key, depth, preliminary);
        if (preliminary >= upper_bound_) {
            return preliminary;
        }
        return strengthen_post_memory_lower_bound(
            key, depth, assigned_weight, assigned_count, preliminary);
    }

    void prepare_load_enumeration() {
        std::fill(load_mask_.begin(), load_mask_.end(), 0U);
        std::fill(load_dff_sums_.begin(), load_dff_sums_.end(), 0);
        std::fill(boundary_eligible_item_mask_.begin(),
                  boundary_eligible_item_mask_.end(), 0U);
        current_load_weight_ = 0;
        current_load_count_ = 0;
        load_assigned_delta_hash_ = 0U;
        current_remaining_capacity_ = precomputed_.capacity;
        root_has_ready_item_ = false;
        bbr12_root_eligible_size_ = 0U;

        const std::uint64_t* assigned = current_key_.data();
        const std::uint64_t* cooldown_one =
            precomputed_.cooldown_levels > 0
                ? current_key_.data() + precomputed_.blocks
                : nullptr;
        for (int item = 0; item < precomputed_.n; ++item) {
            if (bit_is_set(assigned, item) ||
                (cooldown_one != nullptr && bit_is_set(cooldown_one, item)) ||
                !mask_subset(precomputed_.pred_positive_row(item), assigned,
                             precomputed_.blocks)) {
                item_status_[static_cast<std::size_t>(item)] = 3U;
                remaining_zero_predecessors_[static_cast<std::size_t>(item)] = 0;
                continue;
            }
            item_status_[static_cast<std::size_t>(item)] = 0U;
            if (fast_bppp_dominance_masks_enabled_) {
                const unsigned value = static_cast<unsigned>(item);
                boundary_eligible_item_mask_[value >> 6U] |=
                    std::uint64_t{1} << (value & 63U);
            }
            int remaining = 0;
            const std::uint64_t* predecessors =
                precomputed_.pred_zero_row(item);
            for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
                remaining += std::popcount(predecessors[block] & ~assigned[block]);
            }
            remaining_zero_predecessors_[static_cast<std::size_t>(item)] =
                remaining;
            if (remaining == 0) {
                root_has_ready_item_ = true;
            }
        }
        for (const int item : bbr12_item_order_) {
            if (item_status_[static_cast<std::size_t>(item)] == 0U &&
                remaining_zero_predecessors_[static_cast<std::size_t>(item)] ==
                    0) {
                bbr12_eligible_[bbr12_root_eligible_size_++] = item;
            }
        }
    }

    void enumerate_loads() {
        load_dfs_bbr12(0U, bbr12_root_eligible_size_);
    }

    [[nodiscard]] bool stop_load_search() const noexcept {
        return stop_ || current_bound_ >= upper_bound_;
    }

    void load_dfs_bbr12(std::size_t start, std::size_t eligible_size) {
        if (stop_load_search()) {
            return;
        }
        ++statistics_.load_search_nodes;
        if ((statistics_.load_search_nodes & 2047U) == 0U) {
            if (deadline_.expired()) {
                timed_out_ = true;
                stop_ = true;
                return;
            }
        }

        bool extended = false;
        for (std::size_t position = start; position < eligible_size; ++position) {
            const int candidate = bbr12_eligible_[position];
            if (item_status_[static_cast<std::size_t>(candidate)] != 0U ||
                remaining_zero_predecessors_[static_cast<std::size_t>(candidate)] !=
                    0) {
                continue;
            }
            const int weight =
                instance_.items[static_cast<std::size_t>(candidate)].weight;
            if (weight > current_remaining_capacity_) {
                continue;
            }
            extended = true;
            item_status_[static_cast<std::size_t>(candidate)] = 1U;
            set_bit(load_mask_.data(), candidate);
            current_load_weight_ += weight;
            ++current_load_count_;
            load_assigned_delta_hash_ ^=
                precomputed_.item_hash[static_cast<std::size_t>(candidate)];
            current_remaining_capacity_ -= weight;
            const std::int64_t* item_dff =
                precomputed_.transformed_weights_for_item(candidate);
            for (std::size_t transform = 0;
                 transform < load_dff_sums_.size(); ++transform) {
                load_dff_sums_[transform] += item_dff[transform];
            }

            std::size_t next_eligible_size = eligible_size;
            const int successor_begin = bbr12_successor_offsets_[
                static_cast<std::size_t>(candidate)];
            const int successor_end = bbr12_successor_offsets_[
                static_cast<std::size_t>(candidate + 1)];
            for (int successor_position = successor_begin;
                 successor_position < successor_end; ++successor_position) {
                const int successor = bbr12_successors_[
                    static_cast<std::size_t>(successor_position)];
                if (item_status_[static_cast<std::size_t>(successor)] <= 2U) {
                    int& remaining = remaining_zero_predecessors_[
                        static_cast<std::size_t>(successor)];
                    --remaining;
                    if (item_status_[static_cast<std::size_t>(successor)] == 0U &&
                        remaining == 0) {
                        bbr12_eligible_[next_eligible_size++] = successor;
                    }
                }
            }

            load_dfs_bbr12(position + 1U, next_eligible_size);

            for (int successor_position = successor_begin;
                 successor_position < successor_end; ++successor_position) {
                const int successor = bbr12_successors_[
                    static_cast<std::size_t>(successor_position)];
                if (item_status_[static_cast<std::size_t>(successor)] <= 2U) {
                    ++remaining_zero_predecessors_[
                        static_cast<std::size_t>(successor)];
                }
            }
            for (std::size_t transform = 0;
                 transform < load_dff_sums_.size(); ++transform) {
                load_dff_sums_[transform] -= item_dff[transform];
            }
            current_remaining_capacity_ += weight;
            --current_load_count_;
            load_assigned_delta_hash_ ^=
                precomputed_.item_hash[static_cast<std::size_t>(candidate)];
            current_load_weight_ -= weight;
            clear_bit(load_mask_.data(), candidate);
            item_status_[static_cast<std::size_t>(candidate)] = 0U;

            if (stop_load_search()) {
                return;
            }
        }

        if (extended) {
            return;
        }
        for (std::size_t position = 0; position < eligible_size; ++position) {
            const int item = bbr12_eligible_[position];
            if (item_status_[static_cast<std::size_t>(item)] == 0U &&
                remaining_zero_predecessors_[static_cast<std::size_t>(item)] ==
                    0 &&
                instance_.items[static_cast<std::size_t>(item)].weight <=
                    current_remaining_capacity_) {
                ++statistics_.nonmaximal_load_prunes;
                return;
            }
        }
        emit_if_maximal();
    }

    void emit_if_maximal() {
        if (current_load_count_ == 0 && root_has_ready_item_) {
            ++statistics_.nonmaximal_load_prunes;
            return;
        }

        ++statistics_.loads_generated;
        if (config_.bbr_enable_jackson &&
            jackson_dominated()) {
            ++statistics_.jackson_prunes;
            return;
        }
        if (config_.bbr_enable_generalized_item_dominance &&
            (precomputed_.salbp_semantics || precomputed_.bppp_semantics) &&
            generalized_item_dominated()) {
            ++statistics_.generalized_item_dominance_prunes;
            return;
        }
        if (precomputed_.salbp_semantics && config_.bbr_enable_no_successor &&
            no_successor_dominated()) {
            ++statistics_.no_successor_prunes;
            return;
        }
        if (current_load_count_ == 0) {
            ++statistics_.forced_empty_transitions;
        }
        process_load();
    }

    [[nodiscard]] bool jackson_dominated() {
        if (!precomputed_.salbp_semantics ||
            !config_.bbr_enable_bbr12_jackson) {
            return replacement_dominated(precomputed_.jackson_dominators,
                                         false);
        }
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            std::uint64_t selected = load_mask_[block];
            while (selected != 0U) {
                const unsigned bit = std::countr_zero(selected);
                const int dominated = static_cast<int>(block * 64U + bit);
                const int dominated_weight = instance_.items[
                    static_cast<std::size_t>(dominated)].weight;
                for (const int dominator : precomputed_.jackson_dominators[
                         static_cast<std::size_t>(dominated)]) {
                    const unsigned char status =
                        item_status_[static_cast<std::size_t>(dominator)];
                    if ((status == 0U || status == 2U) &&
                        remaining_zero_predecessors_[
                            static_cast<std::size_t>(dominator)] == 0 &&
                        instance_.items[static_cast<std::size_t>(dominator)]
                                    .weight -
                                dominated_weight <=
                            current_remaining_capacity_) {
                        return true;
                    }
                }
                selected &= selected - 1U;
            }
        }
        return false;
    }

    [[nodiscard]] bool generalized_item_dominated() {
        if (fast_bppp_dominance_masks_enabled_) {
            return bppp_generalized_item_dominated();
        }
        return replacement_dominated(precomputed_.generalized_dominators,
                                     true);
    }

    [[nodiscard]] bool bppp_generalized_item_dominated() {
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            std::uint64_t selected = load_mask_[block];
            while (selected != 0U) {
                const unsigned bit = std::countr_zero(selected);
                const int dominated = static_cast<int>(block * 64U + bit);
                const std::uint64_t* dominators =
                    precomputed_.generalized_dominator_row(dominated);
                const int dominated_weight = instance_.items[
                    static_cast<std::size_t>(dominated)].weight;
                const int weight_limit = static_cast<int>(std::min<std::int64_t>(
                    precomputed_.capacity,
                    static_cast<std::int64_t>(current_remaining_capacity_) +
                        dominated_weight));
                const std::uint64_t* fit = fit_item_masks_.data() +
                    static_cast<std::size_t>(weight_limit) *
                        precomputed_.blocks;
                bool available = false;
                for (std::size_t word = 0; word < precomputed_.blocks; ++word) {
                    statistics_.generalized_item_dominance_checks +=
                        std::popcount(dominators[word]);
                    available = available ||
                        (dominators[word] & boundary_eligible_item_mask_[word] &
                         ~load_mask_[word] & fit[word]) != 0U;
                }
                if (available) {
                    return true;
                }
                selected &= selected - 1U;
            }
        }
        return false;
    }

    [[nodiscard]] bool replacement_dominated(
        const std::vector<std::vector<int>>& dominators,
        bool count_checks) {
        const std::uint64_t* assigned = current_key_.data();
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            std::uint64_t selected = load_mask_[block];
            while (selected != 0U) {
                const unsigned bit = std::countr_zero(selected);
                const int dominated =
                    static_cast<int>(block * 64U + bit);
                for (const int dominator : dominators[
                         static_cast<std::size_t>(dominated)]) {
                    if (count_checks) {
                        ++statistics_.generalized_item_dominance_checks;
                    }
                    if (bit_is_set(assigned, dominator) ||
                        bit_is_set(load_mask_.data(), dominator) ||
                        item_status_[static_cast<std::size_t>(dominator)] == 3U) {
                        continue;
                    }
                    const int replacement_weight =
                        current_load_weight_ -
                        instance_.items[
                            static_cast<std::size_t>(dominated)].weight +
                        instance_.items[
                            static_cast<std::size_t>(dominator)].weight;
                    if (replacement_weight > precomputed_.capacity) {
                        continue;
                    }

                    if (precomputed_.bppp_semantics) {
                        if (mask_subset(
                                precomputed_.pred_positive_row(dominator),
                                assigned, precomputed_.blocks)) {
                            return true;
                        }
                        continue;
                    }

                    for (std::size_t word = 0;
                         word < precomputed_.blocks; ++word) {
                        replacement_union_[word] =
                            assigned[word] | load_mask_[word];
                    }
                    clear_bit(replacement_union_.data(), dominated);
                    set_bit(replacement_union_.data(), dominator);

                    bool removed_item_is_required = false;
                    const std::uint64_t* zero_successors =
                        precomputed_.zero_successor_row(dominated);
                    for (std::size_t word = 0;
                         word < precomputed_.blocks; ++word) {
                        if ((zero_successors[word] & load_mask_[word]) != 0U) {
                            removed_item_is_required = true;
                            break;
                        }
                    }
                    if (!removed_item_is_required &&
                        mask_subset(precomputed_.pred_zero_row(dominator),
                                    replacement_union_.data(),
                                    precomputed_.blocks) &&
                        mask_subset(precomputed_.pred_positive_row(dominator),
                                    assigned, precomputed_.blocks)) {
                        return true;
                    }
                }
                selected &= selected - 1U;
            }
        }
        return false;
    }

    [[nodiscard]] bool no_successor_dominated() const noexcept {
        constexpr int kSuccessorMaskThreshold = 192;
        if (precomputed_.n <= kSuccessorMaskThreshold) {
            for (int item = 0; item < precomputed_.n; ++item) {
                if (bit_is_set(load_mask_.data(), item) &&
                    precomputed_.has_successor[
                        static_cast<std::size_t>(item)] != 0U) {
                    return false;
                }
            }
            const std::uint64_t* assigned = current_key_.data();
            for (int item = 0; item < precomputed_.n; ++item) {
                if (!bit_is_set(assigned, item) &&
                    !bit_is_set(load_mask_.data(), item) &&
                    precomputed_.has_successor[
                        static_cast<std::size_t>(item)] != 0U) {
                    return true;
                }
            }
            return false;
        }
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            if ((load_mask_[block] &
                 precomputed_.has_successor_mask[block]) != 0U) {
                return false;
            }
        }
        const std::uint64_t* assigned = current_key_.data();
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            if (((~assigned[block]) & precomputed_.all_mask[block] &
                 precomputed_.has_successor_mask[block]) != 0U) {
                return true;
            }
        }
        return false;
    }

    void process_load() {
        const int child_depth = current_depth_ + 1;
        std::uint64_t* child_assigned = child_key_.data();
        for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
            child_assigned[block] = current_key_[block] | load_mask_[block];
        }
        for (int level = 0; level < precomputed_.cooldown_levels; ++level) {
            std::uint64_t* destination = child_key_.data() +
                static_cast<std::size_t>(level + 1) * precomputed_.blocks;
            const std::uint64_t* shifted =
                level + 1 < precomputed_.cooldown_levels
                    ? current_key_.data() +
                          static_cast<std::size_t>(level + 2) *
                              precomputed_.blocks
                    : nullptr;
            for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
                destination[block] = shifted == nullptr ? 0U : shifted[block];
            }
            precomputed_.for_each_set_bit(load_mask_.data(), [&](int item) {
                const std::uint64_t* successors =
                    precomputed_.successors_above_row(level, item);
                for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
                    destination[block] |= successors[block];
                }
            });
            for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
                destination[block] &= ~child_assigned[block];
            }
        }

        const std::int64_t child_assigned_weight =
            current_assigned_weight_ + current_load_weight_;
        const int child_assigned_count =
            current_assigned_count_ + current_load_count_;
        for (std::size_t transform = 0;
             transform < child_dff_sums_.size(); ++transform) {
            child_dff_sums_[transform] =
                current_dff_sums_[transform] - load_dff_sums_[transform];
        }

        if (child_assigned_count == precomputed_.n) {
            if (child_depth < upper_bound_) {
                Assignment candidate = reconstruct_with_load();
                std::string diagnostic;
                if (!check_assignment(instance_, candidate, &diagnostic)) {
                    throw std::logic_error(
                        "BBR reconstructed an invalid incumbent: " + diagnostic);
                }
                incumbent_ = std::move(candidate);
                refresh_fixed_memory_bytes();
                upper_bound_ = child_depth;
                ++statistics_.incumbent_updates;
            }
            return;
        }
        if (child_depth + 1 >= upper_bound_) {
            ++statistics_.bound_prunes;
            return;
        }
        int child_bound = compute_cheap_lower_bound(
            child_key_.data(), child_depth, child_assigned_weight,
            child_dff_sums_.data());
        if (child_bound >= upper_bound_) {
            ++statistics_.bound_prunes;
            return;
        }
        const std::uint64_t child_assigned_hash =
            current_assigned_hash_ ^ load_assigned_delta_hash_;
        const std::uint64_t child_hash =
            precomputed_.salbp_semantics
                ? child_assigned_hash
                : hash_words(child_key_.data(), precomputed_.key_words);
        const std::uint32_t existing = exact_table_.find(child_key_.data(),
                                                         child_hash);
        if (existing != kInvalidState) {
            const std::size_t index = existing;
            if (store_.depths[index] <= child_depth) {
                ++statistics_.exact_memory_prunes;
                return;
            }
        }
        if (profile_dominance_enabled() &&
            profile_table_.dominated(child_key_.data(), child_assigned_hash,
                                     child_depth)) {
            ++statistics_.profile_dominance_prunes;
            return;
        }
        if (config_.bbr_enable_superset_memory &&
            superset_state_exists(child_key_.data(), child_assigned_hash,
                                  child_depth, child_assigned_weight)) {
            ++statistics_.superset_memory_prunes;
            return;
        }

        child_bound = compute_precedence_lower_bound(
            child_key_.data(), child_depth, child_bound);
        if (child_bound >= upper_bound_) {
            ++statistics_.bound_prunes;
            return;
        }

        child_bound = strengthen_post_memory_lower_bound(
            child_key_.data(), child_depth, child_assigned_weight,
            child_assigned_count, child_bound);
        if (child_bound >= upper_bound_) {
            ++statistics_.bound_prunes;
            return;
        }

        if (existing != kInvalidState) {
            const std::size_t index = existing;
            store_.parents[index] = current_state_;
            store_.depths[index] = child_depth;
            store_.bounds[index] = child_bound;
            store_.assigned_weights[index] = child_assigned_weight;
            store_.assigned_counts[index] = child_assigned_count;
            ++store_.versions[index];
            if (store_.queued[index] == 0U) {
                store_.queued[index] = 1U;
                ++open_states_;
            }
            if (!prepare_heap_push(child_depth)) {
                memory_limited_ = true;
                stop_ = true;
                return;
            }
            push_state_unchecked(existing);
            ++statistics_.states_reopened;
            statistics_.peak_open_states =
                std::max(statistics_.peak_open_states, open_states_);
            return;
        }

        if (!acquire_state_slot()) {
            return;
        }
        if (!prepare_exact_insert(profile_table_.memory_bytes() +
                                  heap_memory_bytes() +
                                  fixed_memory_bytes())) {
            memory_limited_ = true;
            stop_ = true;
            return;
        }
        if (profile_dominance_enabled() &&
            !prepare_profile_insert(
                child_key_.data(), child_assigned_hash,
                exact_table_.memory_bytes() + heap_memory_bytes() +
                    fixed_memory_bytes())) {
            memory_limited_ = true;
            stop_ = true;
            return;
        }
        if (!prepare_store_append(
                exact_table_.memory_bytes() + profile_table_.memory_bytes() +
                    heap_memory_bytes() + fixed_memory_bytes()) ||
            !prepare_heap_push(child_depth)) {
            memory_limited_ = true;
            stop_ = true;
            return;
        }
        const std::uint32_t child = store_.append(
            child_key_.data(), child_hash, child_assigned_hash, current_state_,
            child_depth, child_bound, child_assigned_weight,
            child_assigned_count, child_dff_sums_.data());
        exact_table_.insert(child);
        if (profile_dominance_enabled()) {
            profile_table_.insert(child);
        }
        ++statistics_.states_created;
        ++open_states_;
        statistics_.peak_open_states =
            std::max(statistics_.peak_open_states, open_states_);
        push_state_unchecked(child);
        if ((statistics_.states_created & 4095U) == 0U) {
            update_peak_memory();
        }
    }

    [[nodiscard]] bool superset_state_exists(const std::uint64_t* key,
                                              std::uint64_t assigned_hash,
                                              int depth,
                                              std::int64_t assigned_weight) {
        std::copy(key, key + precomputed_.key_words, probe_key_.begin());
        for (const int item : precomputed_.branch_order) {
            if (bit_is_set(key, item)) {
                continue;
            }
            if (!mask_subset(precomputed_.predecessor_closure_row(item), key,
                             precomputed_.blocks)) {
                continue;
            }
            if (assigned_weight +
                    instance_.items[static_cast<std::size_t>(item)].weight >
                static_cast<std::int64_t>(depth) * precomputed_.capacity) {
                continue;
            }
            set_bit(probe_key_.data(), item);
            const std::uint64_t superset_assigned_hash =
                assigned_hash ^
                precomputed_.item_hash[static_cast<std::size_t>(item)];
            bool dominated = false;
            if (precomputed_.cooldown_levels == 0) {
                const std::uint64_t hash = precomputed_.salbp_semantics
                    ? superset_assigned_hash
                    : hash_words(probe_key_.data(), precomputed_.key_words);
                const std::uint32_t state =
                    exact_table_.find(probe_key_.data(), hash);
                dominated = state != kInvalidState &&
                    store_.depths[static_cast<std::size_t>(state)] <= depth;
            } else if (profile_dominance_enabled()) {
                dominated = profile_table_.dominated(
                    probe_key_.data(), superset_assigned_hash, depth);
            }
            clear_bit(probe_key_.data(), item);
            if (dominated) {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] bool profile_dominance_enabled() const noexcept {
        return config_.bbr_enable_profile_dominance &&
               precomputed_.cooldown_levels > 0;
    }

    [[nodiscard]] Assignment reconstruct_with_load() const {
        Assignment assignment;
        assignment.bin_count = current_depth_ + 1;
        assignment.bin_of_item.assign(
            static_cast<std::size_t>(precomputed_.n), -1);

        std::vector<std::uint32_t> chain;
        for (std::uint32_t state = current_state_; state != kInvalidState;
             state = store_.parents[static_cast<std::size_t>(state)]) {
            chain.push_back(state);
        }
        std::reverse(chain.begin(), chain.end());
        for (std::size_t position = 1; position < chain.size(); ++position) {
            const std::uint32_t parent = chain[position - 1U];
            const std::uint32_t child = chain[position];
            const std::uint64_t* parent_assigned = store_.key(parent);
            const std::uint64_t* child_assigned = store_.key(child);
            const int bin = store_.depths[static_cast<std::size_t>(child)] - 1;
            for (std::size_t block = 0; block < precomputed_.blocks; ++block) {
                std::uint64_t value =
                    child_assigned[block] & ~parent_assigned[block];
                while (value != 0U) {
                    const unsigned bit = std::countr_zero(value);
                    const int item = static_cast<int>(block * 64U + bit);
                    assignment.bin_of_item[static_cast<std::size_t>(item)] = bin;
                    value &= value - 1U;
                }
            }
        }
        precomputed_.for_each_set_bit(load_mask_.data(), [&](int item) {
            assignment.bin_of_item[static_cast<std::size_t>(item)] = current_depth_;
        });
        return assignment;
    }

    [[nodiscard]] int certified_open_lower_bound() const noexcept {
        int lower_bound = upper_bound_;
        bool found_open = false;
        for (std::size_t state = 0; state < store_.size(); ++state) {
            if (store_.queued[state] != 0U) {
                lower_bound = std::min(lower_bound, store_.bounds[state]);
                found_open = true;
            }
        }
        return found_open ? std::max(initial_lower_bound_, lower_bound)
                          : initial_lower_bound_;
    }

    [[nodiscard]] std::uint64_t heap_memory_bytes() const noexcept {
        return heap_memory_bytes_cached_;
    }

    [[nodiscard]] std::uint64_t fixed_memory_bytes() const noexcept {
        return saturated_add(
            fixed_memory_bytes_cached_,
            bin_packing_ == nullptr ? 0U : bin_packing_->memory_bytes());
    }

    void refresh_fixed_memory_bytes() noexcept {
        std::uint64_t bytes = precomputed_.memory_bytes();
        const auto add = [&](std::uint64_t amount) {
            bytes = saturated_add(bytes, amount);
        };
        add(vector_memory_bytes(incumbent_.bin_of_item));
        add(vector_memory_bytes(current_dff_sums_));
        add(vector_memory_bytes(current_key_));
        add(vector_memory_bytes(child_key_));
        add(vector_memory_bytes(probe_key_));
        add(vector_memory_bytes(load_mask_));
        add(vector_memory_bytes(replacement_union_));
        add(vector_memory_bytes(item_status_));
        add(vector_memory_bytes(remaining_zero_predecessors_));
        add(vector_memory_bytes(bbr12_item_order_));
        add(vector_memory_bytes(bbr12_rank_));
        add(vector_memory_bytes(bbr12_eligible_));
        add(vector_memory_bytes(bbr12_successors_));
        add(vector_memory_bytes(bbr12_successor_offsets_));
        add(vector_memory_bytes(fit_item_masks_));
        add(vector_memory_bytes(boundary_eligible_item_mask_));
        add(vector_memory_bytes(earliest_offset_));
        add(vector_memory_bytes(tail_offset_));
        add(vector_memory_bytes(machine_order_));
        add(vector_memory_bytes(machine_bucket_weights_));
        add(vector_memory_bytes(child_dff_sums_));
        add(vector_memory_bytes(load_dff_sums_));
        add(vector_memory_bytes(conflict_remaining_));
        fixed_memory_bytes_cached_ = bytes;
    }

    [[nodiscard]] bool ensure_memory_capacity(
        std::uint64_t required_bytes) noexcept {
        return required_bytes <= memory_limit_bytes_;
    }

    [[nodiscard]] bool grow_memory_capacity() noexcept {
        return false;
    }

    [[nodiscard]] bool synchronize_memory_usage() noexcept {
        if (!ensure_memory_capacity(current_memory_bytes())) {
            memory_limited_ = true;
            stop_ = true;
            return false;
        }
        update_peak_memory();
        return true;
    }

    [[nodiscard]] bool prepare_store_append(
        std::uint64_t other_memory_bytes) {
        while (!store_.prepare_append(other_memory_bytes,
                                      memory_limit_bytes_)) {
            if (!grow_memory_capacity()) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] bool prepare_exact_insert(
        std::uint64_t other_memory_bytes) {
        while (!exact_table_.prepare_insert(other_memory_bytes)) {
            if (!grow_memory_capacity()) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] bool prepare_profile_insert(
        const std::uint64_t* key,
        std::uint64_t assigned_hash,
        std::uint64_t other_memory_bytes) {
        while (!profile_table_.prepare_insert(
            key, assigned_hash, other_memory_bytes)) {
            if (!grow_memory_capacity()) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] std::uint64_t current_memory_bytes() const noexcept {
        std::uint64_t bytes = fixed_memory_bytes();
        bytes = saturated_add(bytes, store_.memory_bytes());
        bytes = saturated_add(bytes, exact_table_.memory_bytes());
        bytes = saturated_add(bytes, profile_table_.memory_bytes());
        return saturated_add(bytes, heap_memory_bytes());
    }

    void update_peak_memory() noexcept {
        statistics_.peak_memory_bytes =
            std::max(statistics_.peak_memory_bytes, current_memory_bytes());
    }

    const PreparedInstance& prepared_;
    const Instance& instance_;
    BbrPrecomputed precomputed_;
    int initial_lower_bound_ = 0;
    const Config& config_;
    Deadline& deadline_;
    std::uint64_t maximum_memory_limit_bytes_ = 0U;
    BbrStatistics statistics_;
    StateStore store_;
    std::uint64_t memory_limit_bytes_ = 0;
    ExactStateTable exact_table_;
    AssignedProfileTable profile_table_;
    std::unique_ptr<BinPackingBound> bin_packing_;
    Assignment incumbent_;
    int upper_bound_ = 0;
    std::vector<std::vector<QueueEntry>> heaps_;
    std::uint64_t fixed_memory_bytes_cached_ = 0;
    std::uint64_t heap_memory_bytes_cached_ = 0;
    std::size_t next_depth_ = 0;
    std::uint64_t open_states_ = 0;

    bool stop_ = false;
    bool optimal_ = false;
    bool timed_out_ = false;
    bool memory_limited_ = false;
    bool binlb_disabled_after_abort_ = false;
    bool binlb_budget_exhausted_ = false;

    std::uint32_t current_state_ = kInvalidState;
    std::uint32_t current_version_ = 0;
    int current_depth_ = 0;
    int current_bound_ = 0;
    std::int64_t current_assigned_weight_ = 0;
    int current_assigned_count_ = 0;
    std::uint64_t current_assigned_hash_ = 0;
    std::vector<std::int64_t> current_dff_sums_;

    std::vector<std::uint64_t> current_key_;
    std::vector<std::uint64_t> child_key_;
    std::vector<std::uint64_t> probe_key_;
    std::vector<std::uint64_t> load_mask_;
    std::vector<std::uint64_t> replacement_union_;
    std::vector<unsigned char> item_status_;
    std::vector<int> remaining_zero_predecessors_;
    std::vector<int> bbr12_item_order_;
    std::vector<int> bbr12_rank_;
    std::vector<int> bbr12_eligible_;
    std::vector<int> bbr12_successors_;
    std::vector<int> bbr12_successor_offsets_;
    std::size_t bbr12_root_eligible_size_ = 0U;
    std::vector<std::uint64_t> fit_item_masks_;
    std::vector<std::uint64_t> boundary_eligible_item_mask_;
    bool fit_item_masks_available_ = false;
    bool fast_bppp_dominance_masks_enabled_ = false;
    std::vector<int> earliest_offset_;
    std::vector<int> tail_offset_;
    std::vector<int> machine_order_;
    std::vector<std::int64_t> machine_bucket_weights_;
    int remaining_machine_count_ = 0;
    int remaining_maximum_tail_ = 0;
    int remaining_maximum_earliest_ = 0;
    int remaining_minimum_tail_ = 0;
    std::vector<std::int64_t> child_dff_sums_;
    std::vector<std::int64_t> load_dff_sums_;
    std::vector<std::uint64_t> conflict_remaining_;
    std::int64_t next_queue_machine_numerator_ = 0;
    int next_queue_longest_path_tie_ = 0;

    int current_load_weight_ = 0;
    int current_load_count_ = 0;
    std::uint64_t load_assigned_delta_hash_ = 0;
    int current_remaining_capacity_ = 0;
    bool root_has_ready_item_ = false;
};

void add_phase_counters(BbrStatistics& target,
                        const BbrStatistics& source) noexcept {
    target.states_created += source.states_created;
    target.states_expanded += source.states_expanded;
    target.states_reopened += source.states_reopened;
    target.loads_generated += source.loads_generated;
    target.load_search_nodes += source.load_search_nodes;
    target.forced_empty_transitions += source.forced_empty_transitions;
    target.hash_lookups += source.hash_lookups;
    target.hash_probes += source.hash_probes;
    target.exact_memory_prunes += source.exact_memory_prunes;
    target.profile_dominance_prunes += source.profile_dominance_prunes;
    target.superset_memory_prunes += source.superset_memory_prunes;
    target.bound_prunes += source.bound_prunes;
    target.nonmaximal_load_prunes += source.nonmaximal_load_prunes;
    target.jackson_prunes += source.jackson_prunes;
    target.generalized_item_dominance_search_nodes +=
        source.generalized_item_dominance_search_nodes;
    target.generalized_item_dominance_checks +=
        source.generalized_item_dominance_checks;
    target.generalized_item_dominance_prunes +=
        source.generalized_item_dominance_prunes;
    target.no_successor_prunes += source.no_successor_prunes;
    target.machine_bound_calls += source.machine_bound_calls;
    target.machine_bound_prunes += source.machine_bound_prunes;
    target.closure_bound_calls += source.closure_bound_calls;
    target.closure_bound_improvements += source.closure_bound_improvements;
    target.closure_bound_prunes += source.closure_bound_prunes;
    target.binlb_calls += source.binlb_calls;
    target.binlb_completed += source.binlb_completed;
    target.binlb_aborted += source.binlb_aborted;
    target.binlb_timeouts += source.binlb_timeouts;
    target.binlb_node_limits += source.binlb_node_limits;
    target.binlb_load_limits += source.binlb_load_limits;
    target.binlb_item_skips += source.binlb_item_skips;
    target.binlb_disabled_skips += source.binlb_disabled_skips;
    target.binlb_budget_skips += source.binlb_budget_skips;
    target.binlb_bound_improvements += source.binlb_bound_improvements;
    target.binlb_prunes += source.binlb_prunes;
    target.binlb_search_nodes += source.binlb_search_nodes;
    target.binlb_loads += source.binlb_loads;
    target.binlb_memo_hits += source.binlb_memo_hits;
    target.binlb_memo_entries = std::max(
        target.binlb_memo_entries, source.binlb_memo_entries);
    target.binlb_target_hits += source.binlb_target_hits;
    target.binlb_conflict_edges = std::max(
        target.binlb_conflict_edges, source.binlb_conflict_edges);
    target.binlb_conflict_calls += source.binlb_conflict_calls;
    target.binlb_conflict_completed += source.binlb_conflict_completed;
    target.binlb_conflict_search_nodes +=
        source.binlb_conflict_search_nodes;
    target.binlb_conflict_loads += source.binlb_conflict_loads;
    target.binlb_conflict_memo_hits += source.binlb_conflict_memo_hits;
    target.binlb_ordinary_memo_hits += source.binlb_ordinary_memo_hits;
    target.incumbent_updates += source.incumbent_updates;
    target.peak_open_states =
        std::max(target.peak_open_states, source.peak_open_states);
    target.peak_memory_bytes =
        std::max(target.peak_memory_bytes, source.peak_memory_bytes);
    target.generalized_item_dominance_pairs = std::max(
        target.generalized_item_dominance_pairs,
        source.generalized_item_dominance_pairs);
    target.dff_transform_count = std::max(
        target.dff_transform_count, source.dff_transform_count);
    target.generalized_item_dominance_seconds +=
        source.generalized_item_dominance_seconds;
    target.binlb_seconds += source.binlb_seconds;
}

void restore_requested_configuration_metadata(
    BbrStatistics& statistics,
    const Config& requested_config) noexcept {
    statistics.item_dominance_enabled = requested_config.bbr_enable_jackson;
    statistics.generalized_item_dominance_enabled =
        requested_config.bbr_enable_generalized_item_dominance;
    statistics.paper_queue_order_enabled =
        requested_config.bbr_enable_paper_queue_order;
    statistics.complete_dff_enabled =
        requested_config.bbr_enable_complete_dff;
    statistics.binlb_enabled = requested_config.bbr_enable_binlb;
    statistics.conflict_binlb_enabled =
        requested_config.bbr_enable_binlb &&
        requested_config.bbr_enable_conflict_binlb;
}

void merge_preliminary_exact_statistics(
    BbrStatistics& exact,
    const BbrStatistics& preliminary,
    const Config& requested_config) noexcept {
    add_phase_counters(exact, preliminary);
    exact.attempted = exact.attempted || preliminary.attempted;
    exact.exact_phase_attempted = true;
    restore_requested_configuration_metadata(exact, requested_config);
    exact.search_seconds += preliminary.search_seconds;
    exact.exact_search_seconds += preliminary.exact_search_seconds;
}

[[nodiscard]] BbrResult run_exact_bbr(
    const PreparedInstance& prepared,
    int initial_lower_bound,
    const Assignment& initial_incumbent,
    const Config& config,
    Deadline& deadline) {
    BbrEngine engine(prepared, initial_lower_bound, initial_incumbent, config,
                     deadline);
    return engine.solve();
}

[[nodiscard]] bool complete_dff_applicable(
    const Instance& instance) noexcept {
    const bool all_zero = std::all_of(
        instance.arcs.begin(), instance.arcs.end(),
        [](const Arc& arc) { return arc.separation == 0; });
    const bool all_one = std::all_of(
        instance.arcs.begin(), instance.arcs.end(),
        [](const Arc& arc) { return arc.separation == 1; });
    return all_zero || all_one;
}

}

BbrResult run_branch_bound_remember(const PreparedInstance& prepared,
                                    int initial_lower_bound,
                                    const Config& config,
                                    Deadline& deadline) {
    if (prepared.fixed_bin_offset < 0) {
        throw std::invalid_argument("invalid prepared BBR fixed-bin offset");
    }
    const int residual_lower_bound =
        std::max(1, initial_lower_bound - prepared.fixed_bin_offset);
    const auto restore_full_bound = [&](BbrResult result) {
        restore_requested_configuration_metadata(result.statistics, config);
        result.certified_lower_bound += prepared.fixed_bin_offset;
        return result;
    };
    constexpr double kPreliminaryDffProbeSeconds = 0.001;
    if (config.bbr_enable_complete_dff &&
        complete_dff_applicable(prepared.search_instance) &&
        deadline.remaining_seconds() > 2.0 * kPreliminaryDffProbeSeconds) {
        Config preliminary_config = config;
        preliminary_config.bbr_enable_complete_dff = false;
        preliminary_config.bbr_enable_generalized_item_dominance = false;
        preliminary_config.bbr_enable_binlb = false;
        Deadline preliminary_deadline(std::min(
            kPreliminaryDffProbeSeconds, deadline.remaining_seconds()));
        BbrResult preliminary;
        {
            BbrEngine engine(prepared, residual_lower_bound,
                             prepared.search_incumbent, preliminary_config,
                             preliminary_deadline);
            preliminary = engine.solve();
        }
        preliminary.statistics.complete_dff_enabled =
            config.bbr_enable_complete_dff;
        preliminary.statistics.generalized_item_dominance_enabled =
            config.bbr_enable_generalized_item_dominance;
        if (preliminary.optimal || deadline.expired()) {
            if (!preliminary.optimal && deadline.expired()) {
                preliminary.timed_out = true;
                preliminary.memory_limited = false;
                preliminary.statistics.timed_out = true;
                preliminary.statistics.memory_limited = false;
            }
            return restore_full_bound(std::move(preliminary));
        }

        BbrResult exact = run_exact_bbr(
            prepared, residual_lower_bound, preliminary.incumbent,
            config, deadline);
        merge_preliminary_exact_statistics(
            exact.statistics, preliminary.statistics, config);
        return restore_full_bound(std::move(exact));
    }

    return restore_full_bound(run_exact_bbr(
        prepared, residual_lower_bound, prepared.search_incumbent, config,
        deadline));
}

}
