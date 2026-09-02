#include "precpack/types.hpp"

#include <algorithm>
#include <bit>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace precpack {
namespace {

[[nodiscard]] std::uint64_t arc_key(int from, int to) noexcept {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(from)) << 32U) |
           static_cast<std::uint32_t>(to);
}

[[nodiscard]] int checked_path_distance(int lhs, int rhs) {
    const std::int64_t value =
        static_cast<std::int64_t>(lhs) + static_cast<std::int64_t>(rhs);
    constexpr std::int64_t kMaximumSupportedBinIndex =
        static_cast<std::int64_t>(std::numeric_limits<int>::max()) - 1;
    if (value > kMaximumSupportedBinIndex) {
        throw std::invalid_argument(
            "precedence path exceeds the supported bin-position range");
    }
    return static_cast<int>(value);
}

}

bool Assignment::complete() const noexcept {
    return bin_count > 0 &&
           std::all_of(bin_of_item.begin(), bin_of_item.end(),
                       [](int bin) { return bin >= 0; });
}

void Instance::initialize() {
    const int n = size();
    if (n <= 0) {
        throw std::invalid_argument("instance has no items");
    }
    if (capacity <= 0) {
        throw std::invalid_argument("capacity must be positive");
    }

    total_weight = 0;
    for (int i = 0; i < n; ++i) {
        Item& item = items[static_cast<std::size_t>(i)];
        if (item.original_index < 0) {
            item.original_index = i;
        }
        if (item.weight <= 0 || item.weight > capacity) {
            throw std::invalid_argument("item " + std::to_string(i + 1) +
                                        " has a nonpositive weight or exceeds capacity");
        }
        total_weight += item.weight;
    }

    std::unordered_map<std::uint64_t, int> strongest;
    strongest.reserve(arcs.size() * 2U + 1U);
    for (const Arc& arc : arcs) {
        if (arc.from < 0 || arc.from >= n || arc.to < 0 || arc.to >= n) {
            throw std::invalid_argument("precedence arc endpoint is out of range");
        }
        if (arc.from == arc.to) {
            throw std::invalid_argument("precedence graph contains a self-loop");
        }
        if (arc.separation < 0) {
            throw std::invalid_argument("precedence separation must be nonnegative");
        }
        const auto key = arc_key(arc.from, arc.to);
        auto [it, inserted] = strongest.emplace(key, arc.separation);
        if (!inserted) {
            it->second = std::max(it->second, arc.separation);
        }
    }

    arcs.clear();
    arcs.reserve(strongest.size());
    for (const auto& [key, separation_value] : strongest) {
        arcs.push_back(Arc{static_cast<int>(key >> 32U),
                           static_cast<int>(key & 0xffffffffU),
                           separation_value});
    }
    std::sort(arcs.begin(), arcs.end(), [](const Arc& lhs, const Arc& rhs) {
        return std::tie(lhs.from, lhs.to, lhs.separation) <
               std::tie(rhs.from, rhs.to, rhs.separation);
    });

    predecessor_arcs.assign(static_cast<std::size_t>(n), {});
    successor_arcs.assign(static_cast<std::size_t>(n), {});
    std::vector<int> indegree(static_cast<std::size_t>(n), 0);
    for (const Arc& arc : arcs) {
        predecessor_arcs[static_cast<std::size_t>(arc.to)].emplace_back(
            arc.from, arc.separation);
        successor_arcs[static_cast<std::size_t>(arc.from)].emplace_back(
            arc.to, arc.separation);
        ++indegree[static_cast<std::size_t>(arc.to)];
    }

    std::priority_queue<int, std::vector<int>, std::greater<>> ready;
    for (int i = 0; i < n; ++i) {
        if (indegree[static_cast<std::size_t>(i)] == 0) {
            ready.push(i);
        }
    }
    topological_order.clear();
    topological_order.reserve(static_cast<std::size_t>(n));
    while (!ready.empty()) {
        const int current = ready.top();
        ready.pop();
        topological_order.push_back(current);
        for (const auto& [next, separation_value] :
             successor_arcs[static_cast<std::size_t>(current)]) {
            static_cast<void>(separation_value);
            if (--indegree[static_cast<std::size_t>(next)] == 0) {
                ready.push(next);
            }
        }
    }
    if (static_cast<int>(topological_order.size()) != n) {
        throw std::invalid_argument("precedence graph is cyclic");
    }

    front.assign(static_cast<std::size_t>(n), 0);
    for (const int current : topological_order) {
        for (const auto& [next, distance] :
             successor_arcs[static_cast<std::size_t>(current)]) {
            front[static_cast<std::size_t>(next)] =
                std::max(front[static_cast<std::size_t>(next)],
                         checked_path_distance(
                             front[static_cast<std::size_t>(current)],
                             distance));
        }
    }

    back.assign(static_cast<std::size_t>(n), 0);
    for (auto it = topological_order.rbegin(); it != topological_order.rend(); ++it) {
        const int current = *it;
        for (const auto& [next, distance] :
             successor_arcs[static_cast<std::size_t>(current)]) {
            back[static_cast<std::size_t>(current)] =
                std::max(back[static_cast<std::size_t>(current)],
                         checked_path_distance(
                             distance,
                             back[static_cast<std::size_t>(next)]));
        }
    }

    reachability_blocks = (static_cast<std::size_t>(n) + 63U) / 64U;
    const std::size_t reachability_words =
        static_cast<std::size_t>(n) * reachability_blocks;
    reachable.assign(reachability_words, 0U);
    reaching.assign(reachability_words, 0U);
    positive_reachable.assign(reachability_words, 0U);
    const auto set_bit = [](std::uint64_t* row, int item) noexcept {
        row[static_cast<std::size_t>(item) / 64U] |=
            std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
    };
    for (auto order = topological_order.rbegin();
         order != topological_order.rend(); ++order) {
        const int current = *order;
        std::uint64_t* reachable_row = reachable.data() +
            static_cast<std::size_t>(current) * reachability_blocks;
        std::uint64_t* positive_row = positive_reachable.data() +
            static_cast<std::size_t>(current) * reachability_blocks;
        for (const auto& [next, separation_value] :
             successor_arcs[static_cast<std::size_t>(current)]) {
            const std::uint64_t* next_reachable = reachable.data() +
                static_cast<std::size_t>(next) * reachability_blocks;
            const std::uint64_t* positive_source =
                separation_value > 0
                    ? next_reachable
                    : positive_reachable.data() +
                          static_cast<std::size_t>(next) * reachability_blocks;
            for (std::size_t block = 0; block < reachability_blocks; ++block) {
                reachable_row[block] |= next_reachable[block];
                positive_row[block] |= positive_source[block];
            }
            set_bit(reachable_row, next);
            if (separation_value > 0) {
                set_bit(positive_row, next);
            }
        }
    }
    for (int from = 0; from < n; ++from) {
        const std::uint64_t* row = reachable.data() +
            static_cast<std::size_t>(from) * reachability_blocks;
        for (std::size_t block = 0; block < reachability_blocks; ++block) {
            std::uint64_t value = row[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                const int to = static_cast<int>(block * 64U + bit);
                if (to < n) {
                    set_bit(reaching.data() +
                                static_cast<std::size_t>(to) *
                                    reachability_blocks,
                            from);
                }
                value &= value - 1U;
            }
        }
    }
    for (int item = 0; item < n; ++item) {
        const std::int64_t required_bin_count =
            static_cast<std::int64_t>(
                front[static_cast<std::size_t>(item)]) +
            back[static_cast<std::size_t>(item)] + 1;
        if (required_bin_count > std::numeric_limits<int>::max()) {
            throw std::invalid_argument(
                "precedence path exceeds the supported bin-position range");
        }
    }
}

const char* to_string(SolveStatus status) noexcept {
    switch (status) {
        case SolveStatus::kNotSolved:
            return "NOT_SOLVED";
        case SolveStatus::kOptimal:
            return "OPTIMAL";
        case SolveStatus::kTimeLimit:
            return "TIME_LIMIT";
        case SolveStatus::kMemoryLimit:
            return "MEMORY_LIMIT";
    }
    return "UNKNOWN";
}

}
