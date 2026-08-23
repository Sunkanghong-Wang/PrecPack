#include "precpack/types.hpp"

#include <algorithm>
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

    predecessors.assign(static_cast<std::size_t>(n), {});
    successors.assign(static_cast<std::size_t>(n), {});
    predecessor_arcs.assign(static_cast<std::size_t>(n), {});
    successor_arcs.assign(static_cast<std::size_t>(n), {});
    std::vector<int> indegree(static_cast<std::size_t>(n), 0);
    for (const Arc& arc : arcs) {
        predecessors[static_cast<std::size_t>(arc.to)].push_back(arc.from);
        successors[static_cast<std::size_t>(arc.from)].push_back(arc.to);
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
        for (const int next : successors[static_cast<std::size_t>(current)]) {
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

    constexpr int kUnreachable = std::numeric_limits<int>::min() / 4;
    longest_separation.assign(static_cast<std::size_t>(n) * n, 0);
    std::vector<int> distance(static_cast<std::size_t>(n), kUnreachable);
    for (int source = 0; source < n; ++source) {
        std::fill(distance.begin(), distance.end(), kUnreachable);
        distance[static_cast<std::size_t>(source)] = 0;
        for (const int current : topological_order) {
            if (distance[static_cast<std::size_t>(current)] == kUnreachable) {
                continue;
            }
            for (const auto& [next, arc_distance] :
                 successor_arcs[static_cast<std::size_t>(current)]) {
                distance[static_cast<std::size_t>(next)] =
                    std::max(distance[static_cast<std::size_t>(next)],
                             checked_path_distance(
                                 distance[static_cast<std::size_t>(current)],
                                 arc_distance));
            }
        }
        for (int target = 0; target < n; ++target) {
            if (distance[static_cast<std::size_t>(target)] != kUnreachable) {
                longest_separation[static_cast<std::size_t>(source) * n + target] =
                    distance[static_cast<std::size_t>(target)];
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
        case SolveStatus::kStateLimit:
            return "STATE_LIMIT";
        case SolveStatus::kMemoryLimit:
            return "MEMORY_LIMIT";
    }
    return "UNKNOWN";
}

}
