#include "bppp_flow_bound.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace precpack::internal {
namespace {

class DinicFlow {
public:
    explicit DinicFlow(int vertex_count)
        : adjacency_(static_cast<std::size_t>(vertex_count)),
          level_(static_cast<std::size_t>(vertex_count)),
          cursor_(static_cast<std::size_t>(vertex_count)) {
        if (vertex_count <= 0) {
            throw std::invalid_argument("flow network must contain vertices");
        }
    }

    void reserve(int vertex, std::size_t degree) {
        adjacency_[static_cast<std::size_t>(vertex)].reserve(degree);
    }

    void add_edge(int from, int to, std::int64_t capacity) {
        if (capacity < 0) {
            throw std::invalid_argument("flow capacity must be nonnegative");
        }
        const int reverse_from =
            static_cast<int>(adjacency_[static_cast<std::size_t>(to)].size());
        const int reverse_to =
            static_cast<int>(adjacency_[static_cast<std::size_t>(from)].size());
        adjacency_[static_cast<std::size_t>(from)].push_back(
            Edge{to, reverse_from, capacity});
        adjacency_[static_cast<std::size_t>(to)].push_back(
            Edge{from, reverse_to, 0});
    }

    [[nodiscard]] std::int64_t maximum_flow(int source, int sink) {
        std::int64_t total = 0;
        while (build_levels(source, sink)) {
            std::fill(cursor_.begin(), cursor_.end(), 0);
            while (const std::int64_t sent =
                       send(source, sink,
                            std::numeric_limits<std::int64_t>::max())) {
                total += sent;
            }
        }
        return total;
    }

private:
    struct Edge {
        int to = -1;
        int reverse = -1;
        std::int64_t capacity = 0;
    };

    [[nodiscard]] bool build_levels(int source, int sink) {
        std::fill(level_.begin(), level_.end(), -1);
        std::queue<int> queue;
        level_[static_cast<std::size_t>(source)] = 0;
        queue.push(source);
        while (!queue.empty()) {
            const int vertex = queue.front();
            queue.pop();
            for (const Edge& edge :
                 adjacency_[static_cast<std::size_t>(vertex)]) {
                if (edge.capacity <= 0 ||
                    level_[static_cast<std::size_t>(edge.to)] >= 0) {
                    continue;
                }
                level_[static_cast<std::size_t>(edge.to)] =
                    level_[static_cast<std::size_t>(vertex)] + 1;
                queue.push(edge.to);
            }
        }
        return level_[static_cast<std::size_t>(sink)] >= 0;
    }

    [[nodiscard]] std::int64_t send(int vertex,
                                    int sink,
                                    std::int64_t available) {
        if (vertex == sink) {
            return available;
        }
        auto& edges = adjacency_[static_cast<std::size_t>(vertex)];
        for (int& position = cursor_[static_cast<std::size_t>(vertex)];
             position < static_cast<int>(edges.size()); ++position) {
            Edge& edge = edges[static_cast<std::size_t>(position)];
            if (edge.capacity <= 0 ||
                level_[static_cast<std::size_t>(edge.to)] !=
                    level_[static_cast<std::size_t>(vertex)] + 1) {
                continue;
            }
            const std::int64_t sent =
                send(edge.to, sink, std::min(available, edge.capacity));
            if (sent <= 0) {
                continue;
            }
            edge.capacity -= sent;
            adjacency_[static_cast<std::size_t>(edge.to)]
                      [static_cast<std::size_t>(edge.reverse)]
                          .capacity += sent;
            return sent;
        }
        return 0;
    }

    std::vector<std::vector<Edge>> adjacency_;
    std::vector<int> level_;
    std::vector<int> cursor_;
};

[[nodiscard]] std::vector<int> longest_strict_path(const Instance& instance) {
    const int n = instance.size();
    int last = 0;
    for (int item = 1; item < n; ++item) {
        if (instance.front[static_cast<std::size_t>(item)] >
            instance.front[static_cast<std::size_t>(last)]) {
            last = item;
        }
    }

    std::vector<int> path;
    path.reserve(static_cast<std::size_t>(
        instance.front[static_cast<std::size_t>(last)] + 1));
    while (true) {
        path.push_back(last);
        if (instance.front[static_cast<std::size_t>(last)] == 0) {
            break;
        }
        int predecessor_on_path = -1;
        for (const auto& [predecessor, separation] :
             instance.predecessor_arcs[static_cast<std::size_t>(last)]) {
            if (separation == 1 &&
                instance.front[static_cast<std::size_t>(predecessor)] + 1 ==
                    instance.front[static_cast<std::size_t>(last)]) {
                predecessor_on_path = predecessor;
                break;
            }
        }
        if (predecessor_on_path < 0) {
            throw std::logic_error("BPP-P longest path reconstruction failed");
        }
        last = predecessor_on_path;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

[[nodiscard]] bool comparable(const Instance& instance, int lhs, int rhs) {
    return instance.has_positive_separation_path(lhs, rhs) ||
           instance.has_positive_separation_path(rhs, lhs);
}

[[nodiscard]] int checked_bound(std::size_t path_size,
                                std::int64_t unpacked_weight,
                                int capacity) {
    const std::int64_t extra_bins =
        unpacked_weight / capacity + (unpacked_weight % capacity != 0 ? 1 : 0);
    const std::int64_t bound =
        static_cast<std::int64_t>(path_size) + extra_bins;
    if (bound > std::numeric_limits<int>::max()) {
        throw std::overflow_error("BPP-P flow lower bound does not fit int");
    }
    return static_cast<int>(bound);
}

}  // namespace

int compute_bppp_flow_lower_bound(const Instance& instance) {
    if (instance.problem_type != "BPP-P") {
        throw std::invalid_argument("BPP-P flow bound requires a BPP-P instance");
    }
    if (instance.size() <= 0 || instance.capacity <= 0 ||
        instance.front.size() != instance.items.size() ||
        instance.reachability_blocks !=
            (instance.items.size() + 63U) / 64U ||
        instance.reachable.size() !=
            instance.items.size() * instance.reachability_blocks ||
        instance.reaching.size() != instance.reachable.size() ||
        instance.positive_reachable.size() != instance.reachable.size()) {
        throw std::invalid_argument(
            "BPP-P flow bound requires an initialized instance");
    }
    if (std::any_of(instance.arcs.begin(), instance.arcs.end(),
                    [](const Arc& arc) { return arc.separation != 1; })) {
        throw std::invalid_argument(
            "BPP-P flow bound requires unit-separation precedence arcs");
    }

    const std::vector<int> path = longest_strict_path(instance);
    std::vector<unsigned char> on_path(
        static_cast<std::size_t>(instance.size()), 0U);
    for (const int item : path) {
        on_path[static_cast<std::size_t>(item)] = 1U;
    }

    std::vector<int> remaining_items;
    remaining_items.reserve(
        static_cast<std::size_t>(instance.size()) - path.size());
    std::int64_t remaining_weight = 0;
    for (int item = 0; item < instance.size(); ++item) {
        if (on_path[static_cast<std::size_t>(item)] != 0U) {
            continue;
        }
        remaining_items.push_back(item);
        remaining_weight += instance.items[static_cast<std::size_t>(item)].weight;
    }

    const int source = 0;
    const int item_start = 1;
    const int bin_start = item_start + static_cast<int>(remaining_items.size());
    const int sink = bin_start + static_cast<int>(path.size());
    DinicFlow flow(sink + 1);
    flow.reserve(source, remaining_items.size());
    flow.reserve(sink, path.size());
    for (std::size_t position = 0; position < remaining_items.size(); ++position) {
        const int item_vertex = item_start + static_cast<int>(position);
        flow.reserve(item_vertex, path.size() + 1U);
    }
    for (std::size_t bin = 0; bin < path.size(); ++bin) {
        flow.reserve(
            bin_start + static_cast<int>(bin), remaining_items.size() + 1U);
    }
    for (std::size_t position = 0; position < remaining_items.size(); ++position) {
        const int item_vertex = item_start + static_cast<int>(position);
        const int item = remaining_items[position];
        const int weight = instance.items[static_cast<std::size_t>(item)].weight;
        flow.add_edge(source, item_vertex, weight);
        for (std::size_t bin = 0; bin < path.size(); ++bin) {
            const int path_item = path[bin];
            const int path_weight =
                instance.items[static_cast<std::size_t>(path_item)].weight;
            if (comparable(instance, item, path_item) ||
                static_cast<std::int64_t>(weight) + path_weight > instance.capacity) {
                continue;
            }
            flow.add_edge(item_vertex, bin_start + static_cast<int>(bin), weight);
        }
    }
    for (std::size_t bin = 0; bin < path.size(); ++bin) {
        const int bin_vertex = bin_start + static_cast<int>(bin);
        const int path_item = path[bin];
        flow.add_edge(bin_vertex, sink,
                      instance.capacity -
                          instance.items[static_cast<std::size_t>(path_item)].weight);
    }

    const std::int64_t packed_weight = flow.maximum_flow(source, sink);
    if (packed_weight < 0 || packed_weight > remaining_weight) {
        throw std::logic_error("BPP-P flow bound produced an invalid flow value");
    }
    return checked_bound(
        path.size(), remaining_weight - packed_weight, instance.capacity);
}

}  // namespace precpack::internal
