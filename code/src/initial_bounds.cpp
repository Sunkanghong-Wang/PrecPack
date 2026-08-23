#include "precpack/initial_bounds.hpp"

#include "precpack/algorithms.hpp"
#include "precpack/bbr.hpp"
#include "precpack/exact_arithmetic.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace precpack {
namespace {

using Clock = std::chrono::steady_clock;

[[nodiscard]] int ceil_div(std::int64_t numerator, int denominator) {
    if (numerator < 0 || denominator <= 0) {
        throw std::invalid_argument("ceil_div requires a nonnegative numerator");
    }
    return static_cast<int>((numerator + denominator - 1) / denominator);
}

[[nodiscard]] int checked_position_distance(int lhs, int rhs) {
    const std::int64_t value =
        static_cast<std::int64_t>(lhs) + static_cast<std::int64_t>(rhs);
    constexpr std::int64_t kMaximumSupportedBinIndex =
        static_cast<std::int64_t>(std::numeric_limits<int>::max()) - 1;
    if (value > kMaximumSupportedBinIndex) {
        throw std::invalid_argument(
            "preprocessed precedence path exceeds the supported "
            "bin-position range");
    }
    return static_cast<int>(value);
}

class JavaRandom {
public:
    explicit JavaRandom(std::int64_t seed)
        : state_((static_cast<std::uint64_t>(seed) ^ kMultiplier) & kMask) {}

    [[nodiscard]] int next_int(int bound) {
        if (bound <= 0) {
            throw std::invalid_argument("JavaRandom bound must be positive");
        }
        if ((bound & -bound) == bound) {
            return static_cast<int>((static_cast<std::int64_t>(bound) * next_bits(31)) >> 31);
        }
        int bits = 0;
        int value = 0;
        do {
            bits = next_bits(31);
            value = bits % bound;
        } while (static_cast<std::int64_t>(bits) - value + (bound - 1) >=
                 (std::int64_t{1} << 31));
        return value;
    }

    void shuffle(std::vector<int>& values) {
        for (int size = static_cast<int>(values.size()); size > 1; --size) {
            std::swap(values[static_cast<std::size_t>(size - 1)],
                      values[static_cast<std::size_t>(next_int(size))]);
        }
    }

private:
    [[nodiscard]] int next_bits(int bits) {
        state_ = (state_ * kMultiplier + kAddend) & kMask;
        return static_cast<int>(state_ >> (48 - bits));
    }

    static constexpr std::uint64_t kMultiplier = 0x5DEECE66DULL;
    static constexpr std::uint64_t kAddend = 0xBULL;
    static constexpr std::uint64_t kMask = (std::uint64_t{1} << 48U) - 1U;
    std::uint64_t state_ = 0;
};

[[nodiscard]] std::int32_t java_string_hash(const std::string& value) noexcept {
    std::uint32_t hash = 0;
    for (const unsigned char character : value) {
        hash = 31U * hash + character;
    }
    return static_cast<std::int32_t>(hash);
}

struct WorkItem {
    int original_index = -1;
    int weight = 0;
};

struct WorkInstance {
    explicit WorkInstance(const Instance& original)
        : problem_type(original.problem_type),
          id(original.id),
          capacity(original.capacity),
          n(original.size()) {
        items.reserve(original.items.size());
        for (int i = 0; i < n; ++i) {
            items.push_back(WorkItem{
                i, original.items[static_cast<std::size_t>(i)].weight});
        }
        arcs.reserve(original.arcs.size());
        for (const Arc& arc : original.arcs) {
            arcs.push_back(arc);
        }
        rebuild_graph();
        if (!topologically_numbered()) {
            static_cast<void>(reorder_items());
        }
    }

    [[nodiscard]] bool topologically_numbered() const noexcept {
        return std::all_of(arcs.begin(), arcs.end(), [](const Arc& arc) {
            return arc.from < arc.to;
        });
    }

    [[nodiscard]] int edge_position(int from, int to) const noexcept {
        return edge_index[static_cast<std::size_t>(from) * n + to];
    }

    [[nodiscard]] int edge_separation(int from, int to) const noexcept {
        const int position = edge_position(from, to);
        return position < 0 ? 0 : arcs[static_cast<std::size_t>(position)].separation;
    }

    void rebuild_graph() {
        edge_index.assign(static_cast<std::size_t>(n) * n, -1);
        std::vector<int> out_count(static_cast<std::size_t>(n), 0);
        std::vector<int> in_count(static_cast<std::size_t>(n), 0);
        for (std::size_t index = 0; index < arcs.size(); ++index) {
            const Arc& arc = arcs[index];
            if (arc.from < 0 || arc.from >= n || arc.to < 0 || arc.to >= n ||
                arc.from == arc.to || arc.separation < 0) {
                throw std::logic_error("invalid arc in initial-bound working instance");
            }
            const std::size_t key = static_cast<std::size_t>(arc.from) * n + arc.to;
            if (edge_index[key] >= 0) {
                throw std::logic_error("duplicate working arc");
            }
            edge_index[key] = static_cast<int>(index);
            ++out_count[static_cast<std::size_t>(arc.from)];
            ++in_count[static_cast<std::size_t>(arc.to)];
        }

        succ_offset.assign(static_cast<std::size_t>(n + 1), 0);
        pred_offset.assign(static_cast<std::size_t>(n + 1), 0);
        for (int i = 0; i < n; ++i) {
            succ_offset[static_cast<std::size_t>(i + 1)] =
                succ_offset[static_cast<std::size_t>(i)] +
                out_count[static_cast<std::size_t>(i)];
            pred_offset[static_cast<std::size_t>(i + 1)] =
                pred_offset[static_cast<std::size_t>(i)] +
                in_count[static_cast<std::size_t>(i)];
        }
        succ_to.resize(arcs.size());
        succ_arc.resize(arcs.size());
        pred_from.resize(arcs.size());
        pred_arc.resize(arcs.size());
        std::vector<int> succ_cursor = succ_offset;
        std::vector<int> pred_cursor = pred_offset;
        for (std::size_t index = 0; index < arcs.size(); ++index) {
            const Arc& arc = arcs[index];
            const int out = succ_cursor[static_cast<std::size_t>(arc.from)]++;
            succ_to[static_cast<std::size_t>(out)] = arc.to;
            succ_arc[static_cast<std::size_t>(out)] = static_cast<int>(index);
            const int in = pred_cursor[static_cast<std::size_t>(arc.to)]++;
            pred_from[static_cast<std::size_t>(in)] = arc.from;
            pred_arc[static_cast<std::size_t>(in)] = static_cast<int>(index);
        }

        std::priority_queue<int, std::vector<int>, std::greater<>> ready;
        std::vector<int> indegree = in_count;
        topological_order.clear();
        topological_order.reserve(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            if (indegree[static_cast<std::size_t>(i)] == 0) {
                ready.push(i);
            }
        }
        while (!ready.empty()) {
            const int item = ready.top();
            ready.pop();
            topological_order.push_back(item);
            for (int p = succ_offset[static_cast<std::size_t>(item)];
                 p < succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
                const int next = succ_to[static_cast<std::size_t>(p)];
                if (--indegree[static_cast<std::size_t>(next)] == 0) {
                    ready.push(next);
                }
            }
        }
        if (static_cast<int>(topological_order.size()) != n) {
            throw std::logic_error("preprocessed precedence graph became cyclic");
        }
        recompute_position_bounds();
        total_weight = 0;
        for (const WorkItem& item : items) {
            total_weight += item.weight;
        }
    }

    void recompute_position_bounds() {
        front.assign(static_cast<std::size_t>(n), 0);
        for (const int item : topological_order) {
            for (int p = succ_offset[static_cast<std::size_t>(item)];
                 p < succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
                const int next = succ_to[static_cast<std::size_t>(p)];
                const int separation =
                    arcs[static_cast<std::size_t>(succ_arc[static_cast<std::size_t>(p)])]
                        .separation;
                front[static_cast<std::size_t>(next)] =
                    std::max(front[static_cast<std::size_t>(next)],
                             checked_position_distance(
                                 front[static_cast<std::size_t>(item)],
                                 separation));
            }
        }
        back.assign(static_cast<std::size_t>(n), 0);
        for (auto it = topological_order.rbegin(); it != topological_order.rend(); ++it) {
            const int item = *it;
            for (int p = succ_offset[static_cast<std::size_t>(item)];
                 p < succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
                const int next = succ_to[static_cast<std::size_t>(p)];
                const int separation =
                    arcs[static_cast<std::size_t>(succ_arc[static_cast<std::size_t>(p)])]
                        .separation;
                back[static_cast<std::size_t>(item)] =
                    std::max(back[static_cast<std::size_t>(item)],
                             checked_position_distance(
                                 separation,
                                 back[static_cast<std::size_t>(next)]));
            }
        }
        for (int item = 0; item < n; ++item) {
            const std::int64_t required_bin_count =
                static_cast<std::int64_t>(
                    front[static_cast<std::size_t>(item)]) +
                back[static_cast<std::size_t>(item)] + 1;
            if (required_bin_count > std::numeric_limits<int>::max()) {
                throw std::invalid_argument(
                    "preprocessed precedence path exceeds the supported "
                    "bin-position range");
            }
        }
    }

    [[nodiscard]] bool reorder_items() {
        std::vector<int> indegree(static_cast<std::size_t>(n), 0);
        for (int item = 0; item < n; ++item) {
            indegree[static_cast<std::size_t>(item)] =
                pred_offset[static_cast<std::size_t>(item + 1)] -
                pred_offset[static_cast<std::size_t>(item)];
        }
        std::vector<int> old_to_new(static_cast<std::size_t>(n), -1);
        int next_index = 0;
        while (next_index < n) {
            std::vector<int> ready;
            ready.reserve(static_cast<std::size_t>(n - next_index));
            for (int item = 0; item < n; ++item) {
                if (old_to_new[static_cast<std::size_t>(item)] < 0 &&
                    indegree[static_cast<std::size_t>(item)] == 0) {
                    ready.push_back(item);
                }
            }
            if (ready.empty()) {
                throw std::logic_error("no ready item while reordering initial instance");
            }
            std::stable_sort(ready.begin(), ready.end(), [&](int lhs, int rhs) {
                if (items[static_cast<std::size_t>(lhs)].weight !=
                    items[static_cast<std::size_t>(rhs)].weight) {
                    return items[static_cast<std::size_t>(lhs)].weight >
                           items[static_cast<std::size_t>(rhs)].weight;
                }
                const int lhs_successors =
                    succ_offset[static_cast<std::size_t>(lhs + 1)] -
                    succ_offset[static_cast<std::size_t>(lhs)];
                const int rhs_successors =
                    succ_offset[static_cast<std::size_t>(rhs + 1)] -
                    succ_offset[static_cast<std::size_t>(rhs)];
                return lhs_successors > rhs_successors;
            });
            for (const int item : ready) {
                old_to_new[static_cast<std::size_t>(item)] = next_index++;
                for (int p = succ_offset[static_cast<std::size_t>(item)];
                     p < succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
                    --indegree[static_cast<std::size_t>(succ_to[static_cast<std::size_t>(p)])];
                }
            }
        }
        bool changed = false;
        for (int item = 0; item < n; ++item) {
            changed = changed || old_to_new[static_cast<std::size_t>(item)] != item;
        }
        if (!changed) {
            return false;
        }
        std::vector<WorkItem> reordered(static_cast<std::size_t>(n));
        for (int old = 0; old < n; ++old) {
            reordered[static_cast<std::size_t>(old_to_new[static_cast<std::size_t>(old)])] =
                items[static_cast<std::size_t>(old)];
        }
        items = std::move(reordered);
        for (Arc& arc : arcs) {
            arc.from = old_to_new[static_cast<std::size_t>(arc.from)];
            arc.to = old_to_new[static_cast<std::size_t>(arc.to)];
        }
        rebuild_graph();
        if (!topologically_numbered()) {
            throw std::logic_error("item reordering failed to produce a topological numbering");
        }
        return true;
    }

    std::string problem_type;
    int id = -1;
    int capacity = 0;
    int n = 0;
    bool reversed = false;
    std::vector<WorkItem> items;
    std::vector<Arc> arcs;
    std::int64_t total_weight = 0;
    std::vector<int> edge_index;
    std::vector<int> succ_offset;
    std::vector<int> succ_to;
    std::vector<int> succ_arc;
    std::vector<int> pred_offset;
    std::vector<int> pred_from;
    std::vector<int> pred_arc;
    std::vector<int> topological_order;
    std::vector<int> front;
    std::vector<int> back;
};

class JavaPathCollector {
public:
    JavaPathCollector(const std::vector<std::vector<int>>& successors,
                      std::vector<std::uint32_t>& visited,
                      std::uint32_t stamp,
                      int end,
                      std::vector<int>& nodes)
        : successors_(successors),
          visited_(visited),
          stamp_(stamp),
          end_(end),
          nodes_(nodes) {}

    [[nodiscard]] bool find(int current) {
        bool reaches_end = false;
        for (const int next : successors_[static_cast<std::size_t>(current)]) {
            if (next > end_) {
                break;
            }
            if (next == end_) {
                reaches_end = true;
                break;
            }
            if (visited_[static_cast<std::size_t>(next)] != stamp_ && find(next)) {
                nodes_.push_back(next);
                reaches_end = true;
            }
        }
        visited_[static_cast<std::size_t>(current)] = stamp_;
        return reaches_end;
    }

private:
    const std::vector<std::vector<int>>& successors_;
    std::vector<std::uint32_t>& visited_;
    std::uint32_t stamp_;
    int end_;
    std::vector<int>& nodes_;
};

[[nodiscard]] std::vector<std::uint64_t> compute_reachability(
    const WorkInstance& instance) {
    const std::size_t blocks = (static_cast<std::size_t>(instance.n) + 63U) / 64U;
    std::vector<std::uint64_t> reachable(static_cast<std::size_t>(instance.n) * blocks,
                                         0U);
    for (auto order = instance.topological_order.rbegin();
         order != instance.topological_order.rend(); ++order) {
        const int item = *order;
        std::uint64_t* row = reachable.data() + static_cast<std::size_t>(item) * blocks;
        for (int p = instance.succ_offset[static_cast<std::size_t>(item)];
             p < instance.succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
            const int next = instance.succ_to[static_cast<std::size_t>(p)];
            row[static_cast<std::size_t>(next) / 64U] |=
                std::uint64_t{1} << (static_cast<unsigned>(next) & 63U);
            const std::uint64_t* next_row =
                reachable.data() + static_cast<std::size_t>(next) * blocks;
            for (std::size_t block = 0; block < blocks; ++block) {
                row[block] |= next_row[block];
            }
        }
    }
    return reachable;
}

[[nodiscard]] bool enhance_precedence_graph(WorkInstance& instance) {
    bool changed = false;
    const std::size_t blocks =
        (static_cast<std::size_t>(instance.n) + 63U) / 64U;
    const std::vector<std::uint64_t> reachable = compute_reachability(instance);
    std::vector<std::vector<int>> successors(static_cast<std::size_t>(instance.n));
    for (int item = 0; item < instance.n; ++item) {
        auto& list = successors[static_cast<std::size_t>(item)];
        list.reserve(static_cast<std::size_t>(
            instance.succ_offset[static_cast<std::size_t>(item + 1)] -
            instance.succ_offset[static_cast<std::size_t>(item)]));
        for (int p = instance.succ_offset[static_cast<std::size_t>(item)];
             p < instance.succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
            list.push_back(instance.succ_to[static_cast<std::size_t>(p)]);
        }
    }
    std::vector<std::uint32_t> visited(static_cast<std::size_t>(instance.n), 0U);
    std::uint32_t stamp = 0U;
    std::vector<int> path_nodes;
    path_nodes.reserve(static_cast<std::size_t>(instance.n));
    for (int i = 0; i < instance.n; ++i) {
        const std::uint64_t* reachable_i =
            reachable.data() + static_cast<std::size_t>(i) * blocks;
        for (int j = i + 1; j < instance.n; ++j) {
            const int position = instance.edge_position(i, j);
            if (position >= 0) {
                Arc& arc = instance.arcs[static_cast<std::size_t>(position)];
                const int strengthened =
                    ceil_div(static_cast<std::int64_t>(
                                 instance.items[static_cast<std::size_t>(i)].weight) +
                                 instance.items[static_cast<std::size_t>(j)].weight,
                             instance.capacity) -
                    1;
                if (strengthened > arc.separation) {
                    arc.separation = strengthened;
                    changed = true;
                }
                continue;
            }
            if (((reachable_i[static_cast<std::size_t>(j) / 64U] >>
                  (static_cast<unsigned>(j) & 63U)) &
                 1U) == 0U) {
                continue;
            }
            if (++stamp == 0U) {
                std::fill(visited.begin(), visited.end(), 0U);
                stamp = 1U;
            }
            path_nodes.clear();
            JavaPathCollector collector(successors, visited, stamp, j, path_nodes);
            static_cast<void>(collector.find(i));
            if (path_nodes.empty()) {
                continue;
            }
            std::int64_t total =
                static_cast<std::int64_t>(instance.items[static_cast<std::size_t>(i)].weight) +
                instance.items[static_cast<std::size_t>(j)].weight;
            for (const int item : path_nodes) {
                total += instance.items[static_cast<std::size_t>(item)].weight;
            }
            const int index = static_cast<int>(instance.arcs.size());
            instance.arcs.push_back(Arc{i, j, ceil_div(total, instance.capacity) - 1});
            instance.edge_index[static_cast<std::size_t>(i) * instance.n + j] = index;
            successors[static_cast<std::size_t>(i)].push_back(j);
            changed = true;
        }
    }
    if (changed) {
        instance.rebuild_graph();
    }
    return changed;
}

struct ConflictGraph {
    std::vector<int> offset;
    std::vector<int> neighbors;
    std::vector<std::uint64_t> bits;
    std::size_t blocks = 0;

    [[nodiscard]] bool conflicts(int lhs, int rhs) const noexcept {
        return ((bits[static_cast<std::size_t>(lhs) * blocks +
                      static_cast<std::size_t>(rhs) / 64U] >>
                 (static_cast<unsigned>(rhs) & 63U)) &
                1U) != 0U;
    }
};

[[nodiscard]] ConflictGraph build_conflict_graph(const WorkInstance& instance) {
    ConflictGraph graph;
    graph.blocks = (static_cast<std::size_t>(instance.n) + 63U) / 64U;
    graph.bits.assign(static_cast<std::size_t>(instance.n) * graph.blocks, 0U);
    std::vector<int> degree(static_cast<std::size_t>(instance.n), 0);
    for (const Arc& arc : instance.arcs) {
        if (arc.separation <= 0) {
            continue;
        }
        graph.bits[static_cast<std::size_t>(arc.from) * graph.blocks +
                   static_cast<std::size_t>(arc.to) / 64U] |=
            std::uint64_t{1} << (static_cast<unsigned>(arc.to) & 63U);
        graph.bits[static_cast<std::size_t>(arc.to) * graph.blocks +
                   static_cast<std::size_t>(arc.from) / 64U] |=
            std::uint64_t{1} << (static_cast<unsigned>(arc.from) & 63U);
        ++degree[static_cast<std::size_t>(arc.from)];
        ++degree[static_cast<std::size_t>(arc.to)];
    }
    graph.offset.assign(static_cast<std::size_t>(instance.n + 1), 0);
    for (int i = 0; i < instance.n; ++i) {
        graph.offset[static_cast<std::size_t>(i + 1)] =
            graph.offset[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
    }
    graph.neighbors.resize(static_cast<std::size_t>(graph.offset.back()));
    std::vector<int> cursor = graph.offset;
    for (const Arc& arc : instance.arcs) {
        if (arc.separation <= 0) {
            continue;
        }
        graph.neighbors[static_cast<std::size_t>(cursor[static_cast<std::size_t>(arc.from)]++)] =
            arc.to;
        graph.neighbors[static_cast<std::size_t>(cursor[static_cast<std::size_t>(arc.to)]++)] =
            arc.from;
    }
    return graph;
}

class MaximumFillSolver {
public:
    MaximumFillSolver(const WorkInstance& instance, const ConflictGraph& conflicts)
        : instance_(instance), conflicts_(conflicts), blocked_(static_cast<std::size_t>(instance.n), 0),
          candidate_position_(static_cast<std::size_t>(instance.n), -1) {}

    [[nodiscard]] int solve(int capacity, const std::vector<int>& input_candidates) {
        if (capacity <= 0 || input_candidates.empty()) {
            return 0;
        }
        order_.clear();
        order_.reserve(input_candidates.size());
        for (const int item : input_candidates) {
            if (instance_.items[static_cast<std::size_t>(item)].weight <= capacity) {
                order_.push_back(item);
            }
        }
        if (order_.empty()) {
            return 0;
        }
        for (std::size_t position = 0; position < order_.size(); ++position) {
            candidate_position_[static_cast<std::size_t>(order_[position])] =
                static_cast<int>(position);
        }
        bool active_conflict = false;
        for (const int item : order_) {
            for (int p = conflicts_.offset[static_cast<std::size_t>(item)];
                 p < conflicts_.offset[static_cast<std::size_t>(item + 1)]; ++p) {
                if (candidate_position_[static_cast<std::size_t>(
                        conflicts_.neighbors[static_cast<std::size_t>(p)])] >= 0) {
                    active_conflict = true;
                    break;
                }
            }
            if (active_conflict) {
                break;
            }
        }
        std::stable_sort(order_.begin(), order_.end(), [&](int lhs, int rhs) {
            const int lhs_degree = conflicts_.offset[static_cast<std::size_t>(lhs + 1)] -
                                   conflicts_.offset[static_cast<std::size_t>(lhs)];
            const int rhs_degree = conflicts_.offset[static_cast<std::size_t>(rhs + 1)] -
                                   conflicts_.offset[static_cast<std::size_t>(rhs)];
            if (lhs_degree != rhs_degree) {
                return lhs_degree > rhs_degree;
            }
            return instance_.items[static_cast<std::size_t>(lhs)].weight >
                   instance_.items[static_cast<std::size_t>(rhs)].weight;
        });
        for (int item = 0; item < instance_.n; ++item) {
            candidate_position_[static_cast<std::size_t>(item)] = -1;
        }
        for (std::size_t position = 0; position < order_.size(); ++position) {
            candidate_position_[static_cast<std::size_t>(order_[position])] =
                static_cast<int>(position);
        }

        capacity_ = capacity;
        build_suffix_subset_sum();
        const int relaxed_optimum = suffix_maximum(0, capacity_);
        if (!active_conflict) {
            clear_candidate_positions();
            return relaxed_optimum;
        }

        std::fill(blocked_.begin(), blocked_.end(), 0);
        best_ = greedy_fill();
        if (best_ < relaxed_optimum) {
            search(0, capacity_, 0);
        }
        clear_candidate_positions();
        return best_;
    }

private:
    void clear_candidate_positions() noexcept {
        for (const int item : order_) {
            candidate_position_[static_cast<std::size_t>(item)] = -1;
        }
    }

    void build_suffix_subset_sum() {
        word_count_ = (static_cast<std::size_t>(capacity_) + 64U) / 64U;
        suffix_reachable_.assign((order_.size() + 1U) * word_count_, 0U);
        suffix_reachable_[order_.size() * word_count_] = 1U;
        const unsigned tail_bits = static_cast<unsigned>(capacity_ + 1) & 63U;
        const std::uint64_t tail_mask =
            tail_bits == 0U ? std::numeric_limits<std::uint64_t>::max()
                            : (std::uint64_t{1} << tail_bits) - 1U;
        for (std::size_t reverse = order_.size(); reverse-- > 0U;) {
            std::uint64_t* row = suffix_reachable_.data() + reverse * word_count_;
            const std::uint64_t* next = row + word_count_;
            std::copy(next, next + word_count_, row);
            const unsigned shift = static_cast<unsigned>(
                instance_.items[static_cast<std::size_t>(order_[reverse])].weight);
            const std::size_t word_shift = shift / 64U;
            const unsigned bit_shift = shift & 63U;
            for (std::size_t destination = word_count_; destination-- > word_shift;) {
                std::uint64_t value = next[destination - word_shift] << bit_shift;
                if (bit_shift != 0U && destination > word_shift) {
                    value |= next[destination - word_shift - 1U] >> (64U - bit_shift);
                }
                row[destination] |= value;
            }
            row[word_count_ - 1U] &= tail_mask;
        }
    }

    [[nodiscard]] int suffix_maximum(std::size_t position, int capacity) const noexcept {
        const std::uint64_t* row =
            suffix_reachable_.data() + position * word_count_;
        std::size_t word = static_cast<std::size_t>(capacity) / 64U;
        const unsigned bit = static_cast<unsigned>(capacity) & 63U;
        std::uint64_t value = row[word];
        if (bit != 63U) {
            value &= (std::uint64_t{1} << (bit + 1U)) - 1U;
        }
        while (value == 0U && word > 0U) {
            value = row[--word];
        }
        if (value == 0U) {
            return 0;
        }
        return static_cast<int>(word * 64U + (63U - std::countl_zero(value)));
    }

    [[nodiscard]] int greedy_fill() {
        int value = 0;
        int remaining = capacity_;
        std::vector<int> selected;
        selected.reserve(order_.size());
        for (const int item : order_) {
            const int weight = instance_.items[static_cast<std::size_t>(item)].weight;
            if (weight > remaining || blocked_[static_cast<std::size_t>(item)] != 0) {
                continue;
            }
            selected.push_back(item);
            value += weight;
            remaining -= weight;
            for (int p = conflicts_.offset[static_cast<std::size_t>(item)];
                 p < conflicts_.offset[static_cast<std::size_t>(item + 1)]; ++p) {
                ++blocked_[static_cast<std::size_t>(
                    conflicts_.neighbors[static_cast<std::size_t>(p)])];
            }
        }
        for (const int item : selected) {
            for (int p = conflicts_.offset[static_cast<std::size_t>(item)];
                 p < conflicts_.offset[static_cast<std::size_t>(item + 1)]; ++p) {
                --blocked_[static_cast<std::size_t>(
                    conflicts_.neighbors[static_cast<std::size_t>(p)])];
            }
        }
        return value;
    }

    void search(std::size_t position, int remaining, int value) {
        if (best_ == capacity_ || position == order_.size()) {
            best_ = std::max(best_, value);
            return;
        }
        const int upper_bound = value + suffix_maximum(position, remaining);
        if (upper_bound <= best_) {
            return;
        }
        const int item = order_[position];
        const int weight = instance_.items[static_cast<std::size_t>(item)].weight;
        if (blocked_[static_cast<std::size_t>(item)] == 0 && weight <= remaining) {
            for (int p = conflicts_.offset[static_cast<std::size_t>(item)];
                 p < conflicts_.offset[static_cast<std::size_t>(item + 1)]; ++p) {
                ++blocked_[static_cast<std::size_t>(
                    conflicts_.neighbors[static_cast<std::size_t>(p)])];
            }
            search(position + 1U, remaining - weight, value + weight);
            for (int p = conflicts_.offset[static_cast<std::size_t>(item)];
                 p < conflicts_.offset[static_cast<std::size_t>(item + 1)]; ++p) {
                --blocked_[static_cast<std::size_t>(
                    conflicts_.neighbors[static_cast<std::size_t>(p)])];
            }
            if (best_ == capacity_) {
                return;
            }
        }
        search(position + 1U, remaining, value);
    }

    const WorkInstance& instance_;
    const ConflictGraph& conflicts_;
    std::vector<int> blocked_;
    std::vector<int> candidate_position_;
    std::vector<int> order_;
    std::vector<std::uint64_t> suffix_reachable_;
    std::size_t word_count_ = 0;
    int capacity_ = 0;
    int best_ = 0;
};

[[nodiscard]] bool shrink_bins_and_enlarge_items(WorkInstance& instance) {
    const ConflictGraph conflicts = build_conflict_graph(instance);
    MaximumFillSolver solver(instance, conflicts);
    std::vector<int> candidates(static_cast<std::size_t>(instance.n));
    std::iota(candidates.begin(), candidates.end(), 0);
    bool changed = false;
    const int tightened_capacity = solver.solve(instance.capacity, candidates);
    if (tightened_capacity <= 0) {
        throw std::logic_error("bin shrinking produced a nonpositive capacity");
    }
    if (tightened_capacity < instance.capacity) {
        instance.capacity = tightened_capacity;
        changed = true;
    }

    candidates.clear();
    candidates.reserve(static_cast<std::size_t>(instance.n));
    for (int item = 0; item < instance.n; ++item) {
        candidates.clear();
        for (int other = 0; other < instance.n; ++other) {
            if (other != item && !conflicts.conflicts(item, other)) {
                candidates.push_back(other);
            }
        }
        const int old_weight = instance.items[static_cast<std::size_t>(item)].weight;
        const int companion_load =
            candidates.empty() ? 0 : solver.solve(instance.capacity - old_weight, candidates);
        const int lifted_weight = instance.capacity - companion_load;
        if (lifted_weight > old_weight) {
            instance.items[static_cast<std::size_t>(item)].weight = lifted_weight;
            changed = true;
        }
    }
    instance.total_weight = 0;
    for (const WorkItem& item : instance.items) {
        instance.total_weight += item.weight;
    }
    return changed;
}

void flip_precedence_graph(WorkInstance& instance) {
    instance.reversed = !instance.reversed;
    for (Arc& arc : instance.arcs) {
        std::swap(arc.from, arc.to);
    }
    instance.rebuild_graph();
}

[[nodiscard]] bool reverse_precedence_graph(WorkInstance& instance) {
    int forward_sum = 0;
    int reverse_sum = 0;
    for (int item = 0; item < instance.n; ++item) {
        for (int level = 1; level <= 5; ++level) {
            forward_sum += instance.front[static_cast<std::size_t>(item)] <= level;
            reverse_sum += instance.back[static_cast<std::size_t>(item)] <= level;
        }
    }
    if (forward_sum <= reverse_sum) {
        return false;
    }
    flip_precedence_graph(instance);
    return true;
}

struct UbMinusOneWindows {
    std::vector<int> earliest;
    std::vector<int> latest;
    bool infeasible = false;
    bool prefer_reverse = false;
};

[[nodiscard]] UbMinusOneWindows compute_ub_minus_one_windows(
    const WorkInstance& instance, int target_bins) {
    UbMinusOneWindows result;
    result.earliest.assign(static_cast<std::size_t>(instance.n), 0);
    result.latest.assign(static_cast<std::size_t>(instance.n), target_bins - 1);
    if (target_bins <= 0) {
        result.infeasible = true;
        return result;
    }

    const std::size_t blocks =
        (static_cast<std::size_t>(instance.n) + 63U) / 64U;
    const std::vector<std::uint64_t> reachability =
        compute_reachability(instance);
    const bool strict_precedence = instance.problem_type == "BPP-P";
    for (int item = 0; item < instance.n; ++item) {
        std::int64_t predecessor_weight = 0;
        std::int64_t successor_weight = 0;
        const std::uint64_t* successors = reachability.data() +
            static_cast<std::size_t>(item) * blocks;
        for (std::size_t block = 0; block < blocks; ++block) {
            std::uint64_t bits = successors[block];
            while (bits != 0U) {
                const unsigned offset = std::countr_zero(bits);
                const int successor =
                    static_cast<int>(block * 64U + offset);
                successor_weight += instance.items[
                    static_cast<std::size_t>(successor)].weight;
                bits &= bits - 1U;
            }
        }
        for (int predecessor = 0; predecessor < instance.n; ++predecessor) {
            const std::uint64_t* predecessor_successors =
                reachability.data() +
                static_cast<std::size_t>(predecessor) * blocks;
            if (((predecessor_successors[
                      static_cast<std::size_t>(item) >> 6U] >>
                  (static_cast<unsigned>(item) & 63U)) &
                 1U) != 0U) {
                predecessor_weight += instance.items[
                    static_cast<std::size_t>(predecessor)].weight;
            }
        }

        const std::int64_t own_weight =
            instance.items[static_cast<std::size_t>(item)].weight;
        const int weight_release = strict_precedence
            ? ceil_div(predecessor_weight, instance.capacity)
            : std::max(0, ceil_div(predecessor_weight + own_weight,
                                   instance.capacity) - 1);
        const int weight_tail = strict_precedence
            ? ceil_div(successor_weight, instance.capacity)
            : std::max(0, ceil_div(successor_weight + own_weight,
                                   instance.capacity) - 1);
        const int earliest = std::max(
            instance.front[static_cast<std::size_t>(item)], weight_release);
        const int tail = std::max(
            instance.back[static_cast<std::size_t>(item)], weight_tail);
        result.earliest[static_cast<std::size_t>(item)] = earliest;
        result.latest[static_cast<std::size_t>(item)] =
            target_bins - 1 - tail;
        result.infeasible = result.infeasible ||
            earliest > result.latest[static_cast<std::size_t>(item)];
    }
    if (result.infeasible) {
        return result;
    }

    const int inspected_bins = std::min(5, target_bins);
    std::uint64_t forward_product = 1U;
    std::uint64_t reverse_product = 1U;
    const auto multiply_saturated = [](std::uint64_t lhs,
                                       std::uint64_t rhs) noexcept {
        if (lhs == 0U || rhs == 0U) {
            return std::uint64_t{0};
        }
        if (lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
            return std::numeric_limits<std::uint64_t>::max();
        }
        return lhs * rhs;
    };
    for (int offset = 0; offset < inspected_bins; ++offset) {
        const int forward_bin = offset;
        const int reverse_bin = target_bins - 1 - offset;
        std::uint64_t forward_candidates = 0;
        std::uint64_t reverse_candidates = 0;
        for (int item = 0; item < instance.n; ++item) {
            const int earliest =
                result.earliest[static_cast<std::size_t>(item)];
            const int latest = result.latest[static_cast<std::size_t>(item)];
            forward_candidates += earliest <= forward_bin;
            reverse_candidates += latest >= reverse_bin;
        }
        forward_product = multiply_saturated(forward_product,
                                             forward_candidates);
        reverse_product = multiply_saturated(reverse_product,
                                             reverse_candidates);
    }
    result.prefer_reverse = reverse_product < forward_product;
    return result;
}

void classify_problem_type(WorkInstance& instance) {
    if (instance.arcs.empty()) {
        instance.problem_type = "BPP";
        return;
    }
    std::int64_t sum = 0;
    int maximum = 0;
    for (const Arc& arc : instance.arcs) {
        sum += arc.separation;
        maximum = std::max(maximum, arc.separation);
    }
    if (sum == 0) {
        instance.problem_type = "SALBP-I";
    } else if (maximum == 1 && sum == static_cast<std::int64_t>(instance.arcs.size())) {
        instance.problem_type = "BPP-P";
    } else {
        instance.problem_type = "BPP-GP";
    }
}

void preprocess(WorkInstance& instance) {
    bool changed = true;
    int rounds = 0;
    while (changed) {
        const bool graph_changed = enhance_precedence_graph(instance);
        const bool weights_changed = shrink_bins_and_enlarge_items(instance);
        const bool reversed = reverse_precedence_graph(instance);
        const bool reordered = instance.reorder_items();
        changed = graph_changed || weights_changed || reversed || reordered;
        if (++rounds > 4 * instance.n + 32) {
            throw std::logic_error("initial-bound preprocessing did not reach a fixed point");
        }
    }
    instance.recompute_position_bounds();
    classify_problem_type(instance);
}

struct CertifiedBounds {
    int capacity = 0;
    int precedence_path = 0;
    int window_dff = 0;
};

[[nodiscard]] CertifiedBounds quick_bounds(const WorkInstance& instance) {
    CertifiedBounds bounds;
    bounds.capacity = ceil_div(instance.total_weight, instance.capacity);
    bounds.precedence_path = 1;
    for (const int value : instance.front) {
        bounds.precedence_path =
            std::max(bounds.precedence_path, value + 1);
    }
    return bounds;
}

struct HeuristicState {
    explicit HeuristicState(int n) : bin(static_cast<std::size_t>(n), -1) {}
    std::vector<int> bin;
    std::vector<int> remaining_capacity;
    int bin_count = 0;
};

[[nodiscard]] bool better_state(const HeuristicState& lhs,
                                const HeuristicState& rhs) noexcept {
    return rhs.bin_count == 0 || lhs.bin_count < rhs.bin_count;
}

class FitBasedHeuristic {
public:
    explicit FitBasedHeuristic(const WorkInstance& instance) : instance_(instance) {}

    [[nodiscard]] HeuristicState solve(int lower_bound,
                                       JavaRandom& random,
                                       int random_trials) const {
        std::vector<int> sequence(static_cast<std::size_t>(instance_.n));
        std::iota(sequence.begin(), sequence.end(), 0);
        HeuristicState best = evaluate(sequence);
        if (best.bin_count == lower_bound) {
            return best;
        }
        std::stable_sort(sequence.begin(), sequence.end(), [&](int lhs, int rhs) {
            return instance_.items[static_cast<std::size_t>(lhs)].weight >
                   instance_.items[static_cast<std::size_t>(rhs)].weight;
        });
        HeuristicState candidate = evaluate(sequence);
        if (better_state(candidate, best)) {
            best = std::move(candidate);
        }
        if (best.bin_count == lower_bound) {
            return best;
        }
        for (int trial = 0; trial < random_trials; ++trial) {
            random.shuffle(sequence);
            candidate = evaluate(sequence);
            if (better_state(candidate, best)) {
                best = std::move(candidate);
            }
            if (best.bin_count == lower_bound) {
                break;
            }
        }
        return best;
    }

private:
    [[nodiscard]] std::vector<int> initial_predecessor_count() const {
        std::vector<int> result(static_cast<std::size_t>(instance_.n), 0);
        for (int item = 0; item < instance_.n; ++item) {
            result[static_cast<std::size_t>(item)] =
                instance_.pred_offset[static_cast<std::size_t>(item + 1)] -
                instance_.pred_offset[static_cast<std::size_t>(item)];
        }
        return result;
    }

    [[nodiscard]] int first_available(const std::vector<int>& sequence,
                                      const HeuristicState& state,
                                      const std::vector<int>& predecessor_count) const {
        for (const int item : sequence) {
            if (state.bin[static_cast<std::size_t>(item)] < 0 &&
                predecessor_count[static_cast<std::size_t>(item)] == 0) {
                return item;
            }
        }
        return -1;
    }

    [[nodiscard]] bool can_pack(const HeuristicState& state,
                                int item,
                                int bin,
                                int remaining_after) const {
        if (remaining_after < 0) {
            return false;
        }
        for (int p = instance_.pred_offset[static_cast<std::size_t>(item)];
             p < instance_.pred_offset[static_cast<std::size_t>(item + 1)]; ++p) {
            const int predecessor = instance_.pred_from[static_cast<std::size_t>(p)];
            const int predecessor_bin = state.bin[static_cast<std::size_t>(predecessor)];
            const int separation =
                instance_.arcs[static_cast<std::size_t>(
                    instance_.pred_arc[static_cast<std::size_t>(p)])]
                    .separation;
            if (predecessor_bin < 0 || bin - predecessor_bin < separation) {
                return false;
            }
        }
        for (int p = instance_.succ_offset[static_cast<std::size_t>(item)];
             p < instance_.succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
            const int successor = instance_.succ_to[static_cast<std::size_t>(p)];
            const int successor_bin = state.bin[static_cast<std::size_t>(successor)];
            const int separation =
                instance_.arcs[static_cast<std::size_t>(
                    instance_.succ_arc[static_cast<std::size_t>(p)])]
                    .separation;
            if (successor_bin >= 0 && successor_bin - bin < separation) {
                return false;
            }
        }
        return true;
    }

    void pack(HeuristicState& state,
              std::vector<int>& predecessor_count,
              int item,
              int bin) const {
        state.bin[static_cast<std::size_t>(item)] = bin;
        for (int p = instance_.succ_offset[static_cast<std::size_t>(item)];
             p < instance_.succ_offset[static_cast<std::size_t>(item + 1)]; ++p) {
            --predecessor_count[static_cast<std::size_t>(
                instance_.succ_to[static_cast<std::size_t>(p)])];
        }
    }

    void open_until_packable(HeuristicState& state,
                             std::vector<int>& predecessor_count,
                             int item) const {
        while (true) {
            const int bin = state.bin_count++;
            state.remaining_capacity.push_back(instance_.capacity);
            if (can_pack(state, item, bin,
                         instance_.capacity -
                             instance_.items[static_cast<std::size_t>(item)].weight)) {
                pack(state, predecessor_count, item, bin);
                state.remaining_capacity[static_cast<std::size_t>(bin)] -=
                    instance_.items[static_cast<std::size_t>(item)].weight;
                return;
            }
        }
    }

    [[nodiscard]] HeuristicState first_fit(const std::vector<int>& sequence) const {
        HeuristicState state(instance_.n);
        state.bin_count = 1;
        state.remaining_capacity.push_back(instance_.capacity);
        std::vector<int> predecessor_count = initial_predecessor_count();
        for (int packed = 0; packed < instance_.n; ++packed) {
            const int item = first_available(sequence, state, predecessor_count);
            if (item < 0) {
                throw std::logic_error("fit heuristic found no precedence-ready item");
            }
            bool placed = false;
            for (int bin = 0; bin < state.bin_count; ++bin) {
                const int after = state.remaining_capacity[static_cast<std::size_t>(bin)] -
                                  instance_.items[static_cast<std::size_t>(item)].weight;
                if (can_pack(state, item, bin, after)) {
                    pack(state, predecessor_count, item, bin);
                    state.remaining_capacity[static_cast<std::size_t>(bin)] = after;
                    placed = true;
                    break;
                }
            }
            if (!placed) {
                open_until_packable(state, predecessor_count, item);
            }
        }
        return state;
    }

    [[nodiscard]] int better_fit_replacement(HeuristicState& state,
                                             std::vector<int>& predecessor_count,
                                             int item) const {
        while (true) {
            int best_bin = -1;
            int best_replaced = -1;
            int minimum_remaining = std::numeric_limits<int>::max();
            const int item_weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            for (int replaced = 0; replaced < instance_.n; ++replaced) {
                const int bin = state.bin[static_cast<std::size_t>(replaced)];
                if (bin < 0) {
                    continue;
                }
                const int new_remaining =
                    state.remaining_capacity[static_cast<std::size_t>(bin)] + item_weight -
                    instance_.items[static_cast<std::size_t>(replaced)].weight;
                if (new_remaining < 0 ||
                    new_remaining >= state.remaining_capacity[static_cast<std::size_t>(bin)] ||
                    new_remaining >= minimum_remaining) {
                    continue;
                }
                bool feasible = instance_.edge_position(replaced, item) < 0;
                for (int packed = replaced + 1; feasible && packed < instance_.n; ++packed) {
                    if (state.bin[static_cast<std::size_t>(packed)] >= 0 &&
                        instance_.edge_position(replaced, packed) >= 0) {
                        feasible = false;
                    }
                }
                for (int p = instance_.pred_offset[static_cast<std::size_t>(item)];
                     feasible && p < instance_.pred_offset[static_cast<std::size_t>(item + 1)];
                     ++p) {
                    const int predecessor = instance_.pred_from[static_cast<std::size_t>(p)];
                    const int predecessor_bin = state.bin[static_cast<std::size_t>(predecessor)];
                    const int separation =
                        instance_.arcs[static_cast<std::size_t>(
                            instance_.pred_arc[static_cast<std::size_t>(p)])]
                            .separation;
                    if (predecessor == replaced || predecessor_bin < 0 ||
                        bin - predecessor_bin < separation) {
                        feasible = false;
                    }
                }
                for (int p = instance_.succ_offset[static_cast<std::size_t>(item)];
                     feasible && p < instance_.succ_offset[static_cast<std::size_t>(item + 1)];
                     ++p) {
                    const int successor = instance_.succ_to[static_cast<std::size_t>(p)];
                    const int successor_bin = state.bin[static_cast<std::size_t>(successor)];
                    const int separation =
                        instance_.arcs[static_cast<std::size_t>(
                            instance_.succ_arc[static_cast<std::size_t>(p)])]
                            .separation;
                    if (successor != replaced && successor_bin >= 0 &&
                        successor_bin - bin < separation) {
                        feasible = false;
                    }
                }
                if (feasible) {
                    best_bin = bin;
                    best_replaced = replaced;
                    minimum_remaining = new_remaining;
                }
            }
            if (best_bin < 0) {
                return item;
            }
            state.bin[static_cast<std::size_t>(best_replaced)] = -1;
            for (int p = instance_.succ_offset[static_cast<std::size_t>(best_replaced)];
                 p < instance_.succ_offset[static_cast<std::size_t>(best_replaced + 1)]; ++p) {
                ++predecessor_count[static_cast<std::size_t>(
                    instance_.succ_to[static_cast<std::size_t>(p)])];
            }
            pack(state, predecessor_count, item, best_bin);
            state.remaining_capacity[static_cast<std::size_t>(best_bin)] = minimum_remaining;
            item = best_replaced;
        }
    }

    [[nodiscard]] HeuristicState better_fit(const std::vector<int>& sequence) const {
        HeuristicState state(instance_.n);
        state.bin_count = 1;
        state.remaining_capacity.push_back(instance_.capacity);
        std::vector<int> predecessor_count = initial_predecessor_count();
        for (int packed = 0; packed < instance_.n; ++packed) {
            int item = first_available(sequence, state, predecessor_count);
            if (item < 0) {
                throw std::logic_error("better-fit heuristic found no ready item");
            }
            item = better_fit_replacement(state, predecessor_count, item);
            int best_bin = -1;
            int best_remaining = -1;
            for (int bin = 0; bin < state.bin_count; ++bin) {
                const int after = state.remaining_capacity[static_cast<std::size_t>(bin)] -
                                  instance_.items[static_cast<std::size_t>(item)].weight;
                if (can_pack(state, item, bin, after) &&
                    (best_bin < 0 || after < best_remaining)) {
                    best_bin = bin;
                    best_remaining = after;
                }
            }
            if (best_bin < 0) {
                open_until_packable(state, predecessor_count, item);
            } else {
                pack(state, predecessor_count, item, best_bin);
                state.remaining_capacity[static_cast<std::size_t>(best_bin)] =
                    best_remaining;
            }
        }
        return state;
    }

    [[nodiscard]] HeuristicState evaluate(const std::vector<int>& sequence) const {
        HeuristicState best = better_fit(sequence);
        HeuristicState candidate = first_fit(sequence);
        if (better_state(candidate, best)) {
            best = std::move(candidate);
        }
        return best;
    }

    const WorkInstance& instance_;
};

template <typename WeightAt>
[[nodiscard]] int exact_dff_lower_bound(std::size_t item_count,
                                        int capacity,
                                        WeightAt&& weight_at) {
    if (capacity <= 0) {
        throw std::invalid_argument("exact DFF requires a positive capacity");
    }

    int lower_bound = 1;

    for (int parameter = 1; parameter <= 100; ++parameter) {
        const std::int64_t denominator =
            static_cast<std::int64_t>(capacity) * parameter;
        exact_arithmetic::NonnegativeRatioSum ratio(denominator);
        for (std::size_t item = 0; item < item_count; ++item) {
            const std::int64_t weight = weight_at(item);
            const std::int64_t scaled =
                static_cast<std::int64_t>(parameter + 1) * weight;
            if (scaled % capacity == 0) {
                ratio.add(weight * parameter);
            } else {
                ratio.add((scaled / capacity) * capacity);
            }
        }
        lower_bound = std::max(lower_bound, ratio.ceil_to_int());
    }

    exact_arithmetic::NonnegativeRatioSum half_ratio(capacity);
    for (std::size_t item = 0; item < item_count; ++item) {
        const std::int64_t weight = weight_at(item);
        if (2 * weight > capacity) {
            half_ratio.add(capacity);
        } else if (2 * weight == capacity) {
            half_ratio.add(weight);
        }
    }
    lower_bound = std::max(lower_bound, half_ratio.ceil_to_int());

    for (std::size_t parameter_item = 0;
         parameter_item < item_count;
         ++parameter_item) {
        const std::int64_t parameter_weight = weight_at(parameter_item);
        if (2 * parameter_weight >= capacity) {
            continue;
        }
        exact_arithmetic::NonnegativeRatioSum threshold_ratio(capacity);
        for (std::size_t item = 0; item < item_count; ++item) {
            const std::int64_t weight = weight_at(item);
            if (weight > capacity - parameter_weight) {
                threshold_ratio.add(capacity);
            } else if (weight >= parameter_weight) {
                threshold_ratio.add(weight);
            }
        }
        lower_bound = std::max(lower_bound, threshold_ratio.ceil_to_int());
    }

    for (int parameter_numerator = 1;
         parameter_numerator <= 500;
         ++parameter_numerator) {
        const std::int64_t transformed_capacity =
            1000 / parameter_numerator;
        const std::int64_t ratio_denominator =
            static_cast<std::int64_t>(capacity) * parameter_numerator;
        exact_arithmetic::NonnegativeRatioSum transformed_ratio(
            transformed_capacity);
        for (std::size_t item = 0; item < item_count; ++item) {
            const std::int64_t weight = weight_at(item);
            if (2 * weight > capacity) {
                const std::int64_t complement_units =
                    ((static_cast<std::int64_t>(capacity) - weight) * 1000) /
                    ratio_denominator;
                transformed_ratio.add(
                    transformed_capacity - complement_units);
            } else if (weight * 1000 >= ratio_denominator) {
                transformed_ratio.add(
                    (weight * 1000) / ratio_denominator);
            }
        }
        lower_bound = std::max(
            lower_bound, transformed_ratio.ceil_to_int());
    }
    return lower_bound;
}

[[nodiscard]] int dff_lower_bound(const std::vector<int>& selected,
                                  const WorkInstance& instance) {
    return exact_dff_lower_bound(
        selected.size(), instance.capacity, [&](std::size_t position) {
            const int item = selected[position];
            return instance.items[static_cast<std::size_t>(item)].weight;
        });
}

[[nodiscard]] int compute_window_dff_lower_bound(
    const WorkInstance& instance) {
    int longest = 0;
    for (const int value : instance.front) {
        longest = std::max(longest, value);
    }
    int lower_bound = 1;
    std::vector<int> front_items;
    std::vector<int> selected;
    front_items.reserve(static_cast<std::size_t>(instance.n));
    selected.reserve(static_cast<std::size_t>(instance.n));
    for (int g = 1; g <= longest; ++g) {
        front_items.clear();
        for (int item = 0; item < instance.n; ++item) {
            if (instance.front[static_cast<std::size_t>(item)] >= g) {
                front_items.push_back(item);
            }
        }
        for (int h = 1; h <= longest - g; ++h) {
            selected.clear();
            for (const int item : front_items) {
                if (instance.back[static_cast<std::size_t>(item)] >= h) {
                    selected.push_back(item);
                }
            }
            const int dff_bound = dff_lower_bound(selected, instance);
            const int candidate = g + h + dff_bound;
            lower_bound = std::max(lower_bound, candidate);
        }
    }
    return lower_bound;
}

class ModifiedHoffmannHeuristic {
public:
    explicit ModifiedHoffmannHeuristic(const WorkInstance& instance,
                                       bool paper_bbr_mode = false,
                                       bool bbr12_dynamic_mode = false,
                                       bool bbr12_portfolio_mode = false,
                                       int bbr12_portfolio_max_items = 200,
                                       int bbr12_full_load_limit = 5'000)
        : instance_(instance),
          paper_bbr_mode_(paper_bbr_mode),
          bbr12_dynamic_mode_(bbr12_dynamic_mode),
          bbr12_portfolio_mode_(bbr12_portfolio_mode),
          bbr12_portfolio_max_items_(bbr12_portfolio_max_items),
          bbr12_full_load_limit_(bbr12_full_load_limit),
          bins_(static_cast<std::size_t>(instance.n), -1) {
        if (bbr12_full_load_limit_ <= 0) {
            throw std::invalid_argument(
                "BBR12 Hoffmann full-load limit must be positive");
        }
        if (bbr12_portfolio_max_items_ <= 0) {
            throw std::invalid_argument(
                "BBR12 Hoffmann portfolio item limit must be positive");
        }
        const std::vector<std::uint64_t> reachability = compute_reachability(instance_);
        const std::size_t blocks =
            (static_cast<std::size_t>(instance_.n) + 63U) / 64U;
        successor_weight_.assign(static_cast<std::size_t>(instance_.n), 0);
        successor_count_.assign(static_cast<std::size_t>(instance_.n), 0);
        direct_successor_count_.assign(
            static_cast<std::size_t>(instance_.n), 0);
        for (int item = 0; item < instance_.n; ++item) {
            direct_successor_count_[static_cast<std::size_t>(item)] =
                instance_.succ_offset[static_cast<std::size_t>(item + 1)] -
                instance_.succ_offset[static_cast<std::size_t>(item)];
            const std::uint64_t* row =
                reachability.data() + static_cast<std::size_t>(item) * blocks;
            for (std::size_t block = 0; block < blocks; ++block) {
                std::uint64_t bits = row[block];
                while (bits != 0U) {
                    const unsigned offset = std::countr_zero(bits);
                    const int successor = static_cast<int>(block * 64U + offset);
                    successor_weight_[static_cast<std::size_t>(item)] +=
                        instance_.items[static_cast<std::size_t>(successor)].weight;
                    ++successor_count_[static_cast<std::size_t>(item)];
                    bits &= bits - 1U;
                }
            }
        }
        if (bbr12_dynamic_mode_) {
            initialize_bbr12_workspace();
        }
    }

    [[nodiscard]] HeuristicState solve(int lower_bound, const Deadline& deadline) {
        if (bbr12_dynamic_mode_) {
            return solve_bbr12(lower_bound, deadline);
        }
        static constexpr std::array<double, 5> java_coefficients{
            0.0, 0.005, 0.01, 0.15, 0.02};
        static constexpr std::array<double, 5> paper_coefficients{
            0.0, 0.005, 0.01, 0.015, 0.02};
        static constexpr std::array<double, 4> gammas{0.0, 0.01, 0.02, 0.03};
        const auto& alphas = paper_bbr_mode_ ? paper_coefficients
                                             : java_coefficients;
        const auto& betas = alphas;
        HeuristicState best(instance_.n);
        for (const double alpha : alphas) {
            for (const double beta : betas) {
                for (const double gamma : gammas) {
                    HeuristicState candidate = construct(alpha, beta, gamma);
                    if (better_state(candidate, best)) {
                        best = std::move(candidate);
                    }
                    if (best.bin_count == lower_bound || deadline.expired()) {
                        return best;
                    }
                }
            }
        }
        return best;
    }

private:
    void initialize_bbr12_workspace() {
        bbr12_item_order_.resize(static_cast<std::size_t>(instance_.n));
        std::iota(bbr12_item_order_.begin(), bbr12_item_order_.end(), 0);
        const auto paper_index = [&](int item) {
            const int original = instance_.items[
                static_cast<std::size_t>(item)].original_index;
            return instance_.reversed ? -original : original;
        };
        std::stable_sort(
            bbr12_item_order_.begin(), bbr12_item_order_.end(),
            [&](int lhs, int rhs) {
                return std::pair{paper_index(lhs), lhs} <
                       std::pair{paper_index(rhs), rhs};
            });
        bbr12_rank_.resize(static_cast<std::size_t>(instance_.n));
        for (int rank = 0; rank < instance_.n; ++rank) {
            bbr12_rank_[static_cast<std::size_t>(bbr12_item_order_[
                static_cast<std::size_t>(rank)])] = rank;
        }

        bbr12_successor_offset_.assign(
            static_cast<std::size_t>(instance_.n + 1), 0);
        for (int item = 0; item < instance_.n; ++item) {
            bbr12_successor_offset_[static_cast<std::size_t>(item + 1)] =
                bbr12_successor_offset_[static_cast<std::size_t>(item)] +
                direct_successor_count_[static_cast<std::size_t>(item)];
        }
        bbr12_successors_.resize(
            static_cast<std::size_t>(bbr12_successor_offset_.back()));
        for (int item = 0; item < instance_.n; ++item) {
            const int begin =
                bbr12_successor_offset_[static_cast<std::size_t>(item)];
            const int end =
                bbr12_successor_offset_[static_cast<std::size_t>(item + 1)];
            int cursor = begin;
            for (int position =
                     instance_.succ_offset[static_cast<std::size_t>(item)];
                 position < instance_.succ_offset[
                                static_cast<std::size_t>(item + 1)];
                 ++position) {
                bbr12_successors_[static_cast<std::size_t>(cursor++)] =
                    instance_.succ_to[static_cast<std::size_t>(position)];
            }
            std::sort(
                bbr12_successors_.begin() + begin,
                bbr12_successors_.begin() + end,
                [&](int lhs, int rhs) {
                    return bbr12_rank_[static_cast<std::size_t>(lhs)] <
                           bbr12_rank_[static_cast<std::size_t>(rhs)];
                });
        }
        bbr12_degrees_.resize(static_cast<std::size_t>(instance_.n));
        bbr12_eligible_.resize(static_cast<std::size_t>(instance_.n));
        bbr12_tasks_.resize(static_cast<std::size_t>(instance_.n));
        bbr12_best_tasks_.resize(static_cast<std::size_t>(instance_.n));
        bbr12_item_score_.resize(static_cast<std::size_t>(instance_.n));
    }

    [[nodiscard]] HeuristicState solve_bbr12(
        int lower_bound, const Deadline& deadline) {
        if (!bbr12_portfolio_mode_) {
            return solve_bbr12_at_current_limit(lower_bound, deadline);
        }

        static constexpr std::array<int, 4> portfolio_limits{
            50, 250, 500, 1'000};
        const int configured_limit = bbr12_full_load_limit_;
        if (instance_.n > bbr12_portfolio_max_items_) {
            bbr12_full_load_limit_ = std::min(configured_limit, 50);
            HeuristicState best =
                solve_bbr12_at_current_limit(lower_bound, deadline);
            bbr12_full_load_limit_ = configured_limit;
            return best;
        }
        HeuristicState best(instance_.n);
        int previous_limit = 0;
        const auto run_limit = [&](int limit) {
            bbr12_full_load_limit_ = limit;
            HeuristicState candidate =
                solve_bbr12_at_current_limit(lower_bound, deadline);
            if (better_state(candidate, best)) {
                best = std::move(candidate);
            }
            previous_limit = limit;
        };
        for (const int limit : portfolio_limits) {
            if (limit > configured_limit) {
                break;
            }
            run_limit(limit);
            if (best.bin_count == lower_bound || deadline.expired()) {
                bbr12_full_load_limit_ = configured_limit;
                return best;
            }
        }
        if (previous_limit != configured_limit && !deadline.expired()) {
            run_limit(configured_limit);
        }
        bbr12_full_load_limit_ = configured_limit;
        return best;
    }

    [[nodiscard]] HeuristicState solve_bbr12_at_current_limit(
        int lower_bound, const Deadline& deadline) {
        static constexpr std::array<int, 5> coefficients_milli{
            0, 5, 10, 15, 20};
        static constexpr std::array<int, 4> gammas_milli{0, 10, 20, 30};
        HeuristicState best(instance_.n);
        for (const int alpha : coefficients_milli) {
            for (const int beta : coefficients_milli) {
                for (const int gamma : gammas_milli) {
                    HeuristicState candidate =
                        construct_bbr12(alpha, beta, gamma);
                    if (better_state(candidate, best)) {
                        best = std::move(candidate);
                    }
                    if (best.bin_count == lower_bound || deadline.expired()) {
                        return best;
                    }
                }
            }
        }
        return best;
    }

    [[nodiscard]] bool bbr12_packable(int item,
                                      int bin,
                                      int remaining) const noexcept {
        if (remaining <
            instance_.items[static_cast<std::size_t>(item)].weight) {
            return false;
        }
        for (int position =
                 instance_.pred_offset[static_cast<std::size_t>(item)];
             position < instance_.pred_offset[
                            static_cast<std::size_t>(item + 1)];
             ++position) {
            const int predecessor =
                instance_.pred_from[static_cast<std::size_t>(position)];
            const int predecessor_bin =
                bins_[static_cast<std::size_t>(predecessor)];
            const int separation = instance_.arcs[static_cast<std::size_t>(
                instance_.pred_arc[static_cast<std::size_t>(position)])]
                                       .separation;
            if (predecessor_bin < 0 ||
                bin - predecessor_bin < separation) {
                return false;
            }
        }
        return true;
    }

    void enumerate_bbr12(int depth,
                         int remaining,
                         int start,
                         int eligible_count,
                         int bin,
                         std::int64_t score) {
        if (bbr12_full_load_count_ >= bbr12_full_load_limit_) {
            return;
        }
        bool full_load = true;
        for (int position = start; position < eligible_count; ++position) {
            const int item =
                bbr12_eligible_[static_cast<std::size_t>(position)];
            if (bbr12_degrees_[static_cast<std::size_t>(item)] != 0 ||
                !bbr12_packable(item, bin, remaining)) {
                continue;
            }
            full_load = false;
            bbr12_tasks_[static_cast<std::size_t>(depth)] = item;
            bins_[static_cast<std::size_t>(item)] = bin;
            bbr12_degrees_[static_cast<std::size_t>(item)] = -1;
            int child_eligible_count = eligible_count;
            const int successor_begin = bbr12_successor_offset_[
                static_cast<std::size_t>(item)];
            const int successor_end = bbr12_successor_offset_[
                static_cast<std::size_t>(item + 1)];
            for (int successor_position = successor_begin;
                 successor_position < successor_end; ++successor_position) {
                const int successor = bbr12_successors_[
                    static_cast<std::size_t>(successor_position)];
                int& degree =
                    bbr12_degrees_[static_cast<std::size_t>(successor)];
                --degree;
                if (degree == 0) {
                    bbr12_eligible_[static_cast<std::size_t>(
                        child_eligible_count++)] = successor;
                }
            }

            const std::int64_t child_score = score +
                bbr12_item_score_[static_cast<std::size_t>(item)];
            if (bbr12_best_task_count_ == 0 ||
                child_score > bbr12_best_score_) {
                bbr12_best_score_ = child_score;
                bbr12_best_task_count_ = depth + 1;
                std::copy_n(bbr12_tasks_.begin(),
                            bbr12_best_task_count_,
                            bbr12_best_tasks_.begin());
            }
            enumerate_bbr12(
                depth + 1,
                remaining -
                    instance_.items[static_cast<std::size_t>(item)].weight,
                position + 1, child_eligible_count, bin, child_score);

            for (int successor_position = successor_begin;
                 successor_position < successor_end; ++successor_position) {
                const int successor = bbr12_successors_[
                    static_cast<std::size_t>(successor_position)];
                ++bbr12_degrees_[static_cast<std::size_t>(successor)];
            }
            bbr12_degrees_[static_cast<std::size_t>(item)] = 0;
            bins_[static_cast<std::size_t>(item)] = -1;
            if (bbr12_full_load_count_ >= bbr12_full_load_limit_) {
                return;
            }
        }
        bbr12_full_load_count_ += full_load ? 1 : 0;
    }

    [[nodiscard]] HeuristicState construct_bbr12(
        int alpha_milli, int beta_milli, int gamma_milli) {
        std::fill(bins_.begin(), bins_.end(), -1);
        for (int item = 0; item < instance_.n; ++item) {
            bbr12_degrees_[static_cast<std::size_t>(item)] =
                instance_.pred_offset[static_cast<std::size_t>(item + 1)] -
                instance_.pred_offset[static_cast<std::size_t>(item)];
            const std::int64_t weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            const std::int64_t positional_weight =
                weight + successor_weight_[static_cast<std::size_t>(item)];
            bbr12_item_score_[static_cast<std::size_t>(item)] =
                1000 * weight + alpha_milli * positional_weight +
                beta_milli * direct_successor_count_[
                    static_cast<std::size_t>(item)] -
                gamma_milli;
        }

        int bin_count = 0;
        int packed = 0;
        while (packed < instance_.n) {
            int eligible_count = 0;
            for (const int item : bbr12_item_order_) {
                if (bbr12_degrees_[static_cast<std::size_t>(item)] == 0) {
                    bbr12_eligible_[static_cast<std::size_t>(eligible_count++)] =
                        item;
                }
            }
            bbr12_best_task_count_ = 0;
            bbr12_full_load_count_ = 0;
            bbr12_best_score_ = 0;
            enumerate_bbr12(
                0, instance_.capacity, 0, eligible_count, bin_count, 0);
            if (bbr12_best_task_count_ == 0) {
                ++bin_count;
                const int maximum_front = *std::max_element(
                    instance_.front.begin(), instance_.front.end());
                if (bin_count > instance_.n + maximum_front) {
                    throw std::logic_error(
                        "BBR12 Hoffmann heuristic cannot reach an eligible item");
                }
                continue;
            }
            for (int position = 0; position < bbr12_best_task_count_;
                 ++position) {
                const int item = bbr12_best_tasks_[
                    static_cast<std::size_t>(position)];
                bins_[static_cast<std::size_t>(item)] = bin_count;
                bbr12_degrees_[static_cast<std::size_t>(item)] = -1;
                const int begin = bbr12_successor_offset_[
                    static_cast<std::size_t>(item)];
                const int end = bbr12_successor_offset_[
                    static_cast<std::size_t>(item + 1)];
                for (int successor_position = begin;
                     successor_position < end; ++successor_position) {
                    const int successor = bbr12_successors_[
                        static_cast<std::size_t>(successor_position)];
                    --bbr12_degrees_[static_cast<std::size_t>(successor)];
                }
                ++packed;
            }
            ++bin_count;
        }
        HeuristicState result(instance_.n);
        result.bin = bins_;
        result.bin_count = bin_count;
        return result;
    }

    [[nodiscard]] bool packable(int item, int bin, int remaining) const {
        if (remaining < instance_.items[static_cast<std::size_t>(item)].weight) {
            return false;
        }
        for (int p = instance_.pred_offset[static_cast<std::size_t>(item)];
             p < instance_.pred_offset[static_cast<std::size_t>(item + 1)]; ++p) {
            const int predecessor = instance_.pred_from[static_cast<std::size_t>(p)];
            const int predecessor_bin = bins_[static_cast<std::size_t>(predecessor)];
            const int separation =
                instance_.arcs[static_cast<std::size_t>(
                    instance_.pred_arc[static_cast<std::size_t>(p)])]
                    .separation;
            if (predecessor_bin < 0 || bin - predecessor_bin < separation) {
                return false;
            }
        }
        return true;
    }

    void enumerate(int start, int remaining, int bin) {
        if (combination_count_ >= maximum_combinations_) {
            return;
        }
        bool extended = false;
        for (int item = start; item < instance_.n; ++item) {
            if (bins_[static_cast<std::size_t>(item)] >= 0 ||
                !packable(item, bin, remaining)) {
                continue;
            }
            extended = true;
            bins_[static_cast<std::size_t>(item)] = bin;
            combination_.push_back(item);
            enumerate(item + 1,
                      remaining - instance_.items[static_cast<std::size_t>(item)].weight,
                      bin);
            combination_.pop_back();
            bins_[static_cast<std::size_t>(item)] = -1;
            if (combination_count_ >= maximum_combinations_) {
                return;
            }
        }
        if (!extended) {
            double score = 0.0;
            for (const int item : combination_) {
                score += item_score_[static_cast<std::size_t>(item)];
            }
            if (best_combination_.empty() || score > best_score_) {
                best_score_ = score;
                best_combination_ = combination_;
            }
            ++combination_count_;
        }
    }

    [[nodiscard]] HeuristicState construct(double alpha, double beta, double gamma) {
        std::fill(bins_.begin(), bins_.end(), -1);
        item_score_.resize(static_cast<std::size_t>(instance_.n));
        for (int item = 0; item < instance_.n; ++item) {
            const double weight =
                instance_.items[static_cast<std::size_t>(item)].weight;
            item_score_[static_cast<std::size_t>(item)] =
                weight + alpha * (weight + successor_weight_[static_cast<std::size_t>(item)]) +
                beta * successor_count_[static_cast<std::size_t>(item)] - gamma;
        }
        int bin_count = 0;
        int packed = 0;
        maximum_combinations_ = paper_bbr_mode_
                                    ? 1'000
                                    : std::max(1, 100000 / instance_.n);
        while (packed < instance_.n) {
            best_combination_.clear();
            combination_.clear();
            combination_count_ = 0;
            enumerate(0, instance_.capacity, bin_count);
            if (best_combination_.empty()) {
                ++bin_count;
                const int maximum_front = *std::max_element(
                    instance_.front.begin(), instance_.front.end());
                if (bin_count > instance_.n + maximum_front) {
                    throw std::logic_error(
                        "Modified Hoffmann heuristic cannot reach an eligible item");
                }
                continue;
            }
            for (const int item : best_combination_) {
                bins_[static_cast<std::size_t>(item)] = bin_count;
                ++packed;
            }
            ++bin_count;
        }
        HeuristicState result(instance_.n);
        result.bin = bins_;
        result.bin_count = bin_count;
        return result;
    }

    const WorkInstance& instance_;
    bool paper_bbr_mode_ = false;
    bool bbr12_dynamic_mode_ = false;
    bool bbr12_portfolio_mode_ = false;
    int bbr12_portfolio_max_items_ = 200;
    int bbr12_full_load_limit_ = 5'000;
    std::vector<int> successor_weight_;
    std::vector<int> successor_count_;
    std::vector<int> direct_successor_count_;
    std::vector<int> bins_;
    std::vector<double> item_score_;
    std::vector<int> combination_;
    std::vector<int> best_combination_;
    int maximum_combinations_ = 0;
    int combination_count_ = 0;
    double best_score_ = 0.0;
    std::vector<int> bbr12_item_order_;
    std::vector<int> bbr12_rank_;
    std::vector<int> bbr12_successor_offset_;
    std::vector<int> bbr12_successors_;
    std::vector<int> bbr12_degrees_;
    std::vector<int> bbr12_eligible_;
    std::vector<int> bbr12_tasks_;
    std::vector<int> bbr12_best_tasks_;
    std::vector<std::int64_t> bbr12_item_score_;
    int bbr12_best_task_count_ = 0;
    int bbr12_full_load_count_ = 0;
    std::int64_t bbr12_best_score_ = 0;
};

class BoundedDpHeuristic {
public:
    static constexpr int kAlpha = 1'000;
    static constexpr int kBeta = 50;
    static constexpr int kGamma = 1'000;

    BoundedDpHeuristic(const WorkInstance& instance,
                       int lower_bound,
                       const Deadline& deadline)
        : instance_(instance),
          lower_bound_(lower_bound),
          deadline_(deadline),
          blocks_((static_cast<std::size_t>(instance.n) + 63U) / 64U),
          status_(static_cast<std::size_t>(instance.n), 3U),
          remaining_zero_predecessors_(static_cast<std::size_t>(instance.n), 0),
          load_mask_(blocks_, 0U),
          probe_bins_(static_cast<std::size_t>(instance.n), -1) {
        int maximum_separation = 0;
        for (const Arc& arc : instance_.arcs) {
            maximum_separation = std::max(maximum_separation, arc.separation);
        }
        cooldown_levels_ = std::max(0, maximum_separation - 1);
        key_words_ = blocks_ * static_cast<std::size_t>(1 + cooldown_levels_);
        probe_key_.assign(key_words_, 0U);
        build_zero_successors();
        build_branch_order();
        branch_rank_.resize(static_cast<std::size_t>(instance_.n));
        for (int rank = 0; rank < instance_.n; ++rank) {
            branch_rank_[static_cast<std::size_t>(branch_order_[
                static_cast<std::size_t>(rank)])] = rank;
        }
        ready_rank_mask_.assign(blocks_, 0U);
        ready_rank_block_summary_.assign((blocks_ + 63U) / 64U, 0U);
        excluded_ready_item_mask_.assign(blocks_, 0U);
        constexpr int kMaximumFitMaskCapacity = 65'536;
        fit_item_masks_available_ =
            instance_.capacity <= kMaximumFitMaskCapacity;
        if (fit_item_masks_available_) {
            fit_item_masks_.assign(
                static_cast<std::size_t>(instance_.capacity + 1) * blocks_,
                0U);
            for (int item = 0; item < instance_.n; ++item) {
                const int weight =
                    instance_.items[static_cast<std::size_t>(item)].weight;
                std::uint64_t* row = fit_item_masks_.data() +
                    static_cast<std::size_t>(weight) * blocks_;
                row[static_cast<std::size_t>(item) >> 6U] |=
                    std::uint64_t{1} <<
                    (static_cast<unsigned>(item) & 63U);
            }
            for (int capacity = 1; capacity <= instance_.capacity;
                 ++capacity) {
                std::uint64_t* row = fit_item_masks_.data() +
                    static_cast<std::size_t>(capacity) * blocks_;
                const std::uint64_t* previous = row - blocks_;
                for (std::size_t block = 0; block < blocks_; ++block) {
                    row[block] |= previous[block];
                }
            }
        }
        tail_order_.resize(static_cast<std::size_t>(instance_.n));
        std::iota(tail_order_.begin(), tail_order_.end(), 0);
        std::stable_sort(tail_order_.begin(), tail_order_.end(),
                         [&](int lhs, int rhs) {
            const int lhs_tail = instance_.back[static_cast<std::size_t>(lhs)];
            const int rhs_tail = instance_.back[static_cast<std::size_t>(rhs)];
            return lhs_tail != rhs_tail ? lhs_tail > rhs_tail : lhs < rhs;
        });

        current_bins_.reserve(static_cast<std::size_t>(kAlpha) * instance_.n);
        next_bins_.reserve(static_cast<std::size_t>(kAlpha) * instance_.n);
        current_weights_.reserve(kAlpha);
        next_weights_.reserve(kAlpha);
        current_counts_.reserve(kAlpha);
        next_counts_.reserve(kAlpha);
        stage_candidates_.reserve(static_cast<std::size_t>(kAlpha) * kBeta);
        stage_loads_.reserve(static_cast<std::size_t>(kAlpha) * kBeta * blocks_);
        local_candidates_.reserve(kBeta);
        local_loads_.assign(static_cast<std::size_t>(kBeta) * blocks_, 0U);
        candidate_order_.reserve(static_cast<std::size_t>(kAlpha) * kBeta);
        selected_keys_.reserve(static_cast<std::size_t>(kAlpha) * key_words_);
        selected_hashes_.reserve(kAlpha);
        selected_slots_.assign(2048U, -1);
    }

    [[nodiscard]] HeuristicState solve(const HeuristicState& initial) {
        incumbent_ = initial;
        current_bins_.assign(static_cast<std::size_t>(instance_.n), -1);
        current_weights_.assign(1U, 0);
        current_counts_.assign(1U, 0);
        int state_count = 1;
        int depth = 0;
        states_kept_ = 1;

        while (state_count > 0 && incumbent_.bin_count > lower_bound_ &&
               !aborted_ && !deadline_.expired()) {
            stage_candidates_.clear();
            stage_loads_.clear();
            for (int parent = 0; parent < state_count && !aborted_; ++parent) {
                enumerate_parent(parent, depth);
                flush_local_candidates();
            }
            if (incumbent_.bin_count <= depth + 1 || stage_candidates_.empty() ||
                aborted_) {
                break;
            }
            state_count = reduce_stage(depth + 1);
            ++depth;
        }
        completed_ = !aborted_ && !deadline_.expired();
        return incumbent_;
    }

    [[nodiscard]] std::uint64_t states_kept() const noexcept {
        return states_kept_;
    }

    [[nodiscard]] std::uint64_t transitions_generated() const noexcept {
        return transitions_generated_;
    }

    [[nodiscard]] bool completed() const noexcept {
        return completed_;
    }

private:
    struct Candidate {
        int parent = -1;
        int bound = 0;
        std::int64_t idle = 0;
        int assigned_count = 0;
        std::int64_t machine_numerator = 0;
        int longest_tail = 0;
        std::uint64_t serial = 0;
    };

    [[nodiscard]] static std::uint64_t mix(std::uint64_t value) noexcept {
        value ^= value >> 30U;
        value *= 0xbf58476d1ce4e5b9ULL;
        value ^= value >> 27U;
        value *= 0x94d049bb133111ebULL;
        value ^= value >> 31U;
        return value;
    }

    [[nodiscard]] std::uint64_t hash_key(const std::uint64_t* key) const noexcept {
        std::uint64_t hash = mix(static_cast<std::uint64_t>(key_words_));
        for (std::size_t word = 0; word < key_words_; ++word) {
            hash = mix(hash ^ mix(key[word] +
                                  0x9e3779b97f4a7c15ULL * (word + 1U)));
        }
        return hash;
    }

    [[nodiscard]] static bool candidate_better(const Candidate& lhs,
                                                const Candidate& rhs) noexcept {
        return std::tuple{lhs.machine_numerator, lhs.bound, lhs.idle,
                          lhs.longest_tail, -lhs.assigned_count, lhs.serial} <
               std::tuple{rhs.machine_numerator, rhs.bound, rhs.idle,
                          rhs.longest_tail, -rhs.assigned_count, rhs.serial};
    }

    [[nodiscard]] bool bit_is_set(const std::uint64_t* words,
                                  int item) const noexcept {
        return ((words[static_cast<std::size_t>(item) >> 6U] >>
                 (static_cast<unsigned>(item) & 63U)) &
                1U) != 0U;
    }

    void set_bit(std::uint64_t* words, int item) const noexcept {
        words[static_cast<std::size_t>(item) >> 6U] |=
            std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
    }

    void clear_bit(std::uint64_t* words, int item) const noexcept {
        words[static_cast<std::size_t>(item) >> 6U] &=
            ~(std::uint64_t{1} << (static_cast<unsigned>(item) & 63U));
    }

    void set_ready_rank(int item) noexcept {
        const unsigned rank = static_cast<unsigned>(
            branch_rank_[static_cast<std::size_t>(item)]);
        const std::size_t block = rank >> 6U;
        std::uint64_t& word = ready_rank_mask_[block];
        if (word == 0U) {
            ready_rank_block_summary_[block >> 6U] |=
                std::uint64_t{1} << (block & 63U);
        }
        word |= std::uint64_t{1} << (rank & 63U);
    }

    void clear_ready_rank(int item) noexcept {
        const unsigned rank = static_cast<unsigned>(
            branch_rank_[static_cast<std::size_t>(item)]);
        const std::size_t block = rank >> 6U;
        std::uint64_t& word = ready_rank_mask_[block];
        word &= ~(std::uint64_t{1} << (rank & 63U));
        if (word == 0U) {
            ready_rank_block_summary_[block >> 6U] &=
                ~(std::uint64_t{1} << (block & 63U));
        }
    }

    [[nodiscard]] int first_ready_candidate() const noexcept {
        for (std::size_t summary_block = 0;
             summary_block < ready_rank_block_summary_.size();
             ++summary_block) {
            const std::uint64_t summary =
                ready_rank_block_summary_[summary_block];
            if (summary == 0U) {
                continue;
            }
            const std::size_t rank_block = summary_block * 64U +
                std::countr_zero(summary);
            const unsigned rank = static_cast<unsigned>(rank_block * 64U) +
                std::countr_zero(ready_rank_mask_[rank_block]);
            return branch_order_[static_cast<std::size_t>(rank)];
        }
        return -1;
    }

    void set_excluded_ready(int item) noexcept {
        const unsigned value = static_cast<unsigned>(item);
        excluded_ready_item_mask_[value >> 6U] |=
            std::uint64_t{1} << (value & 63U);
    }

    void clear_excluded_ready(int item) noexcept {
        const unsigned value = static_cast<unsigned>(item);
        excluded_ready_item_mask_[value >> 6U] &=
            ~(std::uint64_t{1} << (value & 63U));
    }

    void build_zero_successors() {
        std::vector<int> degree(static_cast<std::size_t>(instance_.n), 0);
        for (const Arc& arc : instance_.arcs) {
            if (arc.separation == 0) {
                ++degree[static_cast<std::size_t>(arc.from)];
            }
        }
        zero_successor_offset_.assign(static_cast<std::size_t>(instance_.n + 1),
                                      0);
        for (int item = 0; item < instance_.n; ++item) {
            zero_successor_offset_[static_cast<std::size_t>(item + 1)] =
                zero_successor_offset_[static_cast<std::size_t>(item)] +
                degree[static_cast<std::size_t>(item)];
        }
        zero_successors_.resize(
            static_cast<std::size_t>(zero_successor_offset_.back()));
        std::vector<int> cursor = zero_successor_offset_;
        for (const Arc& arc : instance_.arcs) {
            if (arc.separation == 0) {
                zero_successors_[static_cast<std::size_t>(
                    cursor[static_cast<std::size_t>(arc.from)]++)] = arc.to;
            }
        }
    }

    void build_branch_order() {
        const std::vector<std::uint64_t> reachability =
            compute_reachability(instance_);
        std::vector<int> dominance(static_cast<std::size_t>(instance_.n), 0);
        std::vector<int> successor_count(static_cast<std::size_t>(instance_.n), 0);
        std::vector<std::int64_t> positional_weight(
            static_cast<std::size_t>(instance_.n), 0);
        const auto subset = [&](const std::uint64_t* lhs,
                                const std::uint64_t* rhs) {
            for (std::size_t block = 0; block < blocks_; ++block) {
                if ((lhs[block] & ~rhs[block]) != 0U) {
                    return false;
                }
            }
            return true;
        };
        for (int item = 0; item < instance_.n; ++item) {
            const std::uint64_t* row = reachability.data() +
                static_cast<std::size_t>(item) * blocks_;
            positional_weight[static_cast<std::size_t>(item)] =
                instance_.items[static_cast<std::size_t>(item)].weight;
            for (std::size_t block = 0; block < blocks_; ++block) {
                std::uint64_t bits = row[block];
                successor_count[static_cast<std::size_t>(item)] +=
                    std::popcount(bits);
                while (bits != 0U) {
                    const unsigned offset = std::countr_zero(bits);
                    const int successor =
                        static_cast<int>(block * 64U + offset);
                    positional_weight[static_cast<std::size_t>(item)] +=
                        instance_.items[static_cast<std::size_t>(successor)].weight;
                    bits &= bits - 1U;
                }
            }
        }
        for (int dominator = 0; dominator < instance_.n; ++dominator) {
            const std::uint64_t* dominator_row = reachability.data() +
                static_cast<std::size_t>(dominator) * blocks_;
            for (int dominated = 0; dominated < instance_.n; ++dominated) {
                if (dominator == dominated ||
                    instance_.items[static_cast<std::size_t>(dominator)].weight <
                        instance_.items[static_cast<std::size_t>(dominated)].weight ||
                    bit_is_set(dominator_row, dominated) ||
                    bit_is_set(reachability.data() +
                                   static_cast<std::size_t>(dominated) * blocks_,
                               dominator) ||
                    !subset(reachability.data() +
                                static_cast<std::size_t>(dominated) * blocks_,
                            dominator_row)) {
                    continue;
                }
                ++dominance[static_cast<std::size_t>(dominator)];
            }
        }
        branch_order_.resize(static_cast<std::size_t>(instance_.n));
        std::iota(branch_order_.begin(), branch_order_.end(), 0);
        std::stable_sort(branch_order_.begin(), branch_order_.end(),
                         [&](int lhs, int rhs) {
            return std::tuple{
                       dominance[static_cast<std::size_t>(lhs)],
                       positional_weight[static_cast<std::size_t>(lhs)],
                       instance_.items[static_cast<std::size_t>(lhs)].weight,
                       successor_count[static_cast<std::size_t>(lhs)], -lhs} >
                   std::tuple{
                       dominance[static_cast<std::size_t>(rhs)],
                       positional_weight[static_cast<std::size_t>(rhs)],
                       instance_.items[static_cast<std::size_t>(rhs)].weight,
                       successor_count[static_cast<std::size_t>(rhs)], -rhs};
        });
    }

    void enumerate_parent(int parent, int depth) {
        current_parent_ = parent;
        current_depth_ = depth;
        current_parent_bins_ = current_bins_.data() +
            static_cast<std::size_t>(parent) * instance_.n;
        current_parent_weight_ =
            current_weights_[static_cast<std::size_t>(parent)];
        current_parent_count_ = current_counts_[static_cast<std::size_t>(parent)];
        std::fill(load_mask_.begin(), load_mask_.end(), 0U);
        std::fill(status_.begin(), status_.end(),
                  static_cast<unsigned char>(3));
        std::fill(remaining_zero_predecessors_.begin(),
                  remaining_zero_predecessors_.end(), 0);
        std::fill(ready_rank_mask_.begin(), ready_rank_mask_.end(), 0U);
        std::fill(ready_rank_block_summary_.begin(),
                  ready_rank_block_summary_.end(), 0U);
        std::fill(excluded_ready_item_mask_.begin(),
                  excluded_ready_item_mask_.end(), 0U);
        local_candidates_.clear();
        current_load_weight_ = 0;
        current_load_count_ = 0;
        remaining_capacity_ = instance_.capacity;
        generated_for_parent_ = 0;
        root_has_ready_ = false;

        for (int item = 0; item < instance_.n; ++item) {
            if (current_parent_bins_[item] >= 0) {
                continue;
            }
            bool blocked = false;
            int pending_zero = 0;
            for (int position =
                     instance_.pred_offset[static_cast<std::size_t>(item)];
                 position <
                     instance_.pred_offset[static_cast<std::size_t>(item + 1)];
                 ++position) {
                const int predecessor = instance_.pred_from[
                    static_cast<std::size_t>(position)];
                const int separation = instance_.arcs[static_cast<std::size_t>(
                    instance_.pred_arc[static_cast<std::size_t>(position)])]
                                           .separation;
                const int predecessor_bin = current_parent_bins_[predecessor];
                if (predecessor_bin < 0) {
                    if (separation == 0) {
                        ++pending_zero;
                    } else {
                        blocked = true;
                        break;
                    }
                } else if (depth - predecessor_bin < separation) {
                    blocked = true;
                    break;
                }
            }
            if (!blocked) {
                status_[static_cast<std::size_t>(item)] = 0U;
                remaining_zero_predecessors_[static_cast<std::size_t>(item)] =
                    pending_zero;
                root_has_ready_ = root_has_ready_ || pending_zero == 0;
                if (pending_zero == 0) {
                    set_ready_rank(item);
                }
            }
        }
        enumerate_loads();
    }

    void enumerate_loads() {
        if (aborted_ || generated_for_parent_ >= kGamma) {
            return;
        }
        if (((++load_search_nodes_) & 1023U) == 0U && deadline_.expired()) {
            aborted_ = true;
            return;
        }
        const int candidate = first_ready_candidate();
        if (candidate < 0) {
            emit_maximal_load();
            return;
        }
        clear_ready_rank(candidate);

        const int weight =
            instance_.items[static_cast<std::size_t>(candidate)].weight;
        if (weight <= remaining_capacity_) {
            status_[static_cast<std::size_t>(candidate)] = 1U;
            set_bit(load_mask_.data(), candidate);
            current_load_weight_ += weight;
            ++current_load_count_;
            remaining_capacity_ -= weight;
            for (int position = zero_successor_offset_[
                                      static_cast<std::size_t>(candidate)];
                 position < zero_successor_offset_[
                                static_cast<std::size_t>(candidate + 1)];
                 ++position) {
                const int successor =
                    zero_successors_[static_cast<std::size_t>(position)];
                if (status_[static_cast<std::size_t>(successor)] <= 2U) {
                    --remaining_zero_predecessors_[
                        static_cast<std::size_t>(successor)];
                    if (status_[static_cast<std::size_t>(successor)] == 0U &&
                        remaining_zero_predecessors_[
                            static_cast<std::size_t>(successor)] == 0) {
                        set_ready_rank(successor);
                    }
                }
            }
            enumerate_loads();
            for (int position = zero_successor_offset_[
                                      static_cast<std::size_t>(candidate)];
                 position < zero_successor_offset_[
                                static_cast<std::size_t>(candidate + 1)];
                 ++position) {
                const int successor =
                    zero_successors_[static_cast<std::size_t>(position)];
                if (status_[static_cast<std::size_t>(successor)] <= 2U) {
                    if (status_[static_cast<std::size_t>(successor)] == 0U &&
                        remaining_zero_predecessors_[
                            static_cast<std::size_t>(successor)] == 0) {
                        clear_ready_rank(successor);
                    }
                    ++remaining_zero_predecessors_[
                        static_cast<std::size_t>(successor)];
                }
            }
            remaining_capacity_ += weight;
            --current_load_count_;
            current_load_weight_ -= weight;
            clear_bit(load_mask_.data(), candidate);
        }
        status_[static_cast<std::size_t>(candidate)] = 2U;
        set_excluded_ready(candidate);
        enumerate_loads();
        clear_excluded_ready(candidate);
        status_[static_cast<std::size_t>(candidate)] = 0U;
        set_ready_rank(candidate);
    }

    void emit_maximal_load() {
        if (fit_item_masks_available_) {
            const std::uint64_t* fit = fit_item_masks_.data() +
                static_cast<std::size_t>(remaining_capacity_) * blocks_;
            for (std::size_t block = 0; block < blocks_; ++block) {
                if ((fit[block] & excluded_ready_item_mask_[block]) != 0U) {
                    return;
                }
            }
        } else {
            for (int item = 0; item < instance_.n; ++item) {
                if (status_[static_cast<std::size_t>(item)] == 2U &&
                    remaining_zero_predecessors_[static_cast<std::size_t>(item)] ==
                        0 &&
                    instance_.items[static_cast<std::size_t>(item)].weight <=
                        remaining_capacity_) {
                    return;
                }
            }
        }
        if (current_load_count_ == 0 && root_has_ready_) {
            return;
        }
        ++generated_for_parent_;
        ++transitions_generated_;
        const int child_depth = current_depth_ + 1;
        const int child_count = current_parent_count_ + current_load_count_;
        const std::int64_t child_weight =
            current_parent_weight_ + current_load_weight_;
        if (child_count == instance_.n) {
            if (child_depth < incumbent_.bin_count) {
                incumbent_.bin.assign(current_parent_bins_,
                                      current_parent_bins_ + instance_.n);
                for (int item = 0; item < instance_.n; ++item) {
                    if (bit_is_set(load_mask_.data(), item)) {
                        incumbent_.bin[static_cast<std::size_t>(item)] =
                            current_depth_;
                    }
                }
                incumbent_.bin_count = child_depth;
            }
            return;
        }
        if (child_depth + 1 >= incumbent_.bin_count) {
            return;
        }

        std::int64_t prefix_weight = 0;
        std::int64_t machine_numerator = 0;
        int longest_tail = 0;
        for (const int item : tail_order_) {
            if (current_parent_bins_[item] >= 0 ||
                bit_is_set(load_mask_.data(), item)) {
                continue;
            }
            prefix_weight +=
                instance_.items[static_cast<std::size_t>(item)].weight;
            const int tail = instance_.back[static_cast<std::size_t>(item)];
            longest_tail = std::max(longest_tail, tail);
            machine_numerator = std::max(
                machine_numerator,
                prefix_weight + static_cast<std::int64_t>(tail) *
                                    instance_.capacity);
        }
        const int bound = std::max(
            lower_bound_,
            child_depth + ceil_div(machine_numerator, instance_.capacity));
        if (bound >= incumbent_.bin_count) {
            return;
        }
        Candidate candidate;
        candidate.parent = current_parent_;
        candidate.bound = bound;
        candidate.idle =
            static_cast<std::int64_t>(child_depth) * instance_.capacity -
            child_weight;
        candidate.assigned_count = child_count;
        candidate.machine_numerator = machine_numerator;
        candidate.longest_tail = longest_tail;
        candidate.serial = next_serial_++;
        keep_local_candidate(candidate);
    }

    void keep_local_candidate(const Candidate& candidate) {
        std::size_t destination = local_candidates_.size();
        if (destination < static_cast<std::size_t>(kBeta)) {
            local_candidates_.push_back(candidate);
        } else {
            std::size_t worst = 0;
            for (std::size_t index = 1; index < local_candidates_.size(); ++index) {
                if (candidate_better(local_candidates_[worst],
                                     local_candidates_[index])) {
                    worst = index;
                }
            }
            if (!candidate_better(candidate, local_candidates_[worst])) {
                return;
            }
            destination = worst;
            local_candidates_[destination] = candidate;
        }
        std::copy(load_mask_.begin(), load_mask_.end(),
                  local_loads_.begin() + destination * blocks_);
    }

    void flush_local_candidates() {
        for (std::size_t index = 0; index < local_candidates_.size(); ++index) {
            stage_candidates_.push_back(local_candidates_[index]);
            const std::uint64_t* load =
                local_loads_.data() + index * blocks_;
            stage_loads_.insert(stage_loads_.end(), load, load + blocks_);
        }
    }

    void build_key(const int* bins, int depth) {
        std::fill(probe_key_.begin(), probe_key_.end(), 0U);
        for (int item = 0; item < instance_.n; ++item) {
            if (bins[item] >= 0) {
                set_bit(probe_key_.data(), item);
                continue;
            }
            int remaining_wait = 0;
            for (int position =
                     instance_.pred_offset[static_cast<std::size_t>(item)];
                 position <
                     instance_.pred_offset[static_cast<std::size_t>(item + 1)];
                 ++position) {
                const int predecessor = instance_.pred_from[
                    static_cast<std::size_t>(position)];
                if (bins[predecessor] < 0) {
                    continue;
                }
                const int separation = instance_.arcs[static_cast<std::size_t>(
                    instance_.pred_arc[static_cast<std::size_t>(position)])]
                                           .separation;
                remaining_wait = std::max(
                    remaining_wait,
                    separation - (depth - bins[predecessor]));
            }
            for (int level = 1;
                 level <= std::min(remaining_wait, cooldown_levels_); ++level) {
                set_bit(probe_key_.data() +
                            static_cast<std::size_t>(level) * blocks_,
                        item);
            }
        }
    }

    [[nodiscard]] bool insert_selected_key(std::uint64_t hash) {
        std::size_t slot = static_cast<std::size_t>(hash) &
                           (selected_slots_.size() - 1U);
        while (true) {
            const int selected = selected_slots_[slot];
            if (selected < 0) {
                selected_slots_[slot] =
                    static_cast<int>(selected_hashes_.size());
                selected_hashes_.push_back(hash);
                selected_keys_.insert(selected_keys_.end(), probe_key_.begin(),
                                      probe_key_.end());
                return true;
            }
            if (selected_hashes_[static_cast<std::size_t>(selected)] == hash &&
                std::equal(
                    probe_key_.begin(), probe_key_.end(),
                    selected_keys_.begin() +
                        static_cast<std::size_t>(selected) * key_words_)) {
                return false;
            }
            slot = (slot + 1U) & (selected_slots_.size() - 1U);
        }
    }

    [[nodiscard]] int reduce_stage(int child_depth) {
        candidate_order_.resize(stage_candidates_.size());
        std::iota(candidate_order_.begin(), candidate_order_.end(), 0);
        std::stable_sort(candidate_order_.begin(), candidate_order_.end(),
                         [&](int lhs, int rhs) {
            return candidate_better(
                stage_candidates_[static_cast<std::size_t>(lhs)],
                stage_candidates_[static_cast<std::size_t>(rhs)]);
        });
        next_bins_.clear();
        next_weights_.clear();
        next_counts_.clear();
        selected_keys_.clear();
        selected_hashes_.clear();
        std::fill(selected_slots_.begin(), selected_slots_.end(), -1);

        for (const int candidate_index : candidate_order_) {
            if (next_weights_.size() >= static_cast<std::size_t>(kAlpha)) {
                break;
            }
            const Candidate& candidate =
                stage_candidates_[static_cast<std::size_t>(candidate_index)];
            const int* parent_bins = current_bins_.data() +
                static_cast<std::size_t>(candidate.parent) * instance_.n;
            std::copy(parent_bins, parent_bins + instance_.n,
                      probe_bins_.begin());
            const std::uint64_t* load = stage_loads_.data() +
                static_cast<std::size_t>(candidate_index) * blocks_;
            int load_weight = 0;
            for (int item = 0; item < instance_.n; ++item) {
                if (bit_is_set(load, item)) {
                    probe_bins_[static_cast<std::size_t>(item)] = child_depth - 1;
                    load_weight +=
                        instance_.items[static_cast<std::size_t>(item)].weight;
                }
            }
            build_key(probe_bins_.data(), child_depth);
            if (!insert_selected_key(hash_key(probe_key_.data()))) {
                continue;
            }
            next_bins_.insert(next_bins_.end(), probe_bins_.begin(),
                              probe_bins_.end());
            next_weights_.push_back(
                current_weights_[static_cast<std::size_t>(candidate.parent)] +
                load_weight);
            next_counts_.push_back(candidate.assigned_count);
        }
        states_kept_ += next_weights_.size();
        current_bins_.swap(next_bins_);
        current_weights_.swap(next_weights_);
        current_counts_.swap(next_counts_);
        return static_cast<int>(current_weights_.size());
    }

    const WorkInstance& instance_;
    int lower_bound_ = 0;
    const Deadline& deadline_;
    std::size_t blocks_ = 0;
    int cooldown_levels_ = 0;
    std::size_t key_words_ = 0;
    HeuristicState incumbent_{0};
    std::vector<int> branch_order_;
    std::vector<int> branch_rank_;
    std::vector<int> tail_order_;
    std::vector<int> zero_successor_offset_;
    std::vector<int> zero_successors_;
    std::vector<unsigned char> status_;
    std::vector<int> remaining_zero_predecessors_;
    std::vector<std::uint64_t> load_mask_;
    std::vector<std::uint64_t> ready_rank_mask_;
    std::vector<std::uint64_t> ready_rank_block_summary_;
    std::vector<std::uint64_t> excluded_ready_item_mask_;
    std::vector<std::uint64_t> fit_item_masks_;
    bool fit_item_masks_available_ = false;
    std::vector<int> current_bins_;
    std::vector<int> next_bins_;
    std::vector<std::int64_t> current_weights_;
    std::vector<std::int64_t> next_weights_;
    std::vector<int> current_counts_;
    std::vector<int> next_counts_;
    std::vector<Candidate> local_candidates_;
    std::vector<std::uint64_t> local_loads_;
    std::vector<Candidate> stage_candidates_;
    std::vector<std::uint64_t> stage_loads_;
    std::vector<int> candidate_order_;
    std::vector<int> probe_bins_;
    std::vector<std::uint64_t> probe_key_;
    std::vector<std::uint64_t> selected_keys_;
    std::vector<std::uint64_t> selected_hashes_;
    std::vector<int> selected_slots_;
    const int* current_parent_bins_ = nullptr;
    int current_parent_ = -1;
    int current_depth_ = 0;
    std::int64_t current_parent_weight_ = 0;
    int current_parent_count_ = 0;
    int current_load_weight_ = 0;
    int current_load_count_ = 0;
    int remaining_capacity_ = 0;
    int generated_for_parent_ = 0;
    bool root_has_ready_ = false;
    bool aborted_ = false;
    bool completed_ = false;
    std::uint64_t next_serial_ = 0;
    std::uint64_t load_search_nodes_ = 0;
    std::uint64_t states_kept_ = 0;
    std::uint64_t transitions_generated_ = 0;
};


[[nodiscard]] PreparedInstance make_unreduced_prepared_instance(
    const Instance& original,
    const WorkInstance& working,
    const HeuristicState& incumbent) {
    PreparedInstance prepared;
    prepared.reversed = working.reversed;
    prepared.original_item_count = original.size();
    prepared.search_to_original.resize(static_cast<std::size_t>(working.n));

    Instance& search = prepared.search_instance;
    search.problem_type = original.problem_type;
    search.order_strength = original.order_strength;
    search.id = original.id;
    search.capacity = working.capacity;
    search.items.reserve(static_cast<std::size_t>(working.n));

    prepared.original_to_search.assign(
        static_cast<std::size_t>(original.size()), -1);
    for (int item = 0; item < working.n; ++item) {
        const int original_item =
            working.items[static_cast<std::size_t>(item)].original_index;
        prepared.search_to_original[static_cast<std::size_t>(item)] =
            original_item;
        prepared.original_to_search[static_cast<std::size_t>(original_item)] =
            item;
        search.items.push_back(Item{
            original_item,
            working.items[static_cast<std::size_t>(item)].weight});
    }

    search.arcs.reserve(original.arcs.size());
    for (const Arc& arc : original.arcs) {
        const int original_from = working.reversed ? arc.to : arc.from;
        const int original_to = working.reversed ? arc.from : arc.to;
        search.arcs.push_back(Arc{
            prepared.original_to_search[static_cast<std::size_t>(original_from)],
            prepared.original_to_search[static_cast<std::size_t>(original_to)],
            arc.separation});
    }
    search.initialize();
    prepared.search_incumbent.bin_of_item = incumbent.bin;
    prepared.search_incumbent.bin_count = incumbent.bin_count;
    return prepared;
}

class StructuralBpppGraph {
public:
    explicit StructuralBpppGraph(int item_count)
        : n_(item_count),
          blocks_((static_cast<std::size_t>(item_count) + 63U) / 64U),
          outgoing_(static_cast<std::size_t>(item_count) * blocks_, 0U),
          incoming_(static_cast<std::size_t>(item_count) * blocks_, 0U),
          active_(blocks_, ~std::uint64_t{0}),
          active_count_(item_count) {
        if (n_ <= 0) {
            throw std::invalid_argument("structured BPP-P graph is empty");
        }
        const unsigned tail = static_cast<unsigned>(n_) & 63U;
        if (tail != 0U) {
            active_.back() = (std::uint64_t{1} << tail) - 1U;
        }
    }

    void add_edge(int from, int to) noexcept {
        if (from == to || !is_active(from) || !is_active(to)) {
            return;
        }
        set_bit(outgoing_row(from), to);
        set_bit(incoming_row(to), from);
    }

    [[nodiscard]] bool is_active(int item) const noexcept {
        const unsigned value = static_cast<unsigned>(item);
        return (active_[value >> 6U] &
                (std::uint64_t{1} << (value & 63U))) != 0U;
    }

    [[nodiscard]] int active_count() const noexcept { return active_count_; }

    [[nodiscard]] std::vector<int> active_predecessors(int item) const {
        return active_neighbors(incoming_row(item));
    }

    [[nodiscard]] std::vector<int> active_successors(int item) const {
        return active_neighbors(outgoing_row(item));
    }

    [[nodiscard]] std::vector<int> sources() const {
        std::vector<int> result;
        result.reserve(static_cast<std::size_t>(active_count_));
        for (int item = 0; item < n_; ++item) {
            if (is_active(item) && !has_active_neighbor(incoming_row(item))) {
                result.push_back(item);
            }
        }
        return result;
    }

    [[nodiscard]] std::vector<int> sinks() const {
        std::vector<int> result;
        result.reserve(static_cast<std::size_t>(active_count_));
        for (int item = 0; item < n_; ++item) {
            if (is_active(item) && !has_active_neighbor(outgoing_row(item))) {
                result.push_back(item);
            }
        }
        return result;
    }

    void remove_full_capacity_item(int item) {
        const std::vector<int> predecessors = active_predecessors(item);
        const std::vector<int> successors = active_successors(item);
        for (const int predecessor : predecessors) {
            for (const int successor : successors) {
                add_edge(predecessor, successor);
            }
        }
        remove_item(item, predecessors, successors);
    }

    void remove_group(const std::vector<int>& items) {
        for (const int item : items) {
            if (!is_active(item)) {
                throw std::logic_error("structured BPP-P group contains an inactive item");
            }
            remove_item(item, active_predecessors(item),
                        active_successors(item));
        }
    }

    template <class Callback>
    void for_each_active_edge(Callback&& callback) const {
        for (int from = 0; from < n_; ++from) {
            if (!is_active(from)) {
                continue;
            }
            const std::uint64_t* row = outgoing_row(from);
            for (std::size_t block = 0; block < blocks_; ++block) {
                std::uint64_t value = row[block] & active_[block];
                while (value != 0U) {
                    const unsigned bit = std::countr_zero(value);
                    callback(from,
                             static_cast<int>(block * 64U + bit));
                    value &= value - 1U;
                }
            }
        }
    }

private:
    static void set_bit(std::uint64_t* row, int item) noexcept {
        const unsigned value = static_cast<unsigned>(item);
        row[value >> 6U] |= std::uint64_t{1} << (value & 63U);
    }

    static void clear_bit(std::uint64_t* row, int item) noexcept {
        const unsigned value = static_cast<unsigned>(item);
        row[value >> 6U] &= ~(std::uint64_t{1} << (value & 63U));
    }

    [[nodiscard]] std::uint64_t* outgoing_row(int item) noexcept {
        return outgoing_.data() + static_cast<std::size_t>(item) * blocks_;
    }

    [[nodiscard]] const std::uint64_t* outgoing_row(int item) const noexcept {
        return outgoing_.data() + static_cast<std::size_t>(item) * blocks_;
    }

    [[nodiscard]] std::uint64_t* incoming_row(int item) noexcept {
        return incoming_.data() + static_cast<std::size_t>(item) * blocks_;
    }

    [[nodiscard]] const std::uint64_t* incoming_row(int item) const noexcept {
        return incoming_.data() + static_cast<std::size_t>(item) * blocks_;
    }

    [[nodiscard]] bool has_active_neighbor(const std::uint64_t* row) const noexcept {
        for (std::size_t block = 0; block < blocks_; ++block) {
            if ((row[block] & active_[block]) != 0U) {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] std::vector<int> active_neighbors(
        const std::uint64_t* row) const {
        std::vector<int> result;
        for (std::size_t block = 0; block < blocks_; ++block) {
            std::uint64_t value = row[block] & active_[block];
            while (value != 0U) {
                const unsigned bit = std::countr_zero(value);
                result.push_back(static_cast<int>(block * 64U + bit));
                value &= value - 1U;
            }
        }
        return result;
    }

    void remove_item(int item,
                     const std::vector<int>& predecessors,
                     const std::vector<int>& successors) {
        for (const int predecessor : predecessors) {
            clear_bit(outgoing_row(predecessor), item);
        }
        for (const int successor : successors) {
            clear_bit(incoming_row(successor), item);
        }
        std::fill(outgoing_row(item), outgoing_row(item) + blocks_, 0U);
        std::fill(incoming_row(item), incoming_row(item) + blocks_, 0U);
        clear_bit(active_.data(), item);
        --active_count_;
    }

    int n_ = 0;
    std::size_t blocks_ = 0;
    std::vector<std::uint64_t> outgoing_;
    std::vector<std::uint64_t> incoming_;
    std::vector<std::uint64_t> active_;
    int active_count_ = 0;
};

[[nodiscard]] std::vector<int> map_work_items_to_original(
    const WorkInstance& working,
    const std::vector<int>& work_items) {
    std::vector<int> originals;
    originals.reserve(work_items.size());
    for (const int item : work_items) {
        originals.push_back(
            working.items[static_cast<std::size_t>(item)].original_index);
    }
    return originals;
}

[[nodiscard]] PreparedInstance make_prepared_instance(
    const Instance& original,
    const WorkInstance& working,
    const HeuristicState& incumbent,
    bool enable_structured_preprocessing) {
    const bool bppp_semantics = !original.arcs.empty() &&
        std::all_of(original.arcs.begin(), original.arcs.end(),
                    [](const Arc& arc) { return arc.separation == 1; });
    if (!enable_structured_preprocessing || !bppp_semantics) {
        return make_unreduced_prepared_instance(original, working, incumbent);
    }

    PreparedInstance prepared;
    prepared.reversed = working.reversed;
    prepared.original_item_count = original.size();

    std::vector<int> original_to_work(
        static_cast<std::size_t>(original.size()), -1);
    for (int item = 0; item < working.n; ++item) {
        original_to_work[static_cast<std::size_t>(
            working.items[static_cast<std::size_t>(item)].original_index)] = item;
    }

    StructuralBpppGraph graph(working.n);
    for (const Arc& arc : original.arcs) {
        const int original_from = working.reversed ? arc.to : arc.from;
        const int original_to = working.reversed ? arc.from : arc.to;
        graph.add_edge(
            original_to_work[static_cast<std::size_t>(original_from)],
            original_to_work[static_cast<std::size_t>(original_to)]);
    }

    prepared.structured_preprocessing_enabled = true;
    {
        while (graph.active_count() > 1) {
            int selected = -1;
            for (int item = 0; item < working.n; ++item) {
                if (graph.is_active(item) &&
                    working.items[static_cast<std::size_t>(item)].weight ==
                        working.capacity) {
                    selected = item;
                    break;
                }
            }
            if (selected < 0) {
                break;
            }
            FullCapacityRemoval removal;
            removal.original_item =
                working.items[static_cast<std::size_t>(selected)].original_index;
            removal.predecessor_original_items = map_work_items_to_original(
                working, graph.active_predecessors(selected));
            removal.successor_original_items = map_work_items_to_original(
                working, graph.active_successors(selected));
            prepared.full_capacity_removals.push_back(std::move(removal));
            graph.remove_full_capacity_item(selected);
        }

        while (graph.active_count() > 1) {
            const std::vector<int> sources = graph.sources();
            std::int64_t weight = 0;
            for (const int item : sources) {
                weight += working.items[static_cast<std::size_t>(item)].weight;
            }
            if (sources.empty() ||
                static_cast<int>(sources.size()) == graph.active_count() ||
                weight > working.capacity) {
                break;
            }
            prepared.fixed_prefix_bins.push_back(
                map_work_items_to_original(working, sources));
            graph.remove_group(sources);
        }

        while (graph.active_count() > 1) {
            const std::vector<int> sinks = graph.sinks();
            std::int64_t weight = 0;
            for (const int item : sinks) {
                weight += working.items[static_cast<std::size_t>(item)].weight;
            }
            if (sinks.empty() ||
                static_cast<int>(sinks.size()) == graph.active_count() ||
                weight > working.capacity) {
                break;
            }
            prepared.fixed_suffix_bins.push_back(
                map_work_items_to_original(working, sinks));
            graph.remove_group(sinks);
        }
    }

    prepared.fixed_bin_offset = static_cast<int>(
        prepared.full_capacity_removals.size() +
        prepared.fixed_prefix_bins.size() + prepared.fixed_suffix_bins.size());
    prepared.search_to_original.reserve(
        static_cast<std::size_t>(graph.active_count()));

    Instance& search = prepared.search_instance;
    search.problem_type = original.problem_type;
    search.order_strength = original.order_strength;
    search.id = original.id;
    search.capacity = working.capacity;
    search.items.reserve(static_cast<std::size_t>(graph.active_count()));

    prepared.original_to_search.assign(
        static_cast<std::size_t>(original.size()), -1);
    std::vector<int> work_to_search(static_cast<std::size_t>(working.n), -1);
    for (int item = 0; item < working.n; ++item) {
        if (!graph.is_active(item)) {
            continue;
        }
        const int original_item =
            working.items[static_cast<std::size_t>(item)].original_index;
        const int search_item = static_cast<int>(search.items.size());
        work_to_search[static_cast<std::size_t>(item)] = search_item;
        prepared.search_to_original.push_back(original_item);
        prepared.original_to_search[static_cast<std::size_t>(original_item)] =
            search_item;
        search.items.push_back(Item{
            original_item,
            working.items[static_cast<std::size_t>(item)].weight});
    }

    graph.for_each_active_edge([&](int from, int to) {
        search.arcs.push_back(Arc{
            work_to_search[static_cast<std::size_t>(from)],
            work_to_search[static_cast<std::size_t>(to)], 1});
    });
    search.initialize();

    std::vector<int> compact_bin(static_cast<std::size_t>(incumbent.bin_count), -1);
    if (prepared.fixed_bin_offset == 0) {
        std::iota(compact_bin.begin(), compact_bin.end(), 0);
    }
    for (int item = 0; item < working.n; ++item) {
        if (prepared.fixed_bin_offset > 0 && graph.is_active(item)) {
            compact_bin[static_cast<std::size_t>(
                incumbent.bin[static_cast<std::size_t>(item)])] = 0;
        }
    }
    int residual_bins = incumbent.bin_count;
    if (prepared.fixed_bin_offset > 0) {
        residual_bins = 0;
        for (int& bin : compact_bin) {
            if (bin == 0) {
                bin = residual_bins++;
            }
        }
    }
    prepared.search_incumbent.bin_of_item.assign(search.items.size(), -1);
    for (int item = 0; item < working.n; ++item) {
        const int search_item = work_to_search[static_cast<std::size_t>(item)];
        if (search_item >= 0) {
            prepared.search_incumbent.bin_of_item[
                static_cast<std::size_t>(search_item)] = compact_bin[
                    static_cast<std::size_t>(
                        incumbent.bin[static_cast<std::size_t>(item)])];
        }
    }
    prepared.search_incumbent.bin_count = residual_bins;
    std::string diagnostic;
    if (!check_assignment(search, prepared.search_incumbent, &diagnostic)) {
        throw std::logic_error(
            "structured BPP-P residual incumbent is invalid: " + diagnostic);
    }
    if (residual_bins + prepared.fixed_bin_offset > incumbent.bin_count) {
        throw std::logic_error(
            "structured BPP-P reduction did not preserve the incumbent bound");
    }
    return prepared;
}

[[nodiscard]] Assignment map_to_original(const WorkInstance& instance,
                                         const HeuristicState& state) {
    Assignment assignment;
    assignment.bin_count = state.bin_count;
    assignment.bin_of_item.assign(static_cast<std::size_t>(instance.n), -1);
    for (int item = 0; item < instance.n; ++item) {
        int bin = state.bin[static_cast<std::size_t>(item)];
        if (instance.reversed) {
            bin = state.bin_count - bin - 1;
        }
        assignment.bin_of_item[static_cast<std::size_t>(
            instance.items[static_cast<std::size_t>(item)].original_index)] = bin;
    }
    return assignment;
}

[[nodiscard]] HeuristicState map_from_original(
    const WorkInstance& instance,
    const Assignment& assignment) {
    if (assignment.bin_count <= 0 ||
        static_cast<int>(assignment.bin_of_item.size()) != instance.n) {
        throw std::invalid_argument("invalid assignment for working-instance mapping");
    }
    HeuristicState state(instance.n);
    state.bin_count = assignment.bin_count;
    state.remaining_capacity.assign(
        static_cast<std::size_t>(state.bin_count), instance.capacity);
    for (int item = 0; item < instance.n; ++item) {
        const int original_item =
            instance.items[static_cast<std::size_t>(item)].original_index;
        const int original_bin =
            assignment.bin_of_item[static_cast<std::size_t>(original_item)];
        if (original_bin < 0 || original_bin >= state.bin_count) {
            throw std::invalid_argument(
                "original assignment contains an invalid bin index");
        }
        const int working_bin = instance.reversed
                                    ? state.bin_count - original_bin - 1
                                    : original_bin;
        state.bin[static_cast<std::size_t>(item)] = working_bin;
        int& remaining =
            state.remaining_capacity[static_cast<std::size_t>(working_bin)];
        remaining -= instance.items[static_cast<std::size_t>(item)].weight;
        if (remaining < 0) {
            throw std::invalid_argument(
                "original assignment exceeds the working capacity");
        }
    }
    return state;
}

}

int compute_initial_dff_lower_bound(const std::vector<int>& weights,
                                    int capacity) {
    if (capacity <= 0) {
        throw std::invalid_argument("exact DFF requires a positive capacity");
    }
    for (const int weight : weights) {
        if (weight <= 0 || weight > capacity) {
            throw std::invalid_argument(
                "exact DFF item weight is outside capacity");
        }
    }
    return exact_dff_lower_bound(
        weights.size(), capacity,
        [&](std::size_t item) { return weights[item]; });
}

InitialBoundsResult compute_initial_bounds(
    const Instance& original,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    InitialBoundsResult result;
    std::optional<Deadline> bounded_initialization_deadline;
    Deadline* initialization_deadline = &deadline;
    if (config.initialization_time_limit_seconds > 0.0) {
        bounded_initialization_deadline.emplace(std::min(
            config.initialization_time_limit_seconds,
            deadline.remaining_seconds()));
        initialization_deadline = &*bounded_initialization_deadline;
    }
    const auto preprocessing_start = Clock::now();
    WorkInstance working(original);
    preprocess(working);
    statistics.preprocessing_seconds +=
        std::chrono::duration<double>(Clock::now() - preprocessing_start).count();
    const auto lower_start = Clock::now();
    CertifiedBounds bounds = quick_bounds(working);
    result.lower_bound =
        std::max(bounds.capacity, bounds.precedence_path);
    statistics.lower_bound_seconds +=
        std::chrono::duration<double>(Clock::now() - lower_start).count();

    const std::int64_t java_seed =
        std::int64_t{929} * working.id * working.n *
            static_cast<std::int64_t>(java_string_hash(working.problem_type)) +
        (static_cast<std::int64_t>(config.seed) - 1);
    JavaRandom random(java_seed);
    const auto fit_start = Clock::now();
    FitBasedHeuristic fit(working);
    HeuristicState incumbent = fit.solve(result.lower_bound, random, 20);
    statistics.upper_bound_seconds +=
        std::chrono::duration<double>(Clock::now() - fit_start).count();

    if (incumbent.bin_count > result.lower_bound) {
        const auto start = Clock::now();
        bounds.window_dff = compute_window_dff_lower_bound(working);
        result.lower_bound =
            std::max(result.lower_bound, bounds.window_dff);
        statistics.lower_bound_seconds +=
            std::chrono::duration<double>(Clock::now() - start).count();
    }
    if (incumbent.bin_count > result.lower_bound &&
        !initialization_deadline->expired()) {
        const auto start = Clock::now();
        ModifiedHoffmannHeuristic hoffmann(
            working, config.bbr_initialization_mode,
            config.bbr_initialization_mode &&
                config.bbr_enable_bbr12_mhh,
            config.bbr_initialization_mode &&
                config.bbr_enable_bbr12_mhh_portfolio,
            config.bbr12_mhh_portfolio_max_items,
            config.bbr12_mhh_full_load_limit);
        HeuristicState candidate =
            hoffmann.solve(result.lower_bound, *initialization_deadline);
        statistics.upper_bound_seconds +=
            std::chrono::duration<double>(Clock::now() - start).count();
        if (better_state(candidate, incumbent)) {
            incumbent = std::move(candidate);
        }
    }

    if (config.bbr_initialization_mode &&
        config.bbr_enable_early_exact_probe &&
        incumbent.bin_count > result.lower_bound &&
        !initialization_deadline->expired()) {
        constexpr double kEarlyProbeSeconds = 0.05;
        constexpr std::uint64_t kEarlyProbeStates = 20'000U;
        const double probe_budget = std::min(
            {kEarlyProbeSeconds, initialization_deadline->remaining_seconds(),
             deadline.remaining_seconds()});
        if (probe_budget > 0.0) {
            PreparedInstance probe_prepared = make_prepared_instance(
                original, working, incumbent,
                config.bbr_enable_structured_preprocessing);
            Config probe_config = config;
            probe_config.threads = 1;
            probe_config.bbr_state_limit =
                std::min(config.bbr_state_limit, kEarlyProbeStates);
            probe_config.bbr_enable_root_strengthening = false;
            probe_config.bbr_enable_binlb = false;
            probe_config.bbr_enable_generalized_item_dominance = false;
            Deadline probe_deadline(probe_budget);
            BbrResult probe = run_branch_bound_remember(
                probe_prepared, result.lower_bound, probe_config,
                probe_deadline);
            result.early_bbr_attempted = probe.attempted;
            result.early_bbr_optimal = probe.optimal;
            result.early_bbr_statistics = probe.statistics;

            const Assignment probe_original =
                map_prepared_assignment_to_original(probe_prepared,
                                                    probe.incumbent);
            HeuristicState probe_incumbent =
                map_from_original(working, probe_original);
            const bool probe_improved =
                better_state(probe_incumbent, incumbent);
            if (probe.optimal || probe_improved) {
                incumbent = std::move(probe_incumbent);
            }
            if (probe.optimal) {
                result.lower_bound = incumbent.bin_count;
            }
        }
    }

    if (config.bbr_initialization_mode &&
        config.bbr_enable_initial_bdp &&
        incumbent.bin_count > result.lower_bound &&
        !initialization_deadline->expired()) {
        const auto run_bounded_dp = [&](const Deadline& pass_deadline) {
            const auto start = Clock::now();
            BoundedDpHeuristic bounded_dp(
                working, result.lower_bound, pass_deadline);
            HeuristicState candidate = bounded_dp.solve(incumbent);
            const double elapsed =
                std::chrono::duration<double>(Clock::now() - start).count();
            statistics.initial_bdp_seconds += elapsed;
            statistics.upper_bound_seconds += elapsed;
            statistics.initial_bdp_states += bounded_dp.states_kept();
            statistics.initial_bdp_transitions +=
                bounded_dp.transitions_generated();
            if (better_state(candidate, incumbent)) {
                incumbent = std::move(candidate);
            }
        };

        run_bounded_dp(*initialization_deadline);
    }

    if (config.bbr_initialization_mode &&
        incumbent.bin_count > result.lower_bound) {
        const UbMinusOneWindows windows = compute_ub_minus_one_windows(
            working, incumbent.bin_count - 1);
        if (windows.infeasible) {
            result.lower_bound = incumbent.bin_count;
        } else if (windows.prefer_reverse) {
            const Assignment original_incumbent =
                map_to_original(working, incumbent);
            flip_precedence_graph(working);
            incumbent = map_from_original(working, original_incumbent);
        }
    }

    const auto exact_preparation_start = Clock::now();
    result.prepared = make_prepared_instance(
        original, working, incumbent,
        config.bbr_initialization_mode &&
            config.bbr_enable_structured_preprocessing);
    statistics.preprocessing_seconds += std::chrono::duration<double>(
        Clock::now() - exact_preparation_start).count();
    result.incumbent = map_to_original(working, incumbent);
    std::string diagnostic;
    if (!check_assignment(original, result.incumbent, &diagnostic)) {
        throw std::logic_error("initial assignment is invalid after mapping: " +
                               diagnostic);
    }
    const Assignment expanded_search_incumbent =
        map_prepared_assignment_to_original(result.prepared,
                                            result.prepared.search_incumbent);
    if (!check_assignment(original, expanded_search_incumbent, &diagnostic)) {
        throw std::logic_error(
            "structured BBR incumbent is invalid after expansion: " + diagnostic);
    }
    if (expanded_search_incumbent.bin_count > result.incumbent.bin_count) {
        throw std::logic_error(
            "structured BBR incumbent is worse than the initialization incumbent");
    }
    if (result.lower_bound > result.incumbent.bin_count) {
        throw std::logic_error("initial lower bound exceeds the validated incumbent");
    }
    return result;
}

Assignment map_prepared_assignment_to_original(
    const PreparedInstance& prepared,
    const Assignment& search_assignment) {
    const int search_n = prepared.search_instance.size();
    const int original_n = prepared.original_item_count > 0
                               ? prepared.original_item_count
                               : search_n;
    if (static_cast<int>(prepared.search_to_original.size()) != search_n ||
        static_cast<int>(search_assignment.bin_of_item.size()) != search_n ||
        original_n < search_n ||
        search_assignment.bin_count <= 0) {
        throw std::invalid_argument("invalid prepared assignment mapping input");
    }
    const int expected_offset = static_cast<int>(
        prepared.full_capacity_removals.size() +
        prepared.fixed_prefix_bins.size() + prepared.fixed_suffix_bins.size());
    if (prepared.fixed_bin_offset != expected_offset) {
        throw std::invalid_argument(
            "prepared assignment has an inconsistent fixed-bin offset");
    }

    Assignment original;
    original.bin_count = search_assignment.bin_count;
    original.bin_of_item.assign(static_cast<std::size_t>(original_n), -1);
    for (int item = 0; item < search_n; ++item) {
        int bin = search_assignment.bin_of_item[static_cast<std::size_t>(item)];
        if (bin < 0 || bin >= search_assignment.bin_count) {
            throw std::invalid_argument(
                "prepared assignment contains an invalid bin index");
        }
        const int original_item =
            prepared.search_to_original[static_cast<std::size_t>(item)];
        if (original_item < 0 || original_item >= original_n ||
            original.bin_of_item[static_cast<std::size_t>(original_item)] >= 0) {
            throw std::invalid_argument(
                "prepared assignment contains an invalid original-item mapping");
        }
        original.bin_of_item[static_cast<std::size_t>(original_item)] = bin;
    }

    const int prefix_count =
        static_cast<int>(prepared.fixed_prefix_bins.size());
    if (prefix_count > 0) {
        for (int& bin : original.bin_of_item) {
            if (bin >= 0) {
                bin += prefix_count;
            }
        }
        for (int position = 0; position < prefix_count; ++position) {
            for (const int item : prepared.fixed_prefix_bins[
                     static_cast<std::size_t>(position)]) {
                if (item < 0 || item >= original_n ||
                    original.bin_of_item[static_cast<std::size_t>(item)] >= 0) {
                    throw std::invalid_argument(
                        "prepared prefix reduction contains an invalid item");
                }
                original.bin_of_item[static_cast<std::size_t>(item)] = position;
            }
        }
        original.bin_count += prefix_count;
    }

    for (auto group = prepared.fixed_suffix_bins.rbegin();
         group != prepared.fixed_suffix_bins.rend(); ++group) {
        const int bin = original.bin_count++;
        for (const int item : *group) {
            if (item < 0 || item >= original_n ||
                original.bin_of_item[static_cast<std::size_t>(item)] >= 0) {
                throw std::invalid_argument(
                    "prepared suffix reduction contains an invalid item");
            }
            original.bin_of_item[static_cast<std::size_t>(item)] = bin;
        }
    }

    for (auto removal = prepared.full_capacity_removals.rbegin();
         removal != prepared.full_capacity_removals.rend(); ++removal) {
        int lower = 0;
        int upper = original.bin_count;
        for (const int predecessor : removal->predecessor_original_items) {
            if (predecessor < 0 || predecessor >= original_n ||
                original.bin_of_item[static_cast<std::size_t>(predecessor)] < 0) {
                throw std::invalid_argument(
                    "prepared full-capacity predecessor is unavailable");
            }
            lower = std::max(
                lower,
                original.bin_of_item[static_cast<std::size_t>(predecessor)] + 1);
        }
        for (const int successor : removal->successor_original_items) {
            if (successor < 0 || successor >= original_n ||
                original.bin_of_item[static_cast<std::size_t>(successor)] < 0) {
                throw std::invalid_argument(
                    "prepared full-capacity successor is unavailable");
            }
            upper = std::min(
                upper,
                original.bin_of_item[static_cast<std::size_t>(successor)]);
        }
        if (removal->original_item < 0 ||
            removal->original_item >= original_n || lower > upper ||
            original.bin_of_item[
                static_cast<std::size_t>(removal->original_item)] >= 0) {
            throw std::logic_error(
                "prepared full-capacity item cannot be reinserted safely");
        }
        for (int& bin : original.bin_of_item) {
            if (bin >= lower) {
                ++bin;
            }
        }
        original.bin_of_item[
            static_cast<std::size_t>(removal->original_item)] = lower;
        ++original.bin_count;
    }

    if (original.bin_count !=
        search_assignment.bin_count + prepared.fixed_bin_offset ||
        !original.complete()) {
        throw std::logic_error(
            "prepared assignment expansion did not restore every item and bin");
    }
    if (prepared.reversed) {
        for (int& bin : original.bin_of_item) {
            bin = original.bin_count - bin - 1;
        }
    }
    return original;
}

}
