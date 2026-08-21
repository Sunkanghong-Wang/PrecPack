#include "precpack/algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace precpack {
namespace {

enum class FitRule { kFirst, kBest };

[[nodiscard]] Assignment build_assignment(const Instance& instance,
                                          const std::vector<int>& priority_order,
                                          FitRule fit_rule,
                                          const Deadline* deadline) {
    const int n = instance.size();
    std::vector<int> rank(static_cast<std::size_t>(n), 0);
    for (int position = 0; position < n; ++position) {
        rank[static_cast<std::size_t>(priority_order[static_cast<std::size_t>(position)])] =
            position;
    }

    using ReadyEntry = std::pair<int, int>;
    std::priority_queue<ReadyEntry, std::vector<ReadyEntry>, std::greater<>> ready;
    std::vector<int> remaining_predecessors(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        remaining_predecessors[static_cast<std::size_t>(i)] =
            static_cast<int>(instance.predecessors[static_cast<std::size_t>(i)].size());
        if (remaining_predecessors[static_cast<std::size_t>(i)] == 0) {
            ready.emplace(rank[static_cast<std::size_t>(i)], i);
        }
    }

    Assignment assignment;
    assignment.bin_of_item.assign(static_cast<std::size_t>(n), -1);
    std::vector<int> remaining_capacity;
    int packed = 0;
    while (!ready.empty()) {
        const int item = ready.top().second;
        ready.pop();

        int earliest_bin = 0;
        for (const auto& [predecessor, separation] :
             instance.predecessor_arcs[static_cast<std::size_t>(item)]) {
            earliest_bin = std::max(
                earliest_bin,
                assignment.bin_of_item[static_cast<std::size_t>(predecessor)] +
                    separation);
        }

        int chosen_bin = -1;
        int best_residual = instance.capacity + 1;
        for (int bin = earliest_bin; bin < assignment.bin_count; ++bin) {
            const int residual =
                remaining_capacity[static_cast<std::size_t>(bin)] -
                instance.items[static_cast<std::size_t>(item)].weight;
            if (residual < 0) {
                continue;
            }
            if (fit_rule == FitRule::kFirst) {
                chosen_bin = bin;
                break;
            }
            if (residual < best_residual) {
                best_residual = residual;
                chosen_bin = bin;
            }
        }

        if (chosen_bin < 0) {
            chosen_bin = std::max(assignment.bin_count, earliest_bin);
            while (assignment.bin_count <= chosen_bin) {
                remaining_capacity.push_back(instance.capacity);
                ++assignment.bin_count;
            }
        }
        assignment.bin_of_item[static_cast<std::size_t>(item)] = chosen_bin;
        remaining_capacity[static_cast<std::size_t>(chosen_bin)] -=
            instance.items[static_cast<std::size_t>(item)].weight;
        ++packed;

        for (const int successor : instance.successors[static_cast<std::size_t>(item)]) {
            if (--remaining_predecessors[static_cast<std::size_t>(successor)] == 0) {
                ready.emplace(rank[static_cast<std::size_t>(successor)], successor);
            }
        }
    }
    if (packed != n) {
        throw std::logic_error("constructive heuristic encountered a cyclic graph");
    }

    std::vector<int> load(static_cast<std::size_t>(assignment.bin_count), 0);
    for (int i = 0; i < n; ++i) {
        load[static_cast<std::size_t>(assignment.bin_of_item[static_cast<std::size_t>(i)])] +=
            instance.items[static_cast<std::size_t>(i)].weight;
    }
    bool changed = true;
    while (changed) {
        changed = false;
        for (const int item : instance.topological_order) {
            if (deadline != nullptr && deadline->expired()) {
                changed = false;
                break;
            }
            const int old_bin = assignment.bin_of_item[static_cast<std::size_t>(item)];
            for (int candidate = 0; candidate < old_bin; ++candidate) {
                if (load[static_cast<std::size_t>(candidate)] +
                        instance.items[static_cast<std::size_t>(item)].weight >
                    instance.capacity) {
                    continue;
                }
                bool feasible = true;
                for (const auto& [predecessor, separation] :
                     instance.predecessor_arcs[static_cast<std::size_t>(item)]) {
                    if (candidate - assignment.bin_of_item[static_cast<std::size_t>(
                                          predecessor)] <
                        separation) {
                        feasible = false;
                        break;
                    }
                }
                if (feasible) {
                    for (const auto& [successor, separation] :
                         instance.successor_arcs[static_cast<std::size_t>(item)]) {
                        if (assignment.bin_of_item[static_cast<std::size_t>(successor)] -
                                candidate <
                            separation) {
                            feasible = false;
                            break;
                        }
                    }
                }
                if (!feasible) {
                    continue;
                }
                load[static_cast<std::size_t>(old_bin)] -=
                    instance.items[static_cast<std::size_t>(item)].weight;
                load[static_cast<std::size_t>(candidate)] +=
                    instance.items[static_cast<std::size_t>(item)].weight;
                assignment.bin_of_item[static_cast<std::size_t>(item)] = candidate;
                changed = true;
                break;
            }
        }
    }
    assignment.bin_count =
        1 + *std::max_element(assignment.bin_of_item.begin(), assignment.bin_of_item.end());
    return assignment;
}

[[nodiscard]] std::int64_t load_square_score(const Instance& instance,
                                             const Assignment& assignment) {
    std::vector<std::int64_t> load(static_cast<std::size_t>(assignment.bin_count), 0);
    for (int i = 0; i < instance.size(); ++i) {
        load[static_cast<std::size_t>(assignment.bin_of_item[static_cast<std::size_t>(i)])] +=
            instance.items[static_cast<std::size_t>(i)].weight;
    }
    return std::accumulate(load.begin(), load.end(), std::int64_t{0},
                           [](std::int64_t sum, std::int64_t value) {
                               return sum + value * value;
                           });
}

[[nodiscard]] bool better_assignment(const Instance& instance,
                                     const Assignment& candidate,
                                     const Assignment& incumbent) {
    if (!incumbent.complete()) {
        return true;
    }
    if (candidate.bin_count != incumbent.bin_count) {
        return candidate.bin_count < incumbent.bin_count;
    }
    return load_square_score(instance, candidate) >
           load_square_score(instance, incumbent);
}

}

LowerBounds compute_simple_lower_bounds(const Instance& instance) {
    LowerBounds bounds;
    bounds.capacity = static_cast<int>((instance.total_weight + instance.capacity - 1) /
                                       instance.capacity);
    int longest_path = 0;
    for (int i = 0; i < instance.size(); ++i) {
        longest_path = std::max(longest_path,
                                instance.front[static_cast<std::size_t>(i)] +
                                    instance.back[static_cast<std::size_t>(i)]);
    }
    bounds.precedence = longest_path + 1;
    bounds.combined = std::max(bounds.capacity, bounds.precedence);
    return bounds;
}

Assignment construct_initial_assignment(const Instance& instance,
                                        int known_lower_bound,
                                        std::mt19937& random,
                                        int randomized_trials,
                                        const Deadline* deadline) {
    const int n = instance.size();
    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);

    Assignment best;
    const auto evaluate = [&](const std::vector<int>& candidate_order) {
        for (const FitRule rule : {FitRule::kFirst, FitRule::kBest}) {
            if (best.complete() && deadline != nullptr && deadline->expired()) {
                return true;
            }
            Assignment candidate =
                build_assignment(instance, candidate_order, rule, deadline);
            if (better_assignment(instance, candidate, best)) {
                best = std::move(candidate);
            }
            if (best.bin_count == known_lower_bound) {
                return true;
            }
        }
        return false;
    };

    if (evaluate(order)) {
        return best;
    }
    if (deadline != nullptr && deadline->expired()) {
        return best;
    }
    std::stable_sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        const int lhs_weight =
            instance.items[static_cast<std::size_t>(lhs)].weight;
        const int rhs_weight =
            instance.items[static_cast<std::size_t>(rhs)].weight;
        return lhs_weight != rhs_weight ? lhs_weight > rhs_weight : lhs < rhs;
    });
    if (evaluate(order)) {
        return best;
    }
    for (int trial = 0; trial < randomized_trials; ++trial) {
        if (deadline != nullptr && deadline->expired()) {
            break;
        }
        std::shuffle(order.begin(), order.end(), random);
        if (evaluate(order)) {
            break;
        }
    }

    std::string diagnostic;
    if (!check_assignment(instance, best, &diagnostic)) {
        throw std::logic_error("constructive heuristic returned an invalid assignment: " +
                               diagnostic);
    }
    return best;
}

bool check_assignment(const Instance& instance,
                      const Assignment& assignment,
                      std::string* diagnostic) {
    const auto fail = [&](std::string message) {
        if (diagnostic != nullptr) {
            *diagnostic = std::move(message);
        }
        return false;
    };
    if (assignment.bin_of_item.size() != instance.items.size()) {
        return fail("assignment length differs from item count");
    }
    if (assignment.bin_count <= 0) {
        return fail("bin count is not positive");
    }

    std::vector<std::int64_t> load(static_cast<std::size_t>(assignment.bin_count), 0);
    int maximum_bin = -1;
    for (int i = 0; i < instance.size(); ++i) {
        const int bin = assignment.bin_of_item[static_cast<std::size_t>(i)];
        if (bin < 0 || bin >= assignment.bin_count) {
            return fail("item " + std::to_string(i + 1) + " has an invalid bin index");
        }
        maximum_bin = std::max(maximum_bin, bin);
        load[static_cast<std::size_t>(bin)] +=
            instance.items[static_cast<std::size_t>(i)].weight;
    }
    if (maximum_bin + 1 != assignment.bin_count) {
        return fail("bin count is not the last occupied position plus one");
    }
    for (int bin = 0; bin < assignment.bin_count; ++bin) {
        if (load[static_cast<std::size_t>(bin)] > instance.capacity) {
            return fail("capacity is exceeded in bin " + std::to_string(bin));
        }
    }
    for (const Arc& arc : instance.arcs) {
        if (assignment.bin_of_item[static_cast<std::size_t>(arc.to)] -
                assignment.bin_of_item[static_cast<std::size_t>(arc.from)] <
            arc.separation) {
            std::ostringstream message;
            message << "precedence arc (" << arc.from + 1 << ',' << arc.to + 1 << ','
                    << arc.separation << ") is violated";
            return fail(message.str());
        }
    }
    if (diagnostic != nullptr) {
        diagnostic->clear();
    }
    return true;
}

}
