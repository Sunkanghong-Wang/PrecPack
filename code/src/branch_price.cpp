#include "precpack/branch_price.hpp"

#include "precpack/exact_arithmetic.hpp"

#include <gurobi_c++.h>

#include "gurobi_compat.hpp"
#include "precpack/conflict_graph.hpp"

#include "precpack/algorithms.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace precpack {
namespace {

using Clock = std::chrono::steady_clock;

inline constexpr double kIntegralityTolerance = 1e-6;
inline constexpr double kPhaseOneTolerance = 1e-7;

struct Pattern {
    int bin = -1;
    std::vector<int> items;
    std::vector<std::uint64_t> bits;

    [[nodiscard]] bool contains(int item) const noexcept {
        return ((bits[static_cast<std::size_t>(item) / 64U] >>
                 (static_cast<unsigned>(item) & 63U)) &
                std::uint64_t{1}) != 0U;
    }
};

[[nodiscard]] std::string pattern_key(
    int bin, const std::vector<std::uint64_t>& bits) {
    std::string key(sizeof(bin) + bits.size() * sizeof(std::uint64_t), '\0');
    std::memcpy(key.data(), &bin, sizeof(bin));
    if (!bits.empty()) {
        std::memcpy(key.data() + sizeof(bin), bits.data(),
                    bits.size() * sizeof(std::uint64_t));
    }
    return key;
}

[[nodiscard]] std::uint64_t triple_key(int i, int j, int k) noexcept {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 42U) ^
           (static_cast<std::uint64_t>(static_cast<std::uint32_t>(j)) << 21U) ^
           static_cast<std::uint32_t>(k);
}

[[nodiscard]] int sr_coefficient(
    const Pattern& pattern, const std::array<int, 3>& triple) noexcept {
    const int count = static_cast<int>(pattern.contains(triple[0])) +
                      static_cast<int>(pattern.contains(triple[1])) +
                      static_cast<int>(pattern.contains(triple[2]));
    return count >= 2 ? 1 : 0;
}

class DisjointSet {
public:
    explicit DisjointSet(int n)
        : parent_(static_cast<std::size_t>(n)),
          rank_(static_cast<std::size_t>(n), 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    int find(int value) {
        int root = value;
        while (parent_[static_cast<std::size_t>(root)] != root) {
            root = parent_[static_cast<std::size_t>(root)];
        }
        while (parent_[static_cast<std::size_t>(value)] != value) {
            const int next = parent_[static_cast<std::size_t>(value)];
            parent_[static_cast<std::size_t>(value)] = root;
            value = next;
        }
        return root;
    }

    void unite(int lhs, int rhs) {
        lhs = find(lhs);
        rhs = find(rhs);
        if (lhs == rhs) {
            return;
        }
        if (rank_[static_cast<std::size_t>(lhs)] <
            rank_[static_cast<std::size_t>(rhs)]) {
            std::swap(lhs, rhs);
        }
        parent_[static_cast<std::size_t>(rhs)] = lhs;
        if (rank_[static_cast<std::size_t>(lhs)] ==
            rank_[static_cast<std::size_t>(rhs)]) {
            ++rank_[static_cast<std::size_t>(lhs)];
        }
    }

private:
    std::vector<int> parent_;
    std::vector<unsigned char> rank_;
};

struct PositionDecision {
    int item = -1;
    int bin = -1;
    bool force = false;
};

struct SearchNode {
    std::uint64_t id = 0;
    int depth = 0;
    int lower_bound = 0;
    std::vector<std::pair<int, int>> together;
    std::vector<std::pair<int, int>> separate;
    std::vector<PositionDecision> position_decisions;
    std::shared_ptr<const std::vector<Pattern>> pattern_pool;
    std::shared_ptr<const std::vector<int>> active_arcs;
    std::shared_ptr<const std::vector<std::array<int, 3>>> sr_cuts;
};

struct CompiledNode {
    bool infeasible = false;
    int bin_count = 0;
    std::vector<int> component_of_item;
    std::vector<std::vector<int>> members;
    std::vector<int> weights;
    std::vector<std::vector<int>> conflict_neighbors;
    std::vector<unsigned char> conflict_matrix;
    std::vector<unsigned char> eligible;

    [[nodiscard]] int component_count() const noexcept {
        return static_cast<int>(members.size());
    }

    [[nodiscard]] bool conflicts(int lhs, int rhs) const noexcept {
        const std::size_t count = members.size();
        return conflict_matrix[static_cast<std::size_t>(lhs) * count +
                               static_cast<std::size_t>(rhs)] != 0U;
    }

    [[nodiscard]] bool is_eligible(int component, int bin) const noexcept {
        return eligible[static_cast<std::size_t>(component) *
                            static_cast<std::size_t>(bin_count) +
                        static_cast<std::size_t>(bin)] != 0U;
    }
};

[[nodiscard]] std::vector<std::pair<int, int>> build_base_conflict_edges(
    const Instance& instance) {
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < instance.size(); ++i) {
        for (int j = i + 1; j < instance.size(); ++j) {
            if (instance.separation(i, j) > 0 ||
                instance.separation(j, i) > 0) {
                edges.emplace_back(i, j);
            }
        }
    }
    return edges;
}

[[nodiscard]] CompiledNode compile_node(
    const Instance& instance,
    int bin_count,
    const std::vector<std::pair<int, int>>& base_conflicts,
    const SearchNode& node,
    bool position_indexed = true) {
    const int n = instance.size();
    DisjointSet sets(n);
    for (const auto& [lhs, rhs] : node.together) {
        sets.unite(lhs, rhs);
    }

    CompiledNode result;
    result.bin_count = bin_count;
    result.component_of_item.assign(static_cast<std::size_t>(n), -1);
    std::vector<int> root_to_component(static_cast<std::size_t>(n), -1);
    for (int item = 0; item < n; ++item) {
        const int root = sets.find(item);
        int& component = root_to_component[static_cast<std::size_t>(root)];
        if (component < 0) {
            component = static_cast<int>(result.members.size());
            result.members.emplace_back();
            result.weights.push_back(0);
        }
        result.component_of_item[static_cast<std::size_t>(item)] = component;
        result.members[static_cast<std::size_t>(component)].push_back(item);
        result.weights[static_cast<std::size_t>(component)] +=
            instance.items[static_cast<std::size_t>(item)].weight;
    }

    const int component_count = result.component_count();
    result.conflict_matrix.assign(
        static_cast<std::size_t>(component_count) * component_count, 0U);
    const auto add_conflict = [&](int item_lhs, int item_rhs) {
        const int lhs =
            result.component_of_item[static_cast<std::size_t>(item_lhs)];
        const int rhs =
            result.component_of_item[static_cast<std::size_t>(item_rhs)];
        if (lhs == rhs) {
            result.infeasible = true;
            return;
        }
        result.conflict_matrix[static_cast<std::size_t>(lhs) * component_count +
                               static_cast<std::size_t>(rhs)] = 1U;
        result.conflict_matrix[static_cast<std::size_t>(rhs) * component_count +
                               static_cast<std::size_t>(lhs)] = 1U;
    };
    for (const auto& [lhs, rhs] : base_conflicts) {
        add_conflict(lhs, rhs);
    }
    for (const auto& [lhs, rhs] : node.separate) {
        add_conflict(lhs, rhs);
    }
    result.conflict_neighbors.assign(static_cast<std::size_t>(component_count), {});
    for (int lhs = 0; lhs < component_count; ++lhs) {
        for (int rhs = 0; rhs < component_count; ++rhs) {
            if (result.conflicts(lhs, rhs)) {
                result.conflict_neighbors[static_cast<std::size_t>(lhs)].push_back(rhs);
            }
        }
    }

    result.eligible.assign(
        static_cast<std::size_t>(component_count) * bin_count, 1U);
    for (int component = 0; component < component_count; ++component) {
        if (result.weights[static_cast<std::size_t>(component)] > instance.capacity) {
            result.infeasible = true;
        }
        for (int bin = 0; position_indexed && bin < bin_count; ++bin) {
            for (const int item : result.members[static_cast<std::size_t>(component)]) {
                if (bin < instance.front[static_cast<std::size_t>(item)] ||
                    bin + instance.back[static_cast<std::size_t>(item)] >= bin_count) {
                    result.eligible[static_cast<std::size_t>(component) * bin_count +
                                    static_cast<std::size_t>(bin)] = 0U;
                    break;
                }
            }
        }
    }
    for (const PositionDecision& decision : node.position_decisions) {
        if (!position_indexed) {
            throw std::logic_error(
                "position decision reached the position-free M tree");
        }
        const int component = result.component_of_item[static_cast<std::size_t>(
            decision.item)];
        if (decision.force) {
            for (int bin = 0; bin < bin_count; ++bin) {
                if (bin != decision.bin) {
                    result.eligible[static_cast<std::size_t>(component) * bin_count +
                                    static_cast<std::size_t>(bin)] = 0U;
                }
            }
        } else {
            result.eligible[static_cast<std::size_t>(component) * bin_count +
                            static_cast<std::size_t>(decision.bin)] = 0U;
        }
    }
    for (int component = 0; component < component_count; ++component) {
        bool any = false;
        for (int bin = 0; bin < bin_count; ++bin) {
            any = any || result.is_eligible(component, bin);
        }
        if (!any) {
            result.infeasible = true;
        }
    }
    return result;
}

[[nodiscard]] bool compatible_pattern(const Pattern& pattern,
                                      const CompiledNode& node) {
    if (pattern.bin < 0 || pattern.bin >= node.bin_count) {
        return false;
    }
    std::vector<int> counts(static_cast<std::size_t>(node.component_count()), 0);
    std::vector<int> selected;
    selected.reserve(pattern.items.size());
    for (const int item : pattern.items) {
        const int component =
            node.component_of_item[static_cast<std::size_t>(item)];
        if (counts[static_cast<std::size_t>(component)]++ == 0) {
            selected.push_back(component);
        }
    }
    for (const int component : selected) {
        if (counts[static_cast<std::size_t>(component)] !=
                static_cast<int>(node.members[static_cast<std::size_t>(component)].size()) ||
            !node.is_eligible(component, pattern.bin)) {
            return false;
        }
    }
    for (std::size_t i = 0; i < selected.size(); ++i) {
        for (std::size_t j = i + 1; j < selected.size(); ++j) {
            if (node.conflicts(selected[i], selected[j])) {
                return false;
            }
        }
    }
    return true;
}

class RestrictedMaster {
public:
    struct Duals {
        std::vector<double> item;
        std::vector<double> usage;
        std::vector<double> arc;
        std::vector<double> sr;
    };

    RestrictedMaster(GRBEnv& environment,
                     const Instance& instance,
                     const CompiledNode& compiled,
                     const SearchNode& node,
                     const Config& config,
                     Statistics& statistics)
        : instance_(instance),
          compiled_(compiled),
          bin_count_(compiled.bin_count),
          model_(environment),
          active_arc_position_(instance.arcs.size(), -1),
          set_covering_(config.set_covering_master) {
        model_.set(GRB_IntParam_Threads, 1);
        model_.set(GRB_IntParam_Seed, config.seed);
        model_.set(GRB_IntParam_OutputFlag, config.gurobi_log ? 1 : 0);
        model_.set(GRB_IntParam_Method, 1);
        model_.set(GRB_IntParam_NumericFocus, 2);
        model_.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
        model_.set(GRB_DoubleParam_OptimalityTol, 1e-9);

        y_.reserve(static_cast<std::size_t>(bin_count_));
        for (int bin = 0; bin < bin_count_; ++bin) {
            y_.push_back(model_.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS));
        }
        model_.update();

        item_rows_.reserve(static_cast<std::size_t>(instance.size()));
        for (int item = 0; item < instance.size(); ++item) {
            item_rows_.push_back(
                config.set_covering_master
                    ? model_.addConstr(GRBLinExpr(0.0) >= 1.0)
                    : model_.addConstr(GRBLinExpr(0.0) == 1.0));
        }
        usage_rows_.reserve(static_cast<std::size_t>(bin_count_));
        for (int bin = 0; bin < bin_count_; ++bin) {
            usage_rows_.push_back(model_.addConstr(-y_[static_cast<std::size_t>(bin)] <=
                                                   0.0));
        }
        prefix_rows_.reserve(static_cast<std::size_t>(std::max(0, bin_count_ - 1)));
        for (int bin = 0; bin + 1 < bin_count_; ++bin) {
            prefix_rows_.push_back(model_.addConstr(
                y_[static_cast<std::size_t>(bin)] -
                    y_[static_cast<std::size_t>(bin + 1)] >=
                0.0));
        }
        model_.update();

        for (GRBConstr& row : item_rows_) {
            GRBColumn column;
            column.addTerm(1.0, row);
            artificial_variables_.push_back(
                model_.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, column));
        }
        for (const int arc_index : *node.active_arcs) {
            add_arc_row_internal(arc_index, false);
        }
        for (const auto& triple : *node.sr_cuts) {
            add_sr_cut_internal(triple, false);
        }

        const std::size_t expected = node.pattern_pool->size() + 64U;
        patterns_.reserve(expected);
        variables_.reserve(expected);
        variable_values_.reserve(expected);
        pattern_keys_.reserve(expected * 2U + 1U);
        for (const Pattern& pattern : *node.pattern_pool) {
            if (compatible_pattern(pattern, compiled_)) {
                add_pattern(pattern);
            }
        }
        model_.update();
        ++statistics.phase_one_count;
    }

    [[nodiscard]] bool solve(Clock::time_point end_time, Statistics& statistics) {
        const double remaining =
            std::chrono::duration<double>(end_time - Clock::now()).count();
        if (remaining <= 0.0) {
            timed_out_ = true;
            return false;
        }
        model_.set(GRB_DoubleParam_TimeLimit, remaining);
        const auto start = Clock::now();
        model_.optimize();
        statistics.rmp_seconds +=
            std::chrono::duration<double>(Clock::now() - start).count();
        ++statistics.rmp_count;
        const int status = model_.get(GRB_IntAttr_Status);
        if (status != GRB_OPTIMAL) {
            timed_out_ = status == GRB_TIME_LIMIT || status == GRB_INTERRUPTED ||
                         gurobi_compat::is_work_limit_status(status);
            return false;
        }
        objective_value_ = model_.get(GRB_DoubleAttr_ObjVal);
        variable_values_.resize(variables_.size());
        for (std::size_t i = 0; i < variables_.size(); ++i) {
            variable_values_[i] = variables_[i].get(GRB_DoubleAttr_X);
        }
        y_values_.resize(y_.size());
        for (std::size_t i = 0; i < y_.size(); ++i) {
            y_values_[i] = y_[i].get(GRB_DoubleAttr_X);
        }
        return true;
    }

    void add_pattern(Pattern pattern) {
        const std::string key = pattern_key(pattern.bin, pattern.bits);
        if (!pattern_keys_.insert(key).second) {
            return;
        }
        GRBColumn column;
        for (const int item : pattern.items) {
            column.addTerm(1.0, item_rows_[static_cast<std::size_t>(item)]);
        }
        column.addTerm(1.0, usage_rows_[static_cast<std::size_t>(pattern.bin)]);
        for (std::size_t position = 0; position < active_arc_indices_.size();
             ++position) {
            const Arc& arc = instance_.arcs[static_cast<std::size_t>(
                active_arc_indices_[position])];
            const double coefficient =
                static_cast<double>(pattern.bin + 1) *
                (static_cast<int>(pattern.contains(arc.to)) -
                 static_cast<int>(pattern.contains(arc.from)));
            if (coefficient != 0.0) {
                column.addTerm(coefficient, active_arc_rows_[position]);
            }
        }
        for (std::size_t cut = 0; cut < sr_cuts_.size(); ++cut) {
            if (sr_coefficient(pattern, sr_cuts_[cut]) != 0) {
                column.addTerm(1.0, sr_rows_[cut]);
            }
        }
        variables_.push_back(
            model_.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, column));
        patterns_.push_back(std::move(pattern));
    }

    [[nodiscard]] Duals duals() {
        Duals result;
        result.item.reserve(item_rows_.size());
        for (GRBConstr& row : item_rows_) {
            const double dual = row.get(GRB_DoubleAttr_Pi);
            result.item.push_back(set_covering_ ? std::max(0.0, dual)
                                                : dual);
        }
        result.usage.reserve(usage_rows_.size());
        for (GRBConstr& row : usage_rows_) {
            result.usage.push_back(
                std::min(0.0, row.get(GRB_DoubleAttr_Pi)));
        }
        result.arc.assign(instance_.arcs.size(), 0.0);
        for (std::size_t position = 0; position < active_arc_indices_.size();
             ++position) {
            result.arc[static_cast<std::size_t>(active_arc_indices_[position])] =
                std::max(0.0,
                         active_arc_rows_[position].get(GRB_DoubleAttr_Pi));
        }
        result.sr.reserve(sr_rows_.size());
        for (GRBConstr& row : sr_rows_) {
            result.sr.push_back(std::min(0.0, row.get(GRB_DoubleAttr_Pi)));
        }
        return result;
    }

    [[nodiscard]] double artificial_value() const {
        double value = 0.0;
        for (const GRBVar& variable : artificial_variables_) {
            value += variable.get(GRB_DoubleAttr_X);
        }
        return value;
    }

    void finish_phase_one() {
        for (GRBVar& variable : artificial_variables_) {
            variable.set(GRB_DoubleAttr_UB, 0.0);
            variable.set(GRB_DoubleAttr_Obj, 0.0);
        }
        for (GRBVar& variable : y_) {
            variable.set(GRB_DoubleAttr_Obj, 1.0);
        }
        for (GRBVar& variable : variables_) {
            variable.set(GRB_DoubleAttr_Obj, 0.0);
        }
        phase_one_ = false;
        model_.update();
    }

    void add_arc_rows(const std::vector<int>& arc_indices,
                      Statistics& statistics) {
        bool added = false;
        for (const int index : arc_indices) {
            if (active_arc_position_[static_cast<std::size_t>(index)] < 0) {
                add_arc_row_internal(index, true);
                added = true;
            }
        }
        if (added) {
            start_phase_one(statistics);
        }
    }

    void add_sr_cuts(const std::vector<std::array<int, 3>>& cuts,
                     Statistics& statistics) {
        bool added = false;
        for (const auto& cut : cuts) {
            if (!sr_keys_.contains(triple_key(cut[0], cut[1], cut[2]))) {
                add_sr_cut_internal(cut, true);
                added = true;
            }
        }
        if (added) {
            start_phase_one(statistics);
        }
    }

    [[nodiscard]] std::vector<int> violated_inactive_arcs(
        double tolerance, int maximum_rows) const {
        const std::vector<double> position = item_positions();
        std::vector<std::pair<double, int>> violations;
        for (std::size_t index = 0; index < instance_.arcs.size(); ++index) {
            if (active_arc_position_[index] >= 0) {
                continue;
            }
            const Arc& arc = instance_.arcs[index];
            const double lhs = position[static_cast<std::size_t>(arc.to)] -
                               position[static_cast<std::size_t>(arc.from)];
            const double violation = static_cast<double>(arc.separation) - lhs;
            if (violation > tolerance) {
                violations.emplace_back(violation, static_cast<int>(index));
            }
        }
        std::sort(violations.begin(), violations.end(),
                  [](const auto& lhs, const auto& rhs) {
                      if (lhs.first != rhs.first) {
                          return lhs.first > rhs.first;
                      }
                      return lhs.second < rhs.second;
                  });
        if (static_cast<int>(violations.size()) > maximum_rows) {
            violations.resize(static_cast<std::size_t>(maximum_rows));
        }
        std::vector<int> result;
        result.reserve(violations.size());
        for (const auto& [violation, index] : violations) {
            static_cast<void>(violation);
            result.push_back(index);
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<int, 3>> violated_sr_cuts(
        int maximum_new_cuts, int item_limit, double tolerance) const {
        std::vector<std::array<int, 3>> result;
        const int n = instance_.size();
        if (maximum_new_cuts <= 0 || n > item_limit) {
            return result;
        }
        struct Candidate {
            double violation = 0.0;
            std::array<int, 3> triple{};
        };
        std::vector<Candidate> candidates;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    if (sr_keys_.contains(triple_key(i, j, k))) {
                        continue;
                    }
                    const std::array<int, 3> triple{i, j, k};
                    double lhs = 0.0;
                    for (std::size_t column = 0; column < patterns_.size(); ++column) {
                        if (variable_values_[column] > tolerance &&
                            sr_coefficient(patterns_[column], triple) != 0) {
                            lhs += variable_values_[column];
                        }
                    }
                    if (lhs > 1.0 + tolerance) {
                        candidates.push_back(Candidate{lhs - 1.0, triple});
                    }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end(),
                  [](const Candidate& lhs, const Candidate& rhs) {
                      if (lhs.violation != rhs.violation) {
                          return lhs.violation > rhs.violation;
                      }
                      return lhs.triple < rhs.triple;
                  });
        if (static_cast<int>(candidates.size()) > maximum_new_cuts) {
            candidates.resize(static_cast<std::size_t>(maximum_new_cuts));
        }
        result.reserve(candidates.size());
        for (const Candidate& candidate : candidates) {
            result.push_back(candidate.triple);
        }
        return result;
    }

    [[nodiscard]] std::vector<double> item_bin_values() const {
        std::vector<double> values(
            static_cast<std::size_t>(instance_.size()) * bin_count_, 0.0);
        for (std::size_t column = 0; column < patterns_.size(); ++column) {
            const double lambda = variable_values_[column];
            if (lambda <= kEpsilon) {
                continue;
            }
            const Pattern& pattern = patterns_[column];
            for (const int item : pattern.items) {
                values[static_cast<std::size_t>(item) * bin_count_ +
                       static_cast<std::size_t>(pattern.bin)] += lambda;
            }
        }
        return values;
    }

    [[nodiscard]] std::vector<double> item_positions() const {
        const std::vector<double> item_bin = item_bin_values();
        std::vector<double> positions(static_cast<std::size_t>(instance_.size()),
                                      0.0);
        for (int item = 0; item < instance_.size(); ++item) {
            for (int bin = 0; bin < bin_count_; ++bin) {
                positions[static_cast<std::size_t>(item)] +=
                    static_cast<double>(bin + 1) *
                    item_bin[static_cast<std::size_t>(item) * bin_count_ +
                             static_cast<std::size_t>(bin)];
            }
        }
        return positions;
    }

    [[nodiscard]] double objective_value() const noexcept {
        return objective_value_;
    }
    [[nodiscard]] bool phase_one() const noexcept { return phase_one_; }
    [[nodiscard]] bool timed_out() const noexcept { return timed_out_; }
    [[nodiscard]] int bin_count() const noexcept { return bin_count_; }
    [[nodiscard]] int sr_cut_count() const noexcept {
        return static_cast<int>(sr_cuts_.size());
    }
    [[nodiscard]] const std::vector<Pattern>& patterns() const noexcept {
        return patterns_;
    }
    [[nodiscard]] const std::vector<double>& variable_values() const noexcept {
        return variable_values_;
    }
    [[nodiscard]] const std::vector<int>& active_arcs() const noexcept {
        return active_arc_indices_;
    }
    [[nodiscard]] const std::vector<std::array<int, 3>>& sr_cuts() const noexcept {
        return sr_cuts_;
    }
    [[nodiscard]] const std::unordered_set<std::string>& pattern_keys() const noexcept {
        return pattern_keys_;
    }

private:
    void start_phase_one(Statistics& statistics) {
        for (GRBVar& variable : y_) {
            variable.set(GRB_DoubleAttr_Obj, 0.0);
        }
        for (GRBVar& variable : variables_) {
            variable.set(GRB_DoubleAttr_Obj, 0.0);
        }
        phase_one_ = true;
        model_.update();
        ++statistics.phase_one_count;
    }

    void add_arc_row_internal(int arc_index, bool has_existing_patterns) {
        const Arc& arc = instance_.arcs[static_cast<std::size_t>(arc_index)];
        GRBLinExpr expression = 0.0;
        if (has_existing_patterns) {
            for (std::size_t column = 0; column < patterns_.size(); ++column) {
                const Pattern& pattern = patterns_[column];
                const double coefficient =
                    static_cast<double>(pattern.bin + 1) *
                    (static_cast<int>(pattern.contains(arc.to)) -
                     static_cast<int>(pattern.contains(arc.from)));
                if (coefficient != 0.0) {
                    expression += coefficient * variables_[column];
                }
            }
        }
        active_arc_position_[static_cast<std::size_t>(arc_index)] =
            static_cast<int>(active_arc_indices_.size());
        active_arc_indices_.push_back(arc_index);
        active_arc_rows_.push_back(
            model_.addConstr(expression >= static_cast<double>(arc.separation)));
        model_.update();
        GRBColumn column;
        column.addTerm(1.0, active_arc_rows_.back());
        artificial_variables_.push_back(
            model_.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS,
                          column));
    }

    void add_sr_cut_internal(const std::array<int, 3>& triple,
                             bool has_existing_patterns) {
        GRBLinExpr expression = 0.0;
        if (has_existing_patterns) {
            for (std::size_t column = 0; column < patterns_.size(); ++column) {
                if (sr_coefficient(patterns_[column], triple) != 0) {
                    expression += variables_[column];
                }
            }
        }
        sr_cuts_.push_back(triple);
        sr_rows_.push_back(model_.addConstr(expression <= 1.0));
        sr_keys_.insert(triple_key(triple[0], triple[1], triple[2]));
        model_.update();
        GRBColumn column;
        column.addTerm(-1.0, sr_rows_.back());
        artificial_variables_.push_back(
            model_.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS,
                          column));
    }

    const Instance& instance_;
    const CompiledNode& compiled_;
    int bin_count_ = 0;
    GRBModel model_;
    std::vector<GRBVar> y_;
    std::vector<GRBConstr> item_rows_;
    std::vector<GRBConstr> usage_rows_;
    std::vector<GRBConstr> prefix_rows_;
    std::vector<int> active_arc_position_;
    std::vector<int> active_arc_indices_;
    std::vector<GRBConstr> active_arc_rows_;
    std::vector<std::array<int, 3>> sr_cuts_;
    std::vector<GRBConstr> sr_rows_;
    std::unordered_set<std::uint64_t> sr_keys_;
    std::vector<GRBVar> artificial_variables_;
    std::vector<Pattern> patterns_;
    std::vector<GRBVar> variables_;
    std::vector<double> variable_values_;
    std::vector<double> y_values_;
    std::unordered_set<std::string> pattern_keys_;
    bool phase_one_ = true;
    bool timed_out_ = false;
    bool set_covering_ = false;
    double objective_value_ = 0.0;
};

struct FloatingPricedPattern {
    Pattern pattern;
    double profit = 0.0;
};

struct FloatingPricingResult {
    std::vector<FloatingPricedPattern> patterns;
    bool proven = true;
    std::uint64_t search_nodes = 0;
};

inline constexpr int kMaximumDpChoiceSpecialGroups = 128;

template <typename Profit>
class FractionalPrefixBound {
public:
    template <typename Group>
    void initialize(const std::vector<Group>& groups,
                    int first_count,
                    int last_count) {
        const std::size_t group_count = groups.size();
        if (first_count < 0 || last_count < first_count ||
            static_cast<std::size_t>(last_count) > group_count) {
            throw std::invalid_argument(
                "invalid direct-pricing fractional prefix range");
        }
        first_count_ = first_count;
        last_count_ = last_count;
        for (const Group& group : groups) {
            const bool invalid_profit = !(group.profit > Profit{0}) ||
                (!std::is_integral_v<Profit> &&
                 !std::isfinite(static_cast<double>(group.profit)));
            if (group.weight < 0 || invalid_profit) {
                throw std::invalid_argument(
                    "invalid direct-pricing fractional candidate");
            }
        }
        std::vector<int> density_order(group_count);
        for (std::size_t index = 0; index < group_count; ++index) {
            density_order[index] = static_cast<int>(index);
        }
        std::sort(density_order.begin(), density_order.end(),
                  [&groups](int lhs_index, int rhs_index) {
                      const Group& lhs =
                          groups[static_cast<std::size_t>(lhs_index)];
                      const Group& rhs =
                          groups[static_cast<std::size_t>(rhs_index)];
                      if (lhs.weight == 0 || rhs.weight == 0) {
                          if (lhs.weight != rhs.weight) {
                              return lhs.weight == 0;
                          }
                      } else if constexpr (std::is_integral_v<Profit>) {
                          const int comparison =
                              exact_arithmetic::compare_nonnegative_fractions(
                                  static_cast<std::uint64_t>(lhs.profit),
                                  static_cast<std::uint64_t>(lhs.weight),
                                  static_cast<std::uint64_t>(rhs.profit),
                                  static_cast<std::uint64_t>(rhs.weight));
                          if (comparison != 0) {
                              return comparison > 0;
                          }
                      } else {
                          const long double lhs_cross =
                              static_cast<long double>(lhs.profit) * rhs.weight;
                          const long double rhs_cross =
                              static_cast<long double>(rhs.profit) * lhs.weight;
                          if (lhs_cross != rhs_cross) {
                              return lhs_cross > rhs_cross;
                          }
                      }
                      return lhs_index < rhs_index;
                  });

        std::vector<int> density_position(group_count, -1);
        density_component_.assign(group_count, -1);
        density_weight_.assign(group_count, 0);
        density_profit_.assign(group_count, Profit{0});
        for (std::size_t density = 0; density < group_count; ++density) {
            const int group_index = density_order[density];
            const Group& group =
                groups[static_cast<std::size_t>(group_index)];
            density_position[static_cast<std::size_t>(group_index)] =
                static_cast<int>(density);
            density_component_[density] = group.component;
            density_weight_[density] = group.weight;
            density_profit_[density] = group.profit;
        }

        bit_blocks_ = (group_count + 63U) / 64U;
        const std::size_t row_count = static_cast<std::size_t>(
            last_count_ - first_count_ + 1);
        if (bit_blocks_ != 0U &&
            row_count >
                std::numeric_limits<std::size_t>::max() / bit_blocks_) {
            throw std::length_error(
                "direct-pricing fractional prefix bitset is too large");
        }
        prefix_density_bits_.assign(
            row_count * bit_blocks_, 0U);
        for (int position = 0; position < first_count_; ++position) {
            const std::size_t density = static_cast<std::size_t>(
                density_position[static_cast<std::size_t>(position)]);
            prefix_density_bits_[density / 64U] |=
                std::uint64_t{1} << (density & 63U);
        }
        for (int count = first_count_ + 1; count <= last_count_; ++count) {
            const std::size_t previous_offset =
                static_cast<std::size_t>(count - first_count_ - 1) *
                bit_blocks_;
            const std::size_t current_offset =
                static_cast<std::size_t>(count - first_count_) * bit_blocks_;
            std::copy_n(prefix_density_bits_.data() + previous_offset,
                        bit_blocks_,
                        prefix_density_bits_.data() + current_offset);
            const std::size_t density = static_cast<std::size_t>(
                density_position[static_cast<std::size_t>(count - 1)]);
            prefix_density_bits_[current_offset + density / 64U] |=
                std::uint64_t{1} << (density & 63U);
        }
    }

    [[nodiscard]] Profit upper_bound(
        int count,
        int capacity,
        const std::vector<int>& blocked) const {
        assert(count >= first_count_);
        assert(count <= last_count_);
        assert(capacity >= 0);
        std::int64_t remaining = capacity;
        Profit profit = Profit{0};
        const std::size_t row_offset =
            static_cast<std::size_t>(count - first_count_) * bit_blocks_;
        for (std::size_t block = 0; block < bit_blocks_; ++block) {
            std::uint64_t active = prefix_density_bits_[row_offset + block];
            while (active != 0U) {
                const unsigned bit = std::countr_zero(active);
                const std::size_t density = block * 64U + bit;
                if (blocked[static_cast<std::size_t>(
                        density_component_[density])] != 0) {
                    active &= active - 1U;
                    continue;
                }
                const std::int64_t weight = density_weight_[density];
                const Profit item_profit = density_profit_[density];
                if (weight <= remaining) {
                    remaining -= weight;
                    profit = add(profit, item_profit);
                } else {
                    assert(weight > 0);
                    if constexpr (std::is_integral_v<Profit>) {
                        const std::int64_t fractional =
                            exact_arithmetic::ceil_nonnegative_product_ratio(
                                static_cast<std::int64_t>(item_profit),
                                remaining, weight,
                                "direct-pricing fractional bound overflow");
                        return add(profit, static_cast<Profit>(fractional));
                    } else {
                        return profit +
                               item_profit *
                                   (static_cast<Profit>(remaining) /
                                    static_cast<Profit>(weight));
                    }
                }
                active &= active - 1U;
            }
        }
        return profit;
    }

private:
    [[nodiscard]] static Profit add(Profit lhs, Profit rhs) {
        if constexpr (std::is_integral_v<Profit>) {
            return static_cast<Profit>(exact_arithmetic::checked_add(
                static_cast<std::int64_t>(lhs),
                static_cast<std::int64_t>(rhs),
                "direct-pricing fractional sum overflow"));
        } else {
            return lhs + rhs;
        }
    }

    int first_count_ = 0;
    int last_count_ = 0;
    std::size_t bit_blocks_ = 0U;
    std::vector<int> density_component_;
    std::vector<std::int64_t> density_weight_;
    std::vector<Profit> density_profit_;
    std::vector<std::uint64_t> prefix_density_bits_;
};

class FloatingPricingSearch {
public:
    FloatingPricingSearch(
        const Instance& instance,
        const CompiledNode& compiled,
        const std::vector<std::array<int, 3>>& sr_cuts,
        const std::vector<double>& sr_duals,
        const std::unordered_set<std::string>& existing_patterns,
        int bin,
        std::vector<double> item_profit,
        double base_cost,
        double tolerance,
        int maximum_columns,
        Clock::time_point end_time)
        : instance_(instance),
          compiled_(compiled),
          sr_cuts_(sr_cuts),
          sr_duals_(sr_duals),
          existing_patterns_(existing_patterns),
          bin_(bin),
          base_cost_(base_cost),
          tolerance_(tolerance),
          maximum_columns_(std::max(1, maximum_columns)),
          end_time_(end_time),
          blocked_(static_cast<std::size_t>(compiled.component_count()), 0),
          selected_components_(static_cast<std::size_t>(compiled.component_count()),
                               0U),
          selected_bits_((static_cast<std::size_t>(instance.size()) + 63U) / 64U,
                         0U),
          cut_count_(sr_cuts.size(), 0),
          item_cut_indices_(static_cast<std::size_t>(instance.size())),
          diversity_count_(static_cast<std::size_t>(instance.size()), 0) {
        for (std::size_t cut = 0; cut < sr_cuts_.size(); ++cut) {
            for (const int item : sr_cuts_[cut]) {
                item_cut_indices_[static_cast<std::size_t>(item)].push_back(
                    static_cast<int>(cut));
            }
        }

        std::vector<Group> easy;
        std::vector<Group> special;
        easy.reserve(static_cast<std::size_t>(compiled.component_count()));
        special.reserve(static_cast<std::size_t>(compiled.component_count()));
        for (int component = 0; component < compiled.component_count(); ++component) {
            if (!compiled.is_eligible(component, bin)) {
                continue;
            }
            double profit = 0.0;
            bool touches_cut = false;
            for (const int item :
                 compiled.members[static_cast<std::size_t>(component)]) {
                profit += item_profit[static_cast<std::size_t>(item)];
                touches_cut = touches_cut ||
                              !item_cut_indices_[static_cast<std::size_t>(item)].empty();
            }
            if (profit <= 0.0) {
                continue;
            }
            Group group{component,
                        compiled.weights[static_cast<std::size_t>(component)],
                        profit};
            if (touches_cut ||
                !compiled.conflict_neighbors[static_cast<std::size_t>(component)]
                     .empty()) {
                special.push_back(group);
            } else {
                easy.push_back(group);
            }
        }
        const auto order = [](const Group& lhs, const Group& rhs) {
            if (lhs.weight != rhs.weight) {
                return lhs.weight < rhs.weight;
            }
            if (lhs.profit != rhs.profit) {
                return lhs.profit < rhs.profit;
            }
            return lhs.component < rhs.component;
        };
        std::sort(easy.begin(), easy.end(), order);
        std::sort(special.begin(), special.end(), order);
        easy_count_ = static_cast<int>(easy.size());
        groups_.reserve(easy.size() + special.size());
        groups_.insert(groups_.end(), easy.begin(), easy.end());
        groups_.insert(groups_.end(), special.begin(), special.end());
        use_dp_choice_gate_ =
            static_cast<int>(groups_.size()) - easy_count_ <=
            kMaximumDpChoiceSpecialGroups;
        if (use_dp_choice_gate_) {
            fractional_bound_.initialize(
                groups_, easy_count_, static_cast<int>(groups_.size()));
        }

        const std::size_t width = static_cast<std::size_t>(instance.capacity + 1);
        dp_.assign((groups_.size() + 1U) * width, 0.0);
        if (use_dp_choice_gate_) {
            choice_blocks_ = (groups_.size() + 63U) / 64U;
            dp_choice_bits_.assign(width * choice_blocks_, 0U);
        }
        for (std::size_t count = 1; count <= groups_.size(); ++count) {
            const Group& group = groups_[count - 1U];
            const double* previous = dp_.data() + (count - 1U) * width;
            double* current = dp_.data() + count * width;
            for (int capacity = 0; capacity <= instance.capacity; ++capacity) {
                current[static_cast<std::size_t>(capacity)] =
                    previous[static_cast<std::size_t>(capacity)];
                if (group.weight <= capacity) {
                    const double included =
                        previous[static_cast<std::size_t>(
                            capacity - group.weight)] + group.profit;
                    if (included >
                        current[static_cast<std::size_t>(capacity)]) {
                        current[static_cast<std::size_t>(capacity)] = included;
                        const std::size_t position = count - 1U;
                        if (use_dp_choice_gate_) {
                            dp_choice_bits_[
                                static_cast<std::size_t>(capacity) *
                                    choice_blocks_ +
                                position / 64U] |=
                                std::uint64_t{1} << (position & 63U);
                        }
                    }
                }
            }
        }
        lookup_key_.resize(sizeof(bin_) +
                           selected_bits_.size() * sizeof(std::uint64_t));
        std::memcpy(lookup_key_.data(), &bin_, sizeof(bin_));
    }

    [[nodiscard]] FloatingPricingResult solve() {
        FloatingPricingResult result;
        if (Clock::now() >= end_time_) {
            result.proven = false;
            return result;
        }

        heuristic_mode_ = true;
        heuristic_limit_ = std::max<std::uint64_t>(
            1U, static_cast<std::uint64_t>(std::max(1, instance_.size())) *
                    static_cast<std::uint64_t>(std::max(1, instance_.capacity)) /
                    10U);
        if (use_dp_choice_gate_) {
            search<true>(static_cast<int>(groups_.size()),
                         instance_.capacity, 0.0);
        } else {
            search<false>(static_cast<int>(groups_.size()),
                          instance_.capacity, 0.0);
        }

        heuristic_mode_ = false;
        heuristic_stopped_ = false;
        best_profit_ = base_cost_ + tolerance_;
        best_pattern_.reset();
        if (use_dp_choice_gate_) {
            search<true>(static_cast<int>(groups_.size()),
                         instance_.capacity, 0.0);
        } else {
            search<false>(static_cast<int>(groups_.size()),
                          instance_.capacity, 0.0);
        }

        result.proven = !timed_out_;
        result.search_nodes = nodes_;
        if (timed_out_) {
            return result;
        }

        std::sort(diverse_patterns_.begin(), diverse_patterns_.end(),
                  [](const FloatingPricedPattern& lhs,
                     const FloatingPricedPattern& rhs) {
                      if (lhs.profit != rhs.profit) {
                          return lhs.profit > rhs.profit;
                      }
                      return lhs.pattern.items < rhs.pattern.items;
                  });
        std::vector<int> final_count(static_cast<std::size_t>(instance_.size()), 0);
        result.patterns.reserve(static_cast<std::size_t>(maximum_columns_));

        const auto append_if_new = [&](const FloatingPricedPattern& priced,
                                       std::vector<FloatingPricedPattern>& output,
                                       std::vector<int>& counts) {
            const std::string key = pattern_key(priced.pattern.bin,
                                                priced.pattern.bits);
            for (const FloatingPricedPattern& existing : output) {
                if (pattern_key(existing.pattern.bin, existing.pattern.bits) == key) {
                    return;
                }
            }
            if (!output.empty()) {
                for (const int item : priced.pattern.items) {
                    if (counts[static_cast<std::size_t>(item)] >= 3) {
                        return;
                    }
                }
            }
            for (const int item : priced.pattern.items) {
                ++counts[static_cast<std::size_t>(item)];
            }
            output.push_back(priced);
        };

        if (best_pattern_.has_value()) {
            append_if_new(*best_pattern_, result.patterns, final_count);
        }
        for (const FloatingPricedPattern& priced : diverse_patterns_) {
            if (static_cast<int>(result.patterns.size()) >= maximum_columns_) {
                break;
            }
            append_if_new(priced, result.patterns, final_count);
        }
        return result;
    }

private:
    struct Group {
        int component = -1;
        int weight = 0;
        double profit = 0.0;
    };

    [[nodiscard]] double dp_bound(int count, int capacity) const noexcept {
        const std::size_t width = static_cast<std::size_t>(instance_.capacity + 1);
        return dp_[static_cast<std::size_t>(count) * width +
                   static_cast<std::size_t>(capacity)];
    }

    [[nodiscard]] bool dp_solution_uses_blocked(int count,
                                                int capacity) const noexcept {
        while (count > easy_count_) {
            std::size_t block =
                static_cast<std::size_t>(count - 1) / 64U;
            const std::size_t row =
                static_cast<std::size_t>(capacity) * choice_blocks_;
            std::uint64_t choices = dp_choice_bits_[row + block];
            const unsigned valid_bits = static_cast<unsigned>(count) & 63U;
            if (valid_bits != 0U) {
                choices &= (std::uint64_t{1} << valid_bits) - 1U;
            }
            while (choices == 0U) {
                if (block == 0U) {
                    return false;
                }
                choices = dp_choice_bits_[row + --block];
            }
            const int position = static_cast<int>(
                block * 64U + 63U - std::countl_zero(choices));
            if (position < easy_count_) {
                return false;
            }
            const Group& group =
                groups_[static_cast<std::size_t>(position)];
            if (blocked_[static_cast<std::size_t>(group.component)] != 0) {
                return true;
            }
            capacity -= group.weight;
            count = position;
        }
        return false;
    }

    [[nodiscard]] bool selected_pattern_is_existing() {
        if (!selected_bits_.empty()) {
            std::memcpy(lookup_key_.data() + sizeof(bin_), selected_bits_.data(),
                        selected_bits_.size() * sizeof(std::uint64_t));
        }
        return existing_patterns_.contains(lookup_key_);
    }

    [[nodiscard]] Pattern make_pattern() const {
        Pattern pattern;
        pattern.bin = bin_;
        pattern.bits = selected_bits_;
        for (int component = 0; component < compiled_.component_count(); ++component) {
            if (selected_components_[static_cast<std::size_t>(component)] == 0U) {
                continue;
            }
            const auto& members =
                compiled_.members[static_cast<std::size_t>(component)];
            pattern.items.insert(pattern.items.end(), members.begin(), members.end());
        }
        std::sort(pattern.items.begin(), pattern.items.end());
        return pattern;
    }

    void consider(double profit) {
        if (profit <= base_cost_ + tolerance_ || selected_pattern_is_existing()) {
            return;
        }
        if (heuristic_mode_) {
            if (static_cast<int>(diverse_patterns_.size()) >=
                maximum_columns_ * 4) {
                return;
            }
            Pattern pattern = make_pattern();
            const std::string key = pattern_key(pattern.bin, pattern.bits);
            if (!generated_keys_.insert(key).second) {
                return;
            }
            for (const int item : pattern.items) {
                if (diversity_count_[static_cast<std::size_t>(item)] >= 6) {
                    return;
                }
            }
            for (const int item : pattern.items) {
                ++diversity_count_[static_cast<std::size_t>(item)];
            }
            diverse_patterns_.push_back(
                FloatingPricedPattern{std::move(pattern), profit});
            return;
        }
        if (profit > best_profit_) {
            best_profit_ = profit;
            best_pattern_ = FloatingPricedPattern{make_pattern(), profit};
        }
    }

    template <bool TrackBlockedCount>
    void include_group(const Group& group, double& profit) {
        const int component = group.component;
        selected_components_[static_cast<std::size_t>(component)] = 1U;
        profit += group.profit;
        for (const int item :
             compiled_.members[static_cast<std::size_t>(component)]) {
            selected_bits_[static_cast<std::size_t>(item) / 64U] |=
                std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
            for (const int cut :
                 item_cut_indices_[static_cast<std::size_t>(item)]) {
                if (cut_count_[static_cast<std::size_t>(cut)] == 1) {
                    profit += sr_duals_[static_cast<std::size_t>(cut)];
                }
                ++cut_count_[static_cast<std::size_t>(cut)];
            }
        }
        for (const int conflict :
             compiled_.conflict_neighbors[static_cast<std::size_t>(component)]) {
            if constexpr (TrackBlockedCount) {
                int& block_count = blocked_[static_cast<std::size_t>(conflict)];
                if (block_count++ == 0) {
                    ++blocked_component_count_;
                }
            } else {
                ++blocked_[static_cast<std::size_t>(conflict)];
            }
        }
    }

    template <bool TrackBlockedCount>
    void remove_group(const Group& group, double& profit) {
        const int component = group.component;
        for (const int conflict :
             compiled_.conflict_neighbors[static_cast<std::size_t>(component)]) {
            int& block_count = blocked_[static_cast<std::size_t>(conflict)];
            assert(block_count > 0);
            if constexpr (TrackBlockedCount) {
                if (--block_count == 0) {
                    assert(blocked_component_count_ > 0);
                    --blocked_component_count_;
                }
            } else {
                --block_count;
            }
        }
        for (const int item :
             compiled_.members[static_cast<std::size_t>(component)]) {
            for (const int cut :
                 item_cut_indices_[static_cast<std::size_t>(item)]) {
                --cut_count_[static_cast<std::size_t>(cut)];
                if (cut_count_[static_cast<std::size_t>(cut)] == 1) {
                    profit -= sr_duals_[static_cast<std::size_t>(cut)];
                }
            }
            selected_bits_[static_cast<std::size_t>(item) / 64U] &=
                ~(std::uint64_t{1} << (static_cast<unsigned>(item) & 63U));
        }
        profit -= group.profit;
        selected_components_[static_cast<std::size_t>(component)] = 0U;
    }

    template <bool UseFractionalBound>
    void search(int count, int capacity, double profit) {
        if (timed_out_ || heuristic_stopped_) {
            return;
        }
        ++nodes_;
        if ((nodes_ & 1023U) == 0U && Clock::now() >= end_time_) {
            timed_out_ = true;
            return;
        }
        if (heuristic_mode_ && nodes_ > heuristic_limit_ &&
            !diverse_patterns_.empty()) {
            heuristic_stopped_ = true;
            return;
        }

        consider(profit);
        if (count == 0 || capacity == 0) {
            return;
        }
        const double incumbent = heuristic_mode_ ? base_cost_ + tolerance_
                                                  : best_profit_;
        const double upper_bound = profit + dp_bound(count, capacity);
        const double safety = 1e-11 * std::max(1.0, std::abs(upper_bound));
        if (upper_bound + safety <= incumbent) {
            return;
        }
        if constexpr (UseFractionalBound) {
            if (count > easy_count_ && blocked_component_count_ != 0 &&
                dp_solution_uses_blocked(count, capacity)) {
                const double fractional_upper_bound =
                    profit + fractional_bound_.upper_bound(
                                 count, capacity, blocked_);
                const double fractional_safety =
                    1e-11 *
                    std::max(1.0, std::abs(fractional_upper_bound));
                if (fractional_upper_bound + fractional_safety <= incumbent) {
                    return;
                }
            }
        }

        const Group& group = groups_[static_cast<std::size_t>(count - 1)];
        if (group.weight <= capacity &&
            blocked_[static_cast<std::size_t>(group.component)] == 0) {
            double included_profit = profit;
            include_group<UseFractionalBound>(group, included_profit);
            search<UseFractionalBound>(
                count - 1, capacity - group.weight, included_profit);
            remove_group<UseFractionalBound>(group, included_profit);
        }
        search<UseFractionalBound>(count - 1, capacity, profit);
    }

    const Instance& instance_;
    const CompiledNode& compiled_;
    const std::vector<std::array<int, 3>>& sr_cuts_;
    const std::vector<double>& sr_duals_;
    const std::unordered_set<std::string>& existing_patterns_;
    int bin_ = -1;
    double base_cost_ = 0.0;
    double tolerance_ = 0.0;
    int maximum_columns_ = 1;
    Clock::time_point end_time_;
    std::vector<Group> groups_;
    int easy_count_ = 0;
    std::vector<double> dp_;
    std::size_t choice_blocks_ = 0U;
    std::vector<std::uint64_t> dp_choice_bits_;
    bool use_dp_choice_gate_ = false;
    FractionalPrefixBound<double> fractional_bound_;
    std::vector<int> blocked_;
    int blocked_component_count_ = 0;
    std::vector<unsigned char> selected_components_;
    std::vector<std::uint64_t> selected_bits_;
    std::vector<int> cut_count_;
    std::vector<std::vector<int>> item_cut_indices_;
    std::vector<int> diversity_count_;
    std::vector<FloatingPricedPattern> diverse_patterns_;
    std::optional<FloatingPricedPattern> best_pattern_;
    std::unordered_set<std::string> generated_keys_;
    std::string lookup_key_;
    double best_profit_ = 0.0;
    std::uint64_t nodes_ = 0;
    std::uint64_t heuristic_limit_ = 0;
    bool heuristic_mode_ = false;
    bool heuristic_stopped_ = false;
    bool timed_out_ = false;
};

[[nodiscard]] std::int64_t checked_add(std::int64_t lhs,
                                       std::int64_t rhs) {
    return exact_arithmetic::checked_add(
        lhs, rhs, "fixed-point pricing value overflow");
}

struct IntegerPricingResult {
    bool proven = true;
    std::int64_t maximum_profit = 0;
    std::uint64_t search_nodes = 0;
};

class IntegerPricingSearch {
public:
    IntegerPricingSearch(const Instance& instance,
                         const CompiledNode& compiled,
                         const std::vector<std::array<int, 3>>& sr_cuts,
                         const std::vector<std::int64_t>& sr_duals,
                         int bin,
                         std::vector<std::int64_t> item_profit,
                         Clock::time_point end_time)
        : instance_(instance),
          compiled_(compiled),
          sr_cuts_(sr_cuts),
          sr_duals_(sr_duals),
          end_time_(end_time),
          blocked_(static_cast<std::size_t>(compiled.component_count()), 0),
          cut_count_(sr_cuts.size(), 0),
          item_cut_indices_(static_cast<std::size_t>(instance.size())) {
        for (std::size_t cut = 0; cut < sr_cuts_.size(); ++cut) {
            for (const int item : sr_cuts_[cut]) {
                item_cut_indices_[static_cast<std::size_t>(item)].push_back(
                    static_cast<int>(cut));
            }
        }
        const int component_count = compiled.component_count();
        std::vector<std::int64_t> component_profit(
            static_cast<std::size_t>(component_count), 0);
        std::vector<unsigned char> touches_cut(
            static_cast<std::size_t>(component_count), 0U);
        for (int component = 0; component < component_count; ++component) {
            for (const int item :
                 compiled.members[static_cast<std::size_t>(component)]) {
                component_profit[static_cast<std::size_t>(component)] =
                    checked_add(
                        component_profit[static_cast<std::size_t>(component)],
                        item_profit[static_cast<std::size_t>(item)]);
                touches_cut[static_cast<std::size_t>(component)] =
                    static_cast<unsigned char>(
                        touches_cut[static_cast<std::size_t>(component)] != 0U ||
                        !item_cut_indices_[static_cast<std::size_t>(item)].empty());
            }
        }

        std::vector<Group> easy;
        std::vector<Group> special;
        easy.reserve(static_cast<std::size_t>(component_count));
        special.reserve(static_cast<std::size_t>(component_count));
        for (int component = component_count - 1; component >= 0; --component) {
            if (!compiled.is_eligible(component, bin) ||
                component_profit[static_cast<std::size_t>(component)] <= 0 ||
                touches_cut[static_cast<std::size_t>(component)] != 0U ||
                !compiled.conflict_neighbors[static_cast<std::size_t>(component)]
                     .empty()) {
                continue;
            }
            easy.push_back(Group{
                component,
                compiled.weights[static_cast<std::size_t>(component)],
                component_profit[static_cast<std::size_t>(component)]});
        }
        std::vector<int> special_order;
        special_order.reserve(static_cast<std::size_t>(component_count));
        std::vector<unsigned char> seen(static_cast<std::size_t>(component_count),
                                        0U);
        const auto append_special = [&](int component) {
            if (seen[static_cast<std::size_t>(component)] == 0U) {
                seen[static_cast<std::size_t>(component)] = 1U;
                special_order.push_back(component);
            }
        };
        for (int lhs = 0; lhs < component_count; ++lhs) {
            for (int rhs = lhs + 1; rhs < component_count; ++rhs) {
                if (compiled.conflicts(lhs, rhs)) {
                    append_special(lhs);
                    append_special(rhs);
                }
            }
        }
        for (int component = 0; component < component_count; ++component) {
            if (touches_cut[static_cast<std::size_t>(component)] != 0U) {
                append_special(component);
            }
        }
        for (const int component : special_order) {
            if (!compiled.is_eligible(component, bin) ||
                component_profit[static_cast<std::size_t>(component)] <= 0) {
                continue;
            }
            special.push_back(Group{
                component,
                compiled.weights[static_cast<std::size_t>(component)],
                component_profit[static_cast<std::size_t>(component)]});
        }
        std::stable_sort(special.begin(), special.end(),
                         [](const Group& lhs, const Group& rhs) {
                             return lhs.weight < rhs.weight;
                         });
        easy_count_ = static_cast<int>(easy.size());
        groups_.reserve(easy.size() + special.size());
        groups_.insert(groups_.end(), easy.begin(), easy.end());
        groups_.insert(groups_.end(), special.begin(), special.end());
        use_dp_choice_gate_ =
            static_cast<int>(groups_.size()) - easy_count_ <=
            kMaximumDpChoiceSpecialGroups;
        if (use_dp_choice_gate_) {
            fractional_bound_.initialize(
                groups_, easy_count_, static_cast<int>(groups_.size()));
        }

        const std::size_t width = static_cast<std::size_t>(instance.capacity + 1);
        dp_.assign((groups_.size() + 1U) * width, 0);
        if (use_dp_choice_gate_) {
            choice_blocks_ = (groups_.size() + 63U) / 64U;
            dp_choice_bits_.assign(width * choice_blocks_, 0U);
        }
        for (std::size_t count = 1; count <= groups_.size(); ++count) {
            const Group& group = groups_[count - 1U];
            const std::int64_t* previous =
                dp_.data() + (count - 1U) * width;
            std::int64_t* current = dp_.data() + count * width;
            for (int capacity = 0; capacity <= instance.capacity; ++capacity) {
                current[static_cast<std::size_t>(capacity)] =
                    previous[static_cast<std::size_t>(capacity)];
                if (group.weight <= capacity) {
                    const std::int64_t included = checked_add(
                        previous[static_cast<std::size_t>(
                            capacity - group.weight)],
                        group.profit);
                    if (included >
                        current[static_cast<std::size_t>(capacity)]) {
                        current[static_cast<std::size_t>(capacity)] = included;
                        const std::size_t position = count - 1U;
                        if (use_dp_choice_gate_) {
                            dp_choice_bits_[
                                static_cast<std::size_t>(capacity) *
                                    choice_blocks_ +
                                position / 64U] |=
                                std::uint64_t{1} << (position & 63U);
                        }
                    }
                }
            }
        }
    }

    [[nodiscard]] IntegerPricingResult solve() {
        IntegerPricingResult result;
        if (Clock::now() >= end_time_) {
            result.proven = false;
            return result;
        }
        if (use_dp_choice_gate_) {
            search<true>(static_cast<int>(groups_.size()),
                         instance_.capacity, 0);
        } else {
            search<false>(static_cast<int>(groups_.size()),
                          instance_.capacity, 0);
        }
        result.proven = !timed_out_;
        result.maximum_profit = best_profit_;
        result.search_nodes = nodes_;
        return result;
    }

private:
    struct Group {
        int component = -1;
        int weight = 0;
        std::int64_t profit = 0;
    };

    [[nodiscard]] std::int64_t dp_bound(int count, int capacity) const noexcept {
        const std::size_t width = static_cast<std::size_t>(instance_.capacity + 1);
        return dp_[static_cast<std::size_t>(count) * width +
                   static_cast<std::size_t>(capacity)];
    }

    [[nodiscard]] bool dp_solution_uses_blocked(int count,
                                                int capacity) const noexcept {
        while (count > easy_count_) {
            std::size_t block =
                static_cast<std::size_t>(count - 1) / 64U;
            const std::size_t row =
                static_cast<std::size_t>(capacity) * choice_blocks_;
            std::uint64_t choices = dp_choice_bits_[row + block];
            const unsigned valid_bits = static_cast<unsigned>(count) & 63U;
            if (valid_bits != 0U) {
                choices &= (std::uint64_t{1} << valid_bits) - 1U;
            }
            while (choices == 0U) {
                if (block == 0U) {
                    return false;
                }
                choices = dp_choice_bits_[row + --block];
            }
            const int position = static_cast<int>(
                block * 64U + 63U - std::countl_zero(choices));
            if (position < easy_count_) {
                return false;
            }
            const Group& group =
                groups_[static_cast<std::size_t>(position)];
            if (blocked_[static_cast<std::size_t>(group.component)] != 0) {
                return true;
            }
            capacity -= group.weight;
            count = position;
        }
        return false;
    }

    template <bool TrackBlockedCount>
    void include_group(const Group& group, std::int64_t& profit) {
        profit = checked_add(profit, group.profit);
        for (const int item :
             compiled_.members[static_cast<std::size_t>(group.component)]) {
            for (const int cut :
                 item_cut_indices_[static_cast<std::size_t>(item)]) {
                if (cut_count_[static_cast<std::size_t>(cut)] == 1) {
                    profit = checked_add(
                        profit, sr_duals_[static_cast<std::size_t>(cut)]);
                }
                ++cut_count_[static_cast<std::size_t>(cut)];
            }
        }
        for (const int conflict : compiled_.conflict_neighbors[
                 static_cast<std::size_t>(group.component)]) {
            if constexpr (TrackBlockedCount) {
                int& block_count = blocked_[static_cast<std::size_t>(conflict)];
                if (block_count++ == 0) {
                    ++blocked_component_count_;
                }
            } else {
                ++blocked_[static_cast<std::size_t>(conflict)];
            }
        }
    }

    template <bool TrackBlockedCount>
    void remove_group(const Group& group, std::int64_t& profit) {
        for (const int conflict : compiled_.conflict_neighbors[
                 static_cast<std::size_t>(group.component)]) {
            int& block_count = blocked_[static_cast<std::size_t>(conflict)];
            assert(block_count > 0);
            if constexpr (TrackBlockedCount) {
                if (--block_count == 0) {
                    assert(blocked_component_count_ > 0);
                    --blocked_component_count_;
                }
            } else {
                --block_count;
            }
        }
        for (const int item :
             compiled_.members[static_cast<std::size_t>(group.component)]) {
            for (const int cut :
                 item_cut_indices_[static_cast<std::size_t>(item)]) {
                --cut_count_[static_cast<std::size_t>(cut)];
                if (cut_count_[static_cast<std::size_t>(cut)] == 1) {
                    profit = checked_add(
                        profit, -sr_duals_[static_cast<std::size_t>(cut)]);
                }
            }
        }
        profit = checked_add(profit, -group.profit);
    }

    template <bool UseFractionalBound>
    void search(int count, int capacity, std::int64_t profit) {
        if (timed_out_) {
            return;
        }
        ++nodes_;
        if ((nodes_ & 1023U) == 0U && Clock::now() >= end_time_) {
            timed_out_ = true;
            return;
        }
        best_profit_ = std::max(best_profit_, profit);
        if (count == 0 || capacity == 0) {
            return;
        }
        const std::int64_t upper_bound =
            checked_add(profit, dp_bound(count, capacity));
        if (upper_bound <= best_profit_) {
            return;
        }
        if constexpr (UseFractionalBound) {
            if (count > easy_count_ && blocked_component_count_ != 0 &&
                dp_solution_uses_blocked(count, capacity)) {
                const std::int64_t fractional_upper_bound = checked_add(
                    profit,
                    fractional_bound_.upper_bound(count, capacity, blocked_));
                if (fractional_upper_bound <= best_profit_) {
                    return;
                }
            }
        }
        if (count == easy_count_) {
            best_profit_ = std::max(best_profit_, upper_bound);
            return;
        }
        const Group& group = groups_[static_cast<std::size_t>(count - 1)];
        if (group.weight <= capacity &&
            blocked_[static_cast<std::size_t>(group.component)] == 0) {
            std::int64_t included_profit = profit;
            include_group<UseFractionalBound>(group, included_profit);
            search<UseFractionalBound>(
                count - 1, capacity - group.weight, included_profit);
            remove_group<UseFractionalBound>(group, included_profit);
        }
        search<UseFractionalBound>(count - 1, capacity, profit);
    }

    const Instance& instance_;
    const CompiledNode& compiled_;
    const std::vector<std::array<int, 3>>& sr_cuts_;
    const std::vector<std::int64_t>& sr_duals_;
    Clock::time_point end_time_;
    std::vector<Group> groups_;
    int easy_count_ = 0;
    std::vector<std::int64_t> dp_;
    std::size_t choice_blocks_ = 0U;
    std::vector<std::uint64_t> dp_choice_bits_;
    bool use_dp_choice_gate_ = false;
    FractionalPrefixBound<std::int64_t> fractional_bound_;
    std::vector<int> blocked_;
    int blocked_component_count_ = 0;
    std::vector<int> cut_count_;
    std::vector<std::vector<int>> item_cut_indices_;
    std::int64_t best_profit_ = 0;
    std::uint64_t nodes_ = 0;
    bool timed_out_ = false;
};

inline constexpr double kMMaximumObjectiveMultiplier = 400.0;
inline constexpr double kMScaledObjectiveTarget = 900000.0;
inline constexpr std::int64_t kMFixedPointBudget =
    std::int64_t{1} << 62U;

[[nodiscard]] double m_objective_multiplier(int upper_bound) {
    const double range_limited = std::floor(
        kMScaledObjectiveTarget /
        static_cast<double>(std::max(1, upper_bound)));
    return std::clamp(range_limited, 1.0,
                      kMMaximumObjectiveMultiplier);
}

class MRestrictedMaster {
public:
    struct Duals {
        std::vector<double> item;
        double pattern_count = 0.0;
        std::vector<double> sr;
    };

    MRestrictedMaster(GRBEnv& environment,
                      const Instance& instance,
                      const CompiledNode& compiled,
                      const SearchNode& node,
                      int incumbent_value,
                      const Config& config,
                      bool unbounded_relaxation_variables)
        : instance_(instance),
          compiled_(compiled),
          objective_multiplier_(m_objective_multiplier(incumbent_value)),
          horizon_upper_bound_(unbounded_relaxation_variables
                                   ? GRB_INFINITY
                                   : static_cast<double>(incumbent_value)),
          pattern_upper_bound_(unbounded_relaxation_variables
                                   ? GRB_INFINITY
                                   : 1.0),
          model_(environment) {
        model_.set(GRB_IntParam_Threads, 1);
        model_.set(GRB_IntParam_Seed, config.seed);
        model_.set(GRB_IntParam_OutputFlag, config.gurobi_log ? 1 : 0);
        model_.set(GRB_IntParam_Presolve, 0);
        model_.set(GRB_IntParam_Method, 0);
        model_.set(GRB_IntParam_NumericFocus, 2);
        model_.set(GRB_IntParam_ScaleFlag, 0);
        model_.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
        model_.set(GRB_DoubleParam_OptimalityTol, 1e-9);

        horizon_ = model_.addVar(
            static_cast<double>(node.lower_bound),
            horizon_upper_bound_, objective_multiplier_,
            GRB_CONTINUOUS);
        model_.update();

        item_rows_.reserve(static_cast<std::size_t>(instance.size()));
        for (int item = 0; item < instance.size(); ++item) {
            item_rows_.push_back(model_.addConstr(GRBLinExpr(0.0) >= 1.0));
        }
        pattern_count_row_ = model_.addConstr(GRBLinExpr(horizon_) >= 0.0);
        for (const auto& triple : *node.sr_cuts) {
            add_sr_cut_internal(triple, false);
        }
        model_.update();

        const std::size_t expected =
            node.pattern_pool->size() + compiled.members.size() + 64U;
        patterns_.reserve(expected);
        variables_.reserve(expected);
        variable_values_.reserve(expected);
        pattern_keys_.reserve(expected * 2U + 1U);
        for (const Pattern& inherited : *node.pattern_pool) {
            Pattern pattern = inherited;
            pattern.bin = 0;
            if (compatible_pattern(pattern, compiled_)) {
                static_cast<void>(add_pattern(std::move(pattern)));
            }
        }

        const std::size_t bit_count =
            (static_cast<std::size_t>(instance.size()) + 63U) / 64U;
        for (const std::vector<int>& members : compiled.members) {
            Pattern component_pattern;
            component_pattern.bin = 0;
            component_pattern.items = members;
            component_pattern.bits.assign(bit_count, 0U);
            for (const int item : members) {
                component_pattern.bits[static_cast<std::size_t>(item) / 64U] |=
                    std::uint64_t{1} <<
                    (static_cast<unsigned>(item) & 63U);
            }
            static_cast<void>(add_pattern(std::move(component_pattern)));
        }
        model_.update();
    }

    [[nodiscard]] bool solve(Clock::time_point end_time,
                             Statistics& statistics) {
        const double remaining =
            std::chrono::duration<double>(end_time - Clock::now()).count();
        if (remaining <= 0.0) {
            timed_out_ = true;
            return false;
        }
        model_.set(GRB_DoubleParam_TimeLimit, remaining);
        const auto start = Clock::now();
        model_.optimize();
        statistics.rmp_seconds +=
            std::chrono::duration<double>(Clock::now() - start).count();
        ++statistics.rmp_count;
        const int status = model_.get(GRB_IntAttr_Status);
        if (status != GRB_OPTIMAL) {
            timed_out_ = status != GRB_INFEASIBLE;
            return false;
        }
        objective_value_ =
            model_.get(GRB_DoubleAttr_ObjVal) / objective_multiplier_;
        variable_values_.resize(variables_.size());
        for (std::size_t index = 0; index < variables_.size(); ++index) {
            variable_values_[index] =
                variables_[index].get(GRB_DoubleAttr_X);
        }
        return true;
    }

    [[nodiscard]] bool add_pattern(Pattern pattern) {
        pattern.bin = 0;
        const std::string key = pattern_key(0, pattern.bits);
        if (!pattern_keys_.insert(key).second) {
            return false;
        }
        GRBColumn column;
        for (const int item : pattern.items) {
            column.addTerm(1.0,
                           item_rows_[static_cast<std::size_t>(item)]);
        }
        column.addTerm(-1.0, pattern_count_row_);
        for (std::size_t cut = 0; cut < sr_cuts_.size(); ++cut) {
            if (sr_coefficient(pattern, sr_cuts_[cut]) != 0) {
                column.addTerm(1.0, sr_rows_[cut]);
            }
        }
        variables_.push_back(model_.addVar(
            0.0, pattern_upper_bound_, 0.0, GRB_CONTINUOUS, column));
        patterns_.push_back(std::move(pattern));
        return true;
    }

    void update() { model_.update(); }

    [[nodiscard]] Duals duals() {
        Duals result;
        result.item.reserve(item_rows_.size());
        for (GRBConstr& row : item_rows_) {
            const double raw = row.get(GRB_DoubleAttr_Pi);
            if (!std::isfinite(raw)) {
                throw std::runtime_error(
                    "Gurobi returned a non-finite M-master item dual");
            }
            result.item.push_back(
                std::max(0.0, raw / objective_multiplier_));
        }
        const double raw_count =
            pattern_count_row_.get(GRB_DoubleAttr_Pi);
        if (!std::isfinite(raw_count)) {
            throw std::runtime_error(
                "Gurobi returned a non-finite M-master count dual");
        }
        result.pattern_count =
            std::max(0.0, raw_count / objective_multiplier_);
        result.sr.reserve(sr_rows_.size());
        for (GRBConstr& row : sr_rows_) {
            const double raw = row.get(GRB_DoubleAttr_Pi);
            if (!std::isfinite(raw)) {
                throw std::runtime_error(
                    "Gurobi returned a non-finite M-master SR dual");
            }
            result.sr.push_back(
                std::min(0.0, raw / objective_multiplier_));
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<int, 3>> violated_sr_cuts(
        int maximum_new_cuts, int item_limit, double tolerance) const {
        std::vector<std::array<int, 3>> result;
        const int n = instance_.size();
        if (maximum_new_cuts <= 0 || n > item_limit) {
            return result;
        }
        struct Candidate {
            double violation = 0.0;
            std::array<int, 3> triple{};
        };
        std::vector<Candidate> candidates;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    if (sr_keys_.contains(triple_key(i, j, k))) {
                        continue;
                    }
                    const std::array<int, 3> triple{i, j, k};
                    double lhs = 0.0;
                    for (std::size_t column = 0; column < patterns_.size();
                         ++column) {
                        if (variable_values_[column] > tolerance &&
                            sr_coefficient(patterns_[column], triple) != 0) {
                            lhs += variable_values_[column];
                        }
                    }
                    if (lhs > 1.0 + tolerance) {
                        candidates.push_back(
                            Candidate{lhs - 1.0, triple});
                    }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end(),
                  [](const Candidate& lhs, const Candidate& rhs) {
                      if (lhs.violation != rhs.violation) {
                          return lhs.violation > rhs.violation;
                      }
                      return lhs.triple < rhs.triple;
                  });
        if (static_cast<int>(candidates.size()) > maximum_new_cuts) {
            candidates.resize(static_cast<std::size_t>(maximum_new_cuts));
        }
        result.reserve(candidates.size());
        for (const Candidate& candidate : candidates) {
            result.push_back(candidate.triple);
        }
        return result;
    }

    void add_sr_cuts(const std::vector<std::array<int, 3>>& cuts) {
        for (const auto& cut : cuts) {
            if (!sr_keys_.contains(triple_key(cut[0], cut[1], cut[2]))) {
                add_sr_cut_internal(cut, true);
            }
        }
        model_.update();
    }

    [[nodiscard]] bool timed_out() const noexcept { return timed_out_; }
    [[nodiscard]] double objective_value() const noexcept {
        return objective_value_;
    }
    [[nodiscard]] double objective_multiplier() const noexcept {
        return objective_multiplier_;
    }
    [[nodiscard]] int sr_cut_count() const noexcept {
        return static_cast<int>(sr_cuts_.size());
    }
    [[nodiscard]] const std::vector<Pattern>& patterns() const noexcept {
        return patterns_;
    }
    [[nodiscard]] const std::vector<double>& variable_values() const noexcept {
        return variable_values_;
    }
    [[nodiscard]] const std::vector<std::array<int, 3>>& sr_cuts() const
        noexcept {
        return sr_cuts_;
    }
    [[nodiscard]] const std::unordered_set<std::string>& pattern_keys() const
        noexcept {
        return pattern_keys_;
    }

private:
    void add_sr_cut_internal(const std::array<int, 3>& triple,
                             bool has_existing_patterns) {
        GRBLinExpr expression = 0.0;
        if (has_existing_patterns) {
            for (std::size_t column = 0; column < patterns_.size(); ++column) {
                if (sr_coefficient(patterns_[column], triple) != 0) {
                    expression += variables_[column];
                }
            }
        }
        sr_cuts_.push_back(triple);
        sr_rows_.push_back(model_.addConstr(expression <= 1.0));
        sr_keys_.insert(triple_key(triple[0], triple[1], triple[2]));
    }

    const Instance& instance_;
    const CompiledNode& compiled_;
    double objective_multiplier_ = 1.0;
    double horizon_upper_bound_ = GRB_INFINITY;
    double pattern_upper_bound_ = 1.0;
    GRBModel model_;
    GRBVar horizon_;
    std::vector<GRBConstr> item_rows_;
    GRBConstr pattern_count_row_;
    std::vector<std::array<int, 3>> sr_cuts_;
    std::vector<GRBConstr> sr_rows_;
    std::unordered_set<std::uint64_t> sr_keys_;
    std::vector<Pattern> patterns_;
    std::vector<GRBVar> variables_;
    std::vector<double> variable_values_;
    std::unordered_set<std::string> pattern_keys_;
    bool timed_out_ = false;
    double objective_value_ = 0.0;
};

struct MPricingRoundResult {
    bool proven = true;
    bool added_columns = false;
};

[[nodiscard]] MPricingRoundResult run_m_pricing_round(
    const Instance& instance,
    const CompiledNode& compiled,
    MRestrictedMaster& master,
    const MRestrictedMaster::Duals& duals,
    const Config& config,
    Clock::time_point end_time,
    Statistics& statistics) {
    MPricingRoundResult result;
    const auto start = Clock::now();
    FloatingPricingSearch pricing(
        instance, compiled, master.sr_cuts(), duals.sr,
        master.pattern_keys(), 0, duals.item, duals.pattern_count, 0.0,
        config.max_columns_per_pricing, end_time);
    FloatingPricingResult priced = pricing.solve();
    ++statistics.pricing_count;
    statistics.pricing_search_nodes += priced.search_nodes;
    if (!priced.proven) {
        result.proven = false;
    } else {
        for (FloatingPricedPattern& generated : priced.patterns) {
            if (!master.add_pattern(std::move(generated.pattern))) {
                throw std::logic_error(
                    "M pricing returned a duplicate RLMP column");
            }
            ++statistics.generated_columns;
            result.added_columns = true;
        }
        if (result.added_columns) {
            master.update();
        }
    }
    statistics.pricing_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

struct MScaledDuals {
    std::int64_t scale = 1;
    std::vector<std::int64_t> item;
    std::vector<std::int64_t> sr;
    std::int64_t sum = 0;
};

[[nodiscard]] MScaledDuals scale_m_duals(
    const MRestrictedMaster::Duals& duals) {
    long double magnitude = 0.0L;
    for (const double value : duals.item) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("non-finite M-master item dual");
        }
        magnitude += static_cast<long double>(std::max(0.0, value));
    }
    for (const double value : duals.sr) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("non-finite M-master SR dual");
        }
        magnitude += static_cast<long double>(std::abs(std::min(0.0, value)));
    }
    magnitude = std::max(1.0L, magnitude);
    const long double candidate = std::floor(
        static_cast<long double>(kMFixedPointBudget) / magnitude);
    if (!std::isfinite(candidate) || candidate < 1.0L) {
        throw std::overflow_error(
            "M-master duals are too large for int64 fixed point");
    }
    std::int64_t scale = std::min(
        kMFixedPointBudget, static_cast<std::int64_t>(candidate));
    for (;;) {
        MScaledDuals result;
        result.scale = scale;
        std::int64_t absolute_sum = 0;
        result.item.reserve(duals.item.size());
        for (const double raw : duals.item) {
            const long double scaled =
                static_cast<long double>(std::max(0.0, raw)) *
                static_cast<long double>(scale);
            if (!std::isfinite(scaled) ||
                scaled >= static_cast<long double>(
                              std::numeric_limits<std::int64_t>::max())) {
                throw std::overflow_error("M-master dual scaling overflow");
            }
            const std::int64_t value =
                static_cast<std::int64_t>(scaled);
            result.item.push_back(value);
            result.sum = checked_add(result.sum, value);
            absolute_sum = checked_add(absolute_sum, value);
        }
        result.sr.reserve(duals.sr.size());
        for (const double raw : duals.sr) {
            const long double scaled =
                static_cast<long double>(std::min(0.0, raw)) *
                static_cast<long double>(scale);
            if (!std::isfinite(scaled) ||
                scaled <= static_cast<long double>(
                              std::numeric_limits<std::int64_t>::min())) {
                throw std::overflow_error("M-master SR dual scaling overflow");
            }
            const std::int64_t value =
                static_cast<std::int64_t>(std::floor(scaled));
            result.sr.push_back(value);
            result.sum = checked_add(result.sum, value);
            absolute_sum = exact_arithmetic::checked_subtract(
                absolute_sum, value,
                "M-master fixed-point accumulation overflow");
        }
        if (absolute_sum <= kMFixedPointBudget) {
            return result;
        }
        if (scale == 1) {
            throw std::overflow_error(
                "M-master fixed-point accumulation overflow");
        }
        const long double ratio =
            static_cast<long double>(kMFixedPointBudget) /
            static_cast<long double>(absolute_sum);
        const std::int64_t next = std::max<std::int64_t>(
            1, static_cast<std::int64_t>(std::floor(
                   static_cast<long double>(scale) * ratio)));
        scale = next < scale ? next : scale - 1;
    }
}

[[nodiscard]] int ceil_m_ratio(std::int64_t numerator,
                               std::int64_t denominator) {
    if (denominator <= 0) {
        throw std::invalid_argument(
            "nonpositive M-master certificate denominator");
    }
    return exact_arithmetic::ceil_ratio_to_int(
        numerator, denominator, "M-master bound does not fit int");
}

struct MSafeBoundResult {
    bool proven = true;
    int integer_bound = 0;
};

[[nodiscard]] MSafeBoundResult compute_m_safe_bound(
    const Instance& instance,
    const CompiledNode& compiled,
    const MRestrictedMaster& master,
    const MRestrictedMaster::Duals& floating_duals,
    Clock::time_point end_time,
    Statistics& statistics) {
    MSafeBoundResult result;
    const MScaledDuals duals = scale_m_duals(floating_duals);
    const auto start = Clock::now();
    IntegerPricingSearch pricing(
        instance, compiled, master.sr_cuts(), duals.sr, 0, duals.item,
        end_time);
    const IntegerPricingResult priced = pricing.solve();
    statistics.pricing_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    ++statistics.pricing_count;
    statistics.pricing_search_nodes += priced.search_nodes;
    if (!priced.proven) {
        result.proven = false;
        return result;
    }
    const std::int64_t denominator =
        std::max(duals.scale, priced.maximum_profit);
    result.integer_bound = ceil_m_ratio(duals.sum, denominator);
    return result;
}

enum class MNodeState {
    kConverged,
    kInfeasible,
    kTimedOut,
};

struct MNodeSolveResult {
    MNodeState state = MNodeState::kTimedOut;
    int certified_lower_bound = 0;
    double lp_value = 0.0;
    std::shared_ptr<const std::vector<Pattern>> pattern_pool;
    std::shared_ptr<const std::vector<std::array<int, 3>>> sr_cuts;
    std::vector<double> variable_values;
    std::uint64_t column_count = 0;
};

[[nodiscard]] MNodeSolveResult solve_m_node(
    GRBEnv& environment,
    const Instance& instance,
    const CompiledNode& compiled,
    const SearchNode& node,
    int incumbent_value,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics,
    bool unbounded_relaxation_variables = false) {
    MNodeSolveResult result;
    result.certified_lower_bound = node.lower_bound;
    const auto start = Clock::now();
    MRestrictedMaster master(environment, instance, compiled, node,
                             incumbent_value, config,
                             unbounded_relaxation_variables);
    int sr_rounds = 0;

    for (int iteration = 0; iteration < config.max_cg_iterations;
         ++iteration) {
        ++statistics.cg_iterations;
        if (deadline.expired()) {
            result.state = MNodeState::kTimedOut;
            break;
        }
        if (!master.solve(deadline.end_time(), statistics)) {
            result.state = master.timed_out() ? MNodeState::kTimedOut
                                              : MNodeState::kInfeasible;
            break;
        }
        const MRestrictedMaster::Duals duals = master.duals();
        const MPricingRoundResult pricing = run_m_pricing_round(
            instance, compiled, master, duals, config, deadline.end_time(),
            statistics);
        if (!pricing.proven) {
            result.state = MNodeState::kTimedOut;
            break;
        }
        if (pricing.added_columns) {
            continue;
        }

        if (node.depth == 0 && config.enable_sr_cuts &&
            sr_rounds < config.max_sr_rounds_per_node &&
            master.sr_cut_count() < config.max_sr_cuts) {
            const int maximum_new = std::min(
                20, config.max_sr_cuts - master.sr_cut_count());
            const int item_limit =
                std::max(100, config.sr_enumeration_item_limit);
            const auto cuts = master.violated_sr_cuts(
                maximum_new, item_limit,
                config.row_violation_tolerance);
            if (!cuts.empty()) {
                ++sr_rounds;
                statistics.generated_sr_rows += cuts.size();
                master.add_sr_cuts(cuts);
                continue;
            }
        }

        result.lp_value = master.objective_value();
        const MSafeBoundResult safe = compute_m_safe_bound(
            instance, compiled, master, duals, deadline.end_time(),
            statistics);
        if (!safe.proven) {
            result.state = MNodeState::kTimedOut;
            break;
        }
        result.certified_lower_bound =
            std::max(node.lower_bound, safe.integer_bound);
        result.pattern_pool =
            std::make_shared<const std::vector<Pattern>>(master.patterns());
        result.sr_cuts =
            std::make_shared<const std::vector<std::array<int, 3>>>(
                master.sr_cuts());
        result.variable_values = master.variable_values();
        result.state = MNodeState::kConverged;
        break;
    }

    statistics.cg_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    ++statistics.cg_count;
    result.column_count = master.patterns().size();
    return result;
}

struct ScaledDuals {
    std::int64_t scale = 1;
    std::vector<std::int64_t> item;
    std::vector<std::int64_t> arc;
    std::vector<std::int64_t> sr;
};

[[nodiscard]] ScaledDuals scale_duals(
    const Instance& instance,
    int bin_count,
    const RestrictedMaster::Duals& duals) {
    long double magnitude = static_cast<long double>(bin_count + 1);
    for (const double value : duals.item) {
        magnitude += std::abs(static_cast<long double>(value));
    }
    for (std::size_t index = 0; index < duals.arc.size(); ++index) {
        const long double value =
            std::abs(static_cast<long double>(duals.arc[index]));
        magnitude += value *
                     static_cast<long double>(
                         2 * (bin_count + 1) +
                         std::abs(instance.arcs[index].separation));
    }
    for (const double value : duals.sr) {
        magnitude += static_cast<long double>(bin_count + 1) *
                     std::abs(static_cast<long double>(value));
    }
    magnitude = std::max<long double>(1.0L, magnitude);

    std::int64_t scale = std::int64_t{1} << 40U;
    constexpr long double kLimit =
        static_cast<long double>(std::int64_t{1} << 60U);
    while (scale > 1 && magnitude * static_cast<long double>(scale) > kLimit) {
        scale /= 2;
    }

    const auto rounded = [scale](double value) {
        const long double scaled =
            static_cast<long double>(value) * static_cast<long double>(scale);
        if (scaled >
                static_cast<long double>(std::numeric_limits<std::int64_t>::max()) ||
            scaled <
                static_cast<long double>(std::numeric_limits<std::int64_t>::min())) {
            throw std::overflow_error("dual scaling overflow");
        }
        return static_cast<std::int64_t>(std::llround(scaled));
    };

    ScaledDuals result;
    result.scale = scale;
    result.item.reserve(duals.item.size());
    for (const double value : duals.item) {
        result.item.push_back(rounded(value));
    }
    result.arc.reserve(duals.arc.size());
    for (const double value : duals.arc) {
        result.arc.push_back(std::max<std::int64_t>(0, rounded(value)));
    }
    result.sr.reserve(duals.sr.size());
    for (const double value : duals.sr) {
        result.sr.push_back(std::min<std::int64_t>(0, rounded(value)));
    }
    return result;
}

[[nodiscard]] int ceil_scaled_bound(std::int64_t numerator,
                                    std::int64_t scale) {
    return exact_arithmetic::ceil_ratio_to_int(
        numerator, scale, "certified lower bound does not fit in int");
}

struct SafeBoundResult {
    bool proven = true;
    int integer_bound = 0;
};

[[nodiscard]] SafeBoundResult compute_safe_lagrangian_bound(
    const Instance& instance,
    const CompiledNode& compiled,
    const RestrictedMaster& master,
    const RestrictedMaster::Duals& floating_duals,
    Clock::time_point end_time,
    Statistics& statistics) {
    SafeBoundResult result;
    const ScaledDuals duals = scale_duals(instance, master.bin_count(),
                                          floating_duals);

    std::vector<std::int64_t> arc_potential(
        static_cast<std::size_t>(instance.size()), 0);
    std::int64_t numerator = 0;
    for (int item = 0; item < instance.size(); ++item) {
        numerator = checked_add(
            numerator, duals.item[static_cast<std::size_t>(item)]);
    }
    for (std::size_t index = 0; index < instance.arcs.size(); ++index) {
        const std::int64_t value = duals.arc[index];
        if (value == 0) {
            continue;
        }
        const Arc& arc = instance.arcs[index];
        arc_potential[static_cast<std::size_t>(arc.to)] = checked_add(
            arc_potential[static_cast<std::size_t>(arc.to)], value);
        arc_potential[static_cast<std::size_t>(arc.from)] = checked_add(
            arc_potential[static_cast<std::size_t>(arc.from)], -value);
        numerator = checked_add(
            numerator,
            exact_arithmetic::checked_multiply(
                arc.separation, value,
                "fixed-point arc certificate overflow"));
    }
    for (const std::int64_t value : duals.sr) {
        numerator = checked_add(numerator, value);
    }

    std::int64_t prefix_value = 0;
    std::int64_t best_prefix_value = 0;
    const auto pricing_start = Clock::now();
    for (int bin = 0; bin < master.bin_count(); ++bin) {
        std::vector<std::int64_t> item_profit = duals.item;
        const std::int64_t multiplier = static_cast<std::int64_t>(bin + 1);
        for (int item = 0; item < instance.size(); ++item) {
            const std::int64_t contribution =
                exact_arithmetic::checked_multiply(
                    multiplier,
                    arc_potential[static_cast<std::size_t>(item)],
                    "fixed-point item profit overflow");
            const std::int64_t total = exact_arithmetic::checked_add(
                item_profit[static_cast<std::size_t>(item)], contribution,
                "fixed-point item profit overflow");
            item_profit[static_cast<std::size_t>(item)] =
                total;
        }
        IntegerPricingSearch pricing(instance, compiled, master.sr_cuts(),
                                     duals.sr, bin, std::move(item_profit),
                                     end_time);
        const IntegerPricingResult priced = pricing.solve();
        ++statistics.pricing_count;
        statistics.pricing_search_nodes += priced.search_nodes;
        if (!priced.proven) {
            result.proven = false;
            break;
        }
        prefix_value = checked_add(
            prefix_value,
            exact_arithmetic::checked_subtract(
                duals.scale, priced.maximum_profit,
                "fixed-point prefix certificate overflow"));
        best_prefix_value = std::min(best_prefix_value, prefix_value);
    }
    statistics.pricing_seconds +=
        std::chrono::duration<double>(Clock::now() - pricing_start).count();
    if (!result.proven) {
        return result;
    }
    numerator = checked_add(numerator, best_prefix_value);
    result.integer_bound = ceil_scaled_bound(numerator, duals.scale);
    return result;
}

[[nodiscard]] std::vector<double> floating_arc_potential(
    const Instance& instance,
    const RestrictedMaster::Duals& duals) {
    std::vector<double> potential(static_cast<std::size_t>(instance.size()), 0.0);
    for (std::size_t index = 0; index < instance.arcs.size(); ++index) {
        const double value = duals.arc[index];
        if (value == 0.0) {
            continue;
        }
        const Arc& arc = instance.arcs[index];
        potential[static_cast<std::size_t>(arc.to)] += value;
        potential[static_cast<std::size_t>(arc.from)] -= value;
    }
    return potential;
}

struct PricingRoundResult {
    bool proven = true;
    bool added_columns = false;
};

[[nodiscard]] PricingRoundResult run_pricing_round(
    const Instance& instance,
    const CompiledNode& compiled,
    RestrictedMaster& master,
    const RestrictedMaster::Duals& duals,
    const Config& config,
    Clock::time_point end_time,
    Statistics& statistics) {
    PricingRoundResult result;
    const auto start = Clock::now();
    const std::vector<double> arc_potential =
        floating_arc_potential(instance, duals);
    const double tolerance = master.phase_one()
                                 ? std::min(1e-10,
                                            config.reduced_cost_tolerance)
                                 : config.reduced_cost_tolerance;
    for (int bin = 0; bin < master.bin_count(); ++bin) {
        if (Clock::now() >= end_time) {
            result.proven = false;
            break;
        }
        const double base_cost = -duals.usage[static_cast<std::size_t>(bin)];
        std::vector<double> item_profit = duals.item;
        const double multiplier = static_cast<double>(bin + 1);
        for (int item = 0; item < instance.size(); ++item) {
            item_profit[static_cast<std::size_t>(item)] +=
                multiplier * arc_potential[static_cast<std::size_t>(item)];
        }
        FloatingPricingSearch pricing(
            instance, compiled, master.sr_cuts(), duals.sr,
            master.pattern_keys(), bin, std::move(item_profit), base_cost, tolerance,
            config.max_columns_per_pricing, end_time);
        FloatingPricingResult priced = pricing.solve();
        ++statistics.pricing_count;
        statistics.pricing_search_nodes += priced.search_nodes;
        if (!priced.proven) {
            result.proven = false;
            break;
        }
        for (FloatingPricedPattern& generated : priced.patterns) {
            master.add_pattern(std::move(generated.pattern));
            ++statistics.generated_columns;
            result.added_columns = true;
        }
    }
    statistics.pricing_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

enum class NodeState {
    kConverged,
    kInfeasible,
    kTimedOut,
};

struct NodeSolveResult {
    NodeState state = NodeState::kTimedOut;
    int certified_lower_bound = 0;
    double lp_value = 0.0;
    std::shared_ptr<const std::vector<Pattern>> pattern_pool;
    std::shared_ptr<const std::vector<int>> active_arcs;
    std::shared_ptr<const std::vector<std::array<int, 3>>> sr_cuts;
    std::vector<double> variable_values;
    std::uint64_t column_count = 0;
};

[[nodiscard]] NodeSolveResult solve_node(
    GRBEnv& environment,
    const Instance& instance,
    const CompiledNode& compiled,
    const SearchNode& node,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    NodeSolveResult result;
    result.certified_lower_bound = node.lower_bound;
    const auto start = Clock::now();
    RestrictedMaster master(environment, instance, compiled, node, config,
                            statistics);
    int sr_rounds = 0;

    for (int iteration = 0; iteration < config.max_cg_iterations; ++iteration) {
        ++statistics.cg_iterations;
        if (deadline.expired()) {
            result.state = NodeState::kTimedOut;
            break;
        }
        if (!master.solve(deadline.end_time(), statistics)) {
            result.state = NodeState::kTimedOut;
            break;
        }

        if (master.phase_one()) {
            if (master.artificial_value() <= kPhaseOneTolerance) {
                master.finish_phase_one();
                continue;
            }
            const RestrictedMaster::Duals duals = master.duals();
            const PricingRoundResult pricing = run_pricing_round(
                instance, compiled, master, duals, config,
                deadline.end_time(), statistics);
            if (!pricing.proven) {
                result.state = NodeState::kTimedOut;
                break;
            }
            if (pricing.added_columns) {
                continue;
            }
            result.state = NodeState::kInfeasible;
            break;
        }

        const std::vector<int> violated_arcs = master.violated_inactive_arcs(
            config.row_violation_tolerance,
            config.max_precedence_rows_per_round);
        if (!violated_arcs.empty()) {
            statistics.generated_precedence_rows += violated_arcs.size();
            master.add_arc_rows(violated_arcs, statistics);
            continue;
        }

        const RestrictedMaster::Duals duals = master.duals();
        const PricingRoundResult pricing = run_pricing_round(
            instance, compiled, master, duals, config,
            deadline.end_time(), statistics);
        if (!pricing.proven) {
            result.state = NodeState::kTimedOut;
            break;
        }
        if (pricing.added_columns) {
            continue;
        }

        if (!config.set_covering_master && config.enable_sr_cuts &&
            sr_rounds < config.max_sr_rounds_per_node &&
            master.sr_cut_count() < config.max_sr_cuts) {
            const int maximum_new = std::min(
                20, config.max_sr_cuts - master.sr_cut_count());
            const auto cuts = master.violated_sr_cuts(
                maximum_new, config.sr_enumeration_item_limit,
                config.row_violation_tolerance);
            if (!cuts.empty()) {
                ++sr_rounds;
                statistics.generated_sr_rows += cuts.size();
                master.add_sr_cuts(cuts, statistics);
                continue;
            }
        }

        result.lp_value = master.objective_value();
        const SafeBoundResult safe = compute_safe_lagrangian_bound(
            instance, compiled, master, duals, deadline.end_time(), statistics);
        if (!safe.proven) {
            result.state = NodeState::kTimedOut;
            break;
        }
        result.certified_lower_bound =
            std::max(node.lower_bound, safe.integer_bound);
        result.pattern_pool =
            std::make_shared<const std::vector<Pattern>>(master.patterns());
        result.active_arcs =
            std::make_shared<const std::vector<int>>(master.active_arcs());
        result.sr_cuts =
            std::make_shared<const std::vector<std::array<int, 3>>>(
                master.sr_cuts());
        result.variable_values = master.variable_values();
        result.state = NodeState::kConverged;
        break;
    }

    statistics.cg_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    ++statistics.cg_count;
    result.column_count = master.patterns().size();
    return result;
}

[[nodiscard]] std::optional<Assignment> recover_integral_assignment(
    const Instance& instance,
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    const int n = instance.size();
    const int m = compiled.bin_count;
    std::vector<double> item_bin(static_cast<std::size_t>(n) * m, 0.0);
    for (std::size_t column = 0; column < patterns.size(); ++column) {
        const double value = values[column];
        if (value <= kEpsilon) {
            continue;
        }
        const Pattern& pattern = patterns[column];
        for (const int item : pattern.items) {
            item_bin[static_cast<std::size_t>(item) * m +
                     static_cast<std::size_t>(pattern.bin)] += value;
        }
    }

    Assignment assignment;
    assignment.bin_of_item.assign(static_cast<std::size_t>(n), -1);
    int last_bin = -1;
    for (int item = 0; item < n; ++item) {
        int selected_bin = -1;
        for (int bin = 0; bin < m; ++bin) {
            const double value =
                item_bin[static_cast<std::size_t>(item) * m +
                         static_cast<std::size_t>(bin)];
            if (std::abs(value - std::round(value)) > kIntegralityTolerance) {
                return std::nullopt;
            }
            if (value > 0.5) {
                if (selected_bin >= 0) {
                    return std::nullopt;
                }
                selected_bin = bin;
            }
        }
        if (selected_bin < 0) {
            return std::nullopt;
        }
        assignment.bin_of_item[static_cast<std::size_t>(item)] = selected_bin;
        last_bin = std::max(last_bin, selected_bin);
    }
    assignment.bin_count = last_bin + 1;

    for (const auto& members : compiled.members) {
        const int bin = assignment.bin_of_item[static_cast<std::size_t>(members[0])];
        for (const int item : members) {
            if (assignment.bin_of_item[static_cast<std::size_t>(item)] != bin) {
                return std::nullopt;
            }
        }
    }
    for (int lhs = 0; lhs < compiled.component_count(); ++lhs) {
        const int lhs_bin = assignment.bin_of_item[static_cast<std::size_t>(
            compiled.members[static_cast<std::size_t>(lhs)][0])];
        if (!compiled.is_eligible(lhs, lhs_bin)) {
            return std::nullopt;
        }
        for (const int rhs :
             compiled.conflict_neighbors[static_cast<std::size_t>(lhs)]) {
            if (rhs > lhs &&
                assignment.bin_of_item[static_cast<std::size_t>(
                    compiled.members[static_cast<std::size_t>(rhs)][0])] ==
                    lhs_bin) {
                return std::nullopt;
            }
        }
    }
    std::string diagnostic;
    if (!check_assignment(instance, assignment, &diagnostic)) {
        throw std::logic_error("integral master solution is invalid: " + diagnostic);
    }
    return assignment;
}

enum class BranchKind {
    kRyanFoster,
    kPosition,
};

struct BranchDecision {
    BranchKind kind = BranchKind::kRyanFoster;
    int first = -1;
    int second = -1;
};

[[nodiscard]] std::optional<BranchDecision> choose_branch(
    const Instance& instance,
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    const int component_count = compiled.component_count();
    std::vector<double> affinity(
        static_cast<std::size_t>(component_count) * component_count, 0.0);
    std::vector<int> selected_components;
    std::vector<unsigned char> seen(static_cast<std::size_t>(component_count), 0U);
    for (std::size_t column = 0; column < patterns.size(); ++column) {
        const double value = values[column];
        if (value <= kEpsilon) {
            continue;
        }
        selected_components.clear();
        for (const int item : patterns[column].items) {
            const int component =
                compiled.component_of_item[static_cast<std::size_t>(item)];
            if (seen[static_cast<std::size_t>(component)] == 0U) {
                seen[static_cast<std::size_t>(component)] = 1U;
                selected_components.push_back(component);
            }
        }
        for (std::size_t i = 0; i < selected_components.size(); ++i) {
            for (std::size_t j = i + 1; j < selected_components.size(); ++j) {
                const int lhs = selected_components[i];
                const int rhs = selected_components[j];
                affinity[static_cast<std::size_t>(lhs) * component_count +
                         static_cast<std::size_t>(rhs)] += value;
                affinity[static_cast<std::size_t>(rhs) * component_count +
                         static_cast<std::size_t>(lhs)] += value;
            }
        }
        for (const int component : selected_components) {
            seen[static_cast<std::size_t>(component)] = 0U;
        }
    }

    double best_fractionality = kIntegralityTolerance;
    int best_weight = -1;
    int best_lhs = -1;
    int best_rhs = -1;
    for (int lhs = 0; lhs < component_count; ++lhs) {
        for (int rhs = lhs + 1; rhs < component_count; ++rhs) {
            if (compiled.conflicts(lhs, rhs)) {
                continue;
            }
            const double value =
                affinity[static_cast<std::size_t>(lhs) * component_count +
                         static_cast<std::size_t>(rhs)];
            const double fractionality = std::abs(value - std::round(value));
            const int weight = compiled.weights[static_cast<std::size_t>(lhs)] +
                               compiled.weights[static_cast<std::size_t>(rhs)];
            if (fractionality > best_fractionality + 1e-12 ||
                (std::abs(fractionality - best_fractionality) <= 1e-12 &&
                 weight > best_weight)) {
                best_fractionality = fractionality;
                best_weight = weight;
                best_lhs = lhs;
                best_rhs = rhs;
            }
        }
    }
    if (best_lhs >= 0) {
        return BranchDecision{
            BranchKind::kRyanFoster,
            compiled.members[static_cast<std::size_t>(best_lhs)][0],
            compiled.members[static_cast<std::size_t>(best_rhs)][0]};
    }

    const int n = instance.size();
    const int m = compiled.bin_count;
    std::vector<double> item_bin(static_cast<std::size_t>(n) * m, 0.0);
    for (std::size_t column = 0; column < patterns.size(); ++column) {
        const double value = values[column];
        if (value <= kEpsilon) {
            continue;
        }
        for (const int item : patterns[column].items) {
            item_bin[static_cast<std::size_t>(item) * m +
                     static_cast<std::size_t>(patterns[column].bin)] += value;
        }
    }
    best_fractionality = kIntegralityTolerance;
    int best_item = -1;
    int best_bin = -1;
    for (int item = 0; item < n; ++item) {
        for (int bin = 0; bin < m; ++bin) {
            const double value =
                item_bin[static_cast<std::size_t>(item) * m +
                         static_cast<std::size_t>(bin)];
            const double fractionality = std::abs(value - std::round(value));
            if (fractionality > best_fractionality + 1e-12 ||
                (std::abs(fractionality - best_fractionality) <= 1e-12 &&
                 (best_item < 0 ||
                  instance.items[static_cast<std::size_t>(item)].weight >
                      instance.items[static_cast<std::size_t>(best_item)].weight))) {
                best_fractionality = fractionality;
                best_item = item;
                best_bin = bin;
            }
        }
    }
    if (best_item >= 0) {
        return BranchDecision{BranchKind::kPosition, best_item, best_bin};
    }
    return std::nullopt;
}

[[nodiscard]] std::shared_ptr<const std::vector<Pattern>> initial_patterns(
    const Instance& instance, const Assignment& assignment) {
    const std::size_t block_count =
        (static_cast<std::size_t>(instance.size()) + 63U) / 64U;
    std::vector<Pattern> patterns(static_cast<std::size_t>(assignment.bin_count));
    for (int bin = 0; bin < assignment.bin_count; ++bin) {
        patterns[static_cast<std::size_t>(bin)].bin = bin;
        patterns[static_cast<std::size_t>(bin)].bits.assign(block_count, 0U);
    }
    for (int item = 0; item < instance.size(); ++item) {
        Pattern& pattern = patterns[static_cast<std::size_t>(
            assignment.bin_of_item[static_cast<std::size_t>(item)])];
        pattern.items.push_back(item);
        pattern.bits[static_cast<std::size_t>(item) / 64U] |=
            std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
    }
    patterns.erase(
        std::remove_if(patterns.begin(), patterns.end(),
                       [](const Pattern& pattern) { return pattern.items.empty(); }),
        patterns.end());
    return std::make_shared<const std::vector<Pattern>>(std::move(patterns));
}

[[nodiscard]] bool m_solution_is_integral(
    const std::vector<double>& values) noexcept {
    for (const double value : values) {
        if (std::abs(value - std::round(value)) >
            kIntegralityTolerance) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] std::optional<int> recover_m_integral_bppc_value(
    const Instance& instance,
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    if (!m_solution_is_integral(values)) {
        return std::nullopt;
    }
    std::vector<int> selected_columns;
    selected_columns.reserve(values.size());
    for (std::size_t column = 0; column < values.size(); ++column) {
        const long long rounded = std::llround(values[column]);
        if (rounded < 0 || rounded > 1) {
            return std::nullopt;
        }
        if (rounded == 1) {
            selected_columns.push_back(static_cast<int>(column));
        }
    }
    if (selected_columns.empty()) {
        return std::nullopt;
    }

    const int component_count = compiled.component_count();
    std::vector<int> component_bin(static_cast<std::size_t>(component_count), -1);
    std::vector<std::int64_t> bin_load(selected_columns.size(), 0);
    std::vector<unsigned char> used(selected_columns.size(), 0U);
    for (int component = 0; component < component_count; ++component) {
        const int representative =
            compiled.members[static_cast<std::size_t>(component)][0];
        for (std::size_t bin = 0; bin < selected_columns.size(); ++bin) {
            const Pattern& pattern = patterns[static_cast<std::size_t>(
                selected_columns[bin])];
            if (pattern.contains(representative)) {
                component_bin[static_cast<std::size_t>(component)] =
                    static_cast<int>(bin);
                bin_load[bin] +=
                    compiled.weights[static_cast<std::size_t>(component)];
                used[bin] = 1U;
                break;
            }
        }
        if (component_bin[static_cast<std::size_t>(component)] < 0) {
            return std::nullopt;
        }
    }
    for (const std::int64_t load : bin_load) {
        if (load > instance.capacity) {
            return std::nullopt;
        }
    }
    for (int lhs = 0; lhs < component_count; ++lhs) {
        for (const int rhs :
             compiled.conflict_neighbors[static_cast<std::size_t>(lhs)]) {
            if (rhs > lhs &&
                component_bin[static_cast<std::size_t>(lhs)] ==
                    component_bin[static_cast<std::size_t>(rhs)]) {
                return std::nullopt;
            }
        }
    }
    return static_cast<int>(
        std::count(used.begin(), used.end(), static_cast<unsigned char>(1U)));
}

[[nodiscard]] std::vector<int> pattern_components(
    const Pattern& pattern,
    const CompiledNode& compiled,
    std::vector<unsigned char>& seen) {
    std::vector<int> result;
    result.reserve(pattern.items.size());
    for (const int item : pattern.items) {
        const int component =
            compiled.component_of_item[static_cast<std::size_t>(item)];
        if (seen[static_cast<std::size_t>(component)] == 0U) {
            seen[static_cast<std::size_t>(component)] = 1U;
            result.push_back(component);
        }
    }
    for (const int component : result) {
        seen[static_cast<std::size_t>(component)] = 0U;
    }
    return result;
}

[[nodiscard]] std::optional<BranchDecision> choose_m_fractional_branch(
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    const int component_count = compiled.component_count();
    std::vector<double> affinity(
        static_cast<std::size_t>(component_count) * component_count, 0.0);
    std::vector<unsigned char> seen(static_cast<std::size_t>(component_count),
                                    0U);
    std::vector<int> components;
    for (std::size_t column = 0; column < patterns.size(); ++column) {
        const double value = values[column];
        if (value <= kEpsilon) {
            continue;
        }
        components = pattern_components(patterns[column], compiled, seen);
        for (std::size_t lhs_pos = 0; lhs_pos < components.size(); ++lhs_pos) {
            for (std::size_t rhs_pos = lhs_pos + 1;
                 rhs_pos < components.size(); ++rhs_pos) {
                const int lhs = components[lhs_pos];
                const int rhs = components[rhs_pos];
                affinity[static_cast<std::size_t>(lhs) * component_count +
                         static_cast<std::size_t>(rhs)] += value;
                affinity[static_cast<std::size_t>(rhs) * component_count +
                         static_cast<std::size_t>(lhs)] += value;
            }
        }
    }

    double best_fractionality = kIntegralityTolerance;
    int best_weight = -1;
    int best_lhs = -1;
    int best_rhs = -1;
    for (int lhs = 0; lhs < component_count; ++lhs) {
        for (int rhs = lhs + 1; rhs < component_count; ++rhs) {
            if (compiled.conflicts(lhs, rhs)) {
                continue;
            }
            const double value =
                affinity[static_cast<std::size_t>(lhs) * component_count +
                         static_cast<std::size_t>(rhs)];
            const double fractionality =
                std::abs(value - std::round(value));
            const int weight = compiled.weights[static_cast<std::size_t>(lhs)] +
                               compiled.weights[static_cast<std::size_t>(rhs)];
            if (fractionality > best_fractionality + 1e-12 ||
                (std::abs(fractionality - best_fractionality) <= 1e-12 &&
                 weight > best_weight)) {
                best_fractionality = fractionality;
                best_weight = weight;
                best_lhs = lhs;
                best_rhs = rhs;
            }
        }
    }
    if (best_lhs >= 0) {
        return BranchDecision{
            BranchKind::kRyanFoster,
            compiled.members[static_cast<std::size_t>(best_lhs)][0],
            compiled.members[static_cast<std::size_t>(best_rhs)][0]};
    }

    for (std::size_t column = 0; column < patterns.size(); ++column) {
        if (std::abs(values[column] - std::round(values[column])) <=
            kIntegralityTolerance) {
            continue;
        }
        components = pattern_components(patterns[column], compiled, seen);
        if (components.size() >= 2U) {
            return BranchDecision{
                BranchKind::kRyanFoster,
                compiled.members[static_cast<std::size_t>(components[0])][0],
                compiled.members[static_cast<std::size_t>(components[1])][0]};
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<BranchDecision> choose_m_integral_refinement(
    const Instance& instance,
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    std::vector<unsigned char> seen(
        static_cast<std::size_t>(compiled.component_count()), 0U);
    for (std::size_t column = 0; column < patterns.size(); ++column) {
        if (values[column] <= 0.5) {
            continue;
        }
        const std::vector<int> components =
            pattern_components(patterns[column], compiled, seen);
        if (components.size() >= 2U) {
            return BranchDecision{
                BranchKind::kRyanFoster,
                compiled.members[static_cast<std::size_t>(components[0])][0],
                compiled.members[static_cast<std::size_t>(components[1])][0]};
        }
    }

    const int component_count = compiled.component_count();
    int best_lhs = -1;
    int best_rhs = -1;
    int best_weight = -1;
    for (int lhs = 0; lhs < component_count; ++lhs) {
        for (int rhs = lhs + 1; rhs < component_count; ++rhs) {
            const int weight =
                compiled.weights[static_cast<std::size_t>(lhs)] +
                compiled.weights[static_cast<std::size_t>(rhs)];
            if (compiled.conflicts(lhs, rhs) || weight > instance.capacity) {
                continue;
            }
            if (weight > best_weight) {
                best_lhs = lhs;
                best_rhs = rhs;
                best_weight = weight;
            }
        }
    }
    if (best_lhs >= 0) {
        return BranchDecision{
            BranchKind::kRyanFoster,
            compiled.members[static_cast<std::size_t>(best_lhs)][0],
            compiled.members[static_cast<std::size_t>(best_rhs)][0]};
    }
    return std::nullopt;
}

[[nodiscard]] std::vector<std::pair<int, int>>
m_integral_refinement_forest(
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    DisjointSet forest(compiled.component_count());
    std::vector<unsigned char> seen(
        static_cast<std::size_t>(compiled.component_count()), 0U);
    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(static_cast<std::size_t>(
        std::max(0, compiled.component_count() - 1)));
    for (std::size_t column = 0; column < patterns.size(); ++column) {
        if (values[column] <= 0.5) {
            continue;
        }
        const std::vector<int> components =
            pattern_components(patterns[column], compiled, seen);
        if (components.empty()) {
            continue;
        }
        const int anchor = components[0];
        for (std::size_t index = 1; index < components.size(); ++index) {
            const int component = components[index];
            if (forest.find(anchor) == forest.find(component)) {
                continue;
            }
            forest.unite(anchor, component);
            pairs.emplace_back(
                compiled.members[static_cast<std::size_t>(anchor)][0],
                compiled.members[static_cast<std::size_t>(component)][0]);
        }
    }
    return pairs;
}

struct MPrecedenceCheckResult {
    bool proven = true;
    std::optional<Assignment> assignment;
};

struct MPrecedenceCheckCacheEntry {
    bool infeasible = false;
    std::optional<Assignment> optimal_assignment;
};

struct MPrecedenceCheckCache {
    std::unordered_map<std::string, MPrecedenceCheckCacheEntry> entries;
};

[[nodiscard]] std::optional<std::string> exact_cover_cache_key(
    const Instance& instance,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) {
    std::vector<unsigned char> coverage(
        static_cast<std::size_t>(instance.size()), 0U);
    std::vector<std::string> selected_bits;
    for (std::size_t column = 0; column < values.size(); ++column) {
        if (values[column] <= 0.5) {
            continue;
        }
        for (const int item : patterns[column].items) {
            unsigned char& count =
                coverage[static_cast<std::size_t>(item)];
            if (++count > 1U) {
                return std::nullopt;
            }
        }
        const auto& bits = patterns[column].bits;
        selected_bits.emplace_back(
            reinterpret_cast<const char*>(bits.data()),
            bits.size() * sizeof(std::uint64_t));
    }
    if (std::any_of(coverage.begin(), coverage.end(),
                    [](unsigned char value) { return value != 1U; })) {
        return std::nullopt;
    }
    std::sort(selected_bits.begin(), selected_bits.end());
    const std::uint32_t count =
        static_cast<std::uint32_t>(selected_bits.size());
    std::string key(sizeof(count), '\0');
    std::memcpy(key.data(), &count, sizeof(count));
    for (const std::string& bits : selected_bits) {
        key.append(bits);
    }
    return key;
}

[[nodiscard]] bool m_pair_in_positive_pattern(
    int first,
    int second,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values) noexcept {
    for (std::size_t column = 0; column < values.size(); ++column) {
        if (values[column] > 0.5 && patterns[column].contains(first) &&
            patterns[column].contains(second)) {
            return true;
        }
    }
    return false;
}

struct MFastExactCoverCheck {
    bool decided = false;
    bool infeasible = false;
    bool optimal_assignment = false;
    std::optional<Assignment> assignment;
};

[[nodiscard]] MFastExactCoverCheck fast_exact_cover_check(
    const Instance& instance,
    const std::vector<Pattern>& patterns,
    const std::vector<int>& selected_columns,
    int maximum_horizon) {
    MFastExactCoverCheck result;
    const int group_count = static_cast<int>(selected_columns.size());
    if (group_count == 0) {
        return result;
    }
    if (group_count > maximum_horizon) {
        result.decided = true;
        result.infeasible = true;
        return result;
    }

    std::vector<int> group_of_item(
        static_cast<std::size_t>(instance.size()), -1);
    for (int group = 0; group < group_count; ++group) {
        const Pattern& pattern = patterns[static_cast<std::size_t>(
            selected_columns[static_cast<std::size_t>(group)])];
        for (const int item : pattern.items) {
            int& owner = group_of_item[static_cast<std::size_t>(item)];
            if (owner >= 0) {
                return result;
            }
            owner = group;
        }
    }
    if (std::any_of(group_of_item.begin(), group_of_item.end(),
                    [](int group) { return group < 0; })) {
        return result;
    }

    std::vector<int> lag(
        static_cast<std::size_t>(group_count) * group_count, -1);
    int maximum_lag = 0;
    for (const Arc& arc : instance.arcs) {
        const int from =
            group_of_item[static_cast<std::size_t>(arc.from)];
        const int to = group_of_item[static_cast<std::size_t>(arc.to)];
        if (from == to) {
            if (arc.separation > 0) {
                result.decided = true;
                result.infeasible = true;
                return result;
            }
            continue;
        }
        int& value = lag[static_cast<std::size_t>(from) * group_count +
                         static_cast<std::size_t>(to)];
        value = std::max(value, arc.separation);
        maximum_lag = std::max(maximum_lag, arc.separation);
    }

    std::vector<int> indegree(static_cast<std::size_t>(group_count), 0);
    for (int from = 0; from < group_count; ++from) {
        for (int to = 0; to < group_count; ++to) {
            if (lag[static_cast<std::size_t>(from) * group_count +
                    static_cast<std::size_t>(to)] >= 0) {
                ++indegree[static_cast<std::size_t>(to)];
            }
        }
    }
    std::vector<int> topological_order;
    topological_order.reserve(static_cast<std::size_t>(group_count));
    std::vector<int> queue;
    queue.reserve(static_cast<std::size_t>(group_count));
    for (int group = 0; group < group_count; ++group) {
        if (indegree[static_cast<std::size_t>(group)] == 0) {
            queue.push_back(group);
        }
    }
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const int from = queue[head];
        topological_order.push_back(from);
        for (int to = 0; to < group_count; ++to) {
            if (lag[static_cast<std::size_t>(from) * group_count +
                    static_cast<std::size_t>(to)] < 0) {
                continue;
            }
            if (--indegree[static_cast<std::size_t>(to)] == 0) {
                queue.push_back(to);
            }
        }
    }
    if (static_cast<int>(topological_order.size()) != group_count) {
        result.decided = true;
        result.infeasible = true;
        return result;
    }

    std::vector<int> position(static_cast<std::size_t>(group_count), -1);
    if (maximum_lag <= 1) {
        for (int index = 0; index < group_count; ++index) {
            position[static_cast<std::size_t>(
                topological_order[static_cast<std::size_t>(index)])] = index;
        }
        result.optimal_assignment = true;
    } else {
        std::vector<int> earliest(static_cast<std::size_t>(group_count), 0);
        for (const int from : topological_order) {
            for (int to = 0; to < group_count; ++to) {
                const int value =
                    lag[static_cast<std::size_t>(from) * group_count +
                        static_cast<std::size_t>(to)];
                if (value >= 0) {
                    earliest[static_cast<std::size_t>(to)] = std::max(
                        earliest[static_cast<std::size_t>(to)],
                        earliest[static_cast<std::size_t>(from)] + value);
                }
            }
        }
        const int span_lower_bound =
            1 + *std::max_element(earliest.begin(), earliest.end());
        const int lower_bound = std::max(group_count, span_lower_bound);
        if (lower_bound > maximum_horizon) {
            result.decided = true;
            result.infeasible = true;
            return result;
        }

        std::vector<int> tail(static_cast<std::size_t>(group_count), 0);
        for (auto iterator = topological_order.rbegin();
             iterator != topological_order.rend(); ++iterator) {
            const int from = *iterator;
            for (int to = 0; to < group_count; ++to) {
                const int value =
                    lag[static_cast<std::size_t>(from) * group_count +
                        static_cast<std::size_t>(to)];
                if (value >= 0) {
                    tail[static_cast<std::size_t>(from)] = std::max(
                        tail[static_cast<std::size_t>(from)],
                        value + tail[static_cast<std::size_t>(to)]);
                }
            }
        }

        indegree.assign(static_cast<std::size_t>(group_count), 0);
        for (int from = 0; from < group_count; ++from) {
            for (int to = 0; to < group_count; ++to) {
                if (lag[static_cast<std::size_t>(from) * group_count +
                        static_cast<std::size_t>(to)] >= 0) {
                    ++indegree[static_cast<std::size_t>(to)];
                }
            }
        }
        std::vector<int> release(static_cast<std::size_t>(group_count), 0);
        int scheduled = 0;
        int last_position = -1;
        for (int bin = 0; bin < maximum_horizon && scheduled < group_count;
             ++bin) {
            int selected = -1;
            for (int group = 0; group < group_count; ++group) {
                if (position[static_cast<std::size_t>(group)] >= 0 ||
                    indegree[static_cast<std::size_t>(group)] != 0 ||
                    release[static_cast<std::size_t>(group)] > bin) {
                    continue;
                }
                if (selected < 0 ||
                    tail[static_cast<std::size_t>(group)] >
                        tail[static_cast<std::size_t>(selected)] ||
                    (tail[static_cast<std::size_t>(group)] ==
                         tail[static_cast<std::size_t>(selected)] &&
                     group < selected)) {
                    selected = group;
                }
            }
            if (selected < 0) {
                continue;
            }
            position[static_cast<std::size_t>(selected)] = bin;
            last_position = bin;
            ++scheduled;
            for (int to = 0; to < group_count; ++to) {
                const int value =
                    lag[static_cast<std::size_t>(selected) * group_count +
                        static_cast<std::size_t>(to)];
                if (value < 0) {
                    continue;
                }
                release[static_cast<std::size_t>(to)] = std::max(
                    release[static_cast<std::size_t>(to)], bin + value);
                --indegree[static_cast<std::size_t>(to)];
            }
        }
        if (scheduled != group_count) {
            return result;
        }
        result.optimal_assignment = last_position + 1 == lower_bound;
    }

    Assignment assignment;
    assignment.bin_of_item.assign(static_cast<std::size_t>(instance.size()),
                                  -1);
    int last_bin = -1;
    for (int item = 0; item < instance.size(); ++item) {
        const int bin = position[static_cast<std::size_t>(
            group_of_item[static_cast<std::size_t>(item)])];
        assignment.bin_of_item[static_cast<std::size_t>(item)] = bin;
        last_bin = std::max(last_bin, bin);
    }
    assignment.bin_count = last_bin + 1;
    std::string diagnostic;
    if (!check_assignment(instance, assignment, &diagnostic)) {
        throw std::logic_error(
            "fast M precedence check returned an invalid assignment: " +
            diagnostic);
    }
    result.decided = true;
    result.assignment = std::move(assignment);
    return result;
}

[[nodiscard]] MPrecedenceCheckResult check_m_integral_cover(
    GRBEnv& environment,
    const Instance& instance,
    const CompiledNode& compiled,
    const std::vector<Pattern>& patterns,
    const std::vector<double>& values,
    int incumbent_value,
    const Config& config,
    Clock::time_point end_time,
    Statistics& statistics,
    MPrecedenceCheckCache& cache) {
    MPrecedenceCheckResult result;
    const auto check_start = Clock::now();
    ++statistics.precedence_check_count;
    const auto finish_statistics = [&]() {
        statistics.precedence_check_seconds +=
            std::chrono::duration<double>(Clock::now() - check_start).count();
    };

    int maximum_horizon = incumbent_value - 1;
    if (maximum_horizon <= 0) {
        finish_statistics();
        return result;
    }
    std::vector<int> selected_columns;
    for (std::size_t column = 0; column < values.size(); ++column) {
        if (values[column] > 0.5) {
            selected_columns.push_back(static_cast<int>(column));
        }
    }
    if (selected_columns.empty()) {
        throw std::logic_error("integral M cover selected no pattern");
    }

    const std::optional<std::string> cache_key =
        exact_cover_cache_key(instance, patterns, values);
    if (cache_key.has_value()) {
        const auto found = cache.entries.find(*cache_key);
        if (found != cache.entries.end()) {
            ++statistics.precedence_check_cache_hits;
            if (found->second.optimal_assignment.has_value() &&
                found->second.optimal_assignment->bin_count <
                    incumbent_value) {
                result.assignment = found->second.optimal_assignment;
            }
            result.proven = true;
            finish_statistics();
            return result;
        }

        const MFastExactCoverCheck fast = fast_exact_cover_check(
            instance, patterns, selected_columns, maximum_horizon);
        if (fast.decided) {
            result.assignment = fast.assignment;
            if (fast.infeasible || fast.optimal_assignment) {
                result.proven = true;
                MPrecedenceCheckCacheEntry entry;
                entry.infeasible = fast.infeasible;
                entry.optimal_assignment = fast.assignment;
                cache.entries.emplace(*cache_key, std::move(entry));
                finish_statistics();
                return result;
            }
            maximum_horizon = fast.assignment->bin_count - 1;
        }
    }

    const double remaining =
        std::chrono::duration<double>(end_time - Clock::now()).count();
    if (remaining <= 0.0) {
        result.proven = false;
        finish_statistics();
        return result;
    }

    GRBModel model(environment);
    model.set(GRB_IntParam_Threads, 1);
    model.set(GRB_IntParam_Seed, config.seed);
    model.set(GRB_IntParam_OutputFlag, config.gurobi_log ? 1 : 0);
    model.set(GRB_IntParam_NumericFocus, 2);
    model.set(GRB_DoubleParam_MIPGap, 0.0);
    model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
    model.set(GRB_DoubleParam_IntFeasTol, 1e-9);
    model.set(GRB_DoubleParam_TimeLimit, remaining);

    GRBVar horizon = model.addVar(1.0, static_cast<double>(maximum_horizon),
                                  1.0, GRB_INTEGER);
    const int selected_count = static_cast<int>(selected_columns.size());
    std::vector<GRBVar> placement;
    placement.reserve(static_cast<std::size_t>(selected_count) *
                      maximum_horizon);
    for (int pattern = 0; pattern < selected_count; ++pattern) {
        for (int bin = 0; bin < maximum_horizon; ++bin) {
            placement.push_back(
                model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
        }
    }
    model.update();

    const int component_count = compiled.component_count();
    std::vector<GRBLinExpr> component_assignment(
        static_cast<std::size_t>(component_count), GRBLinExpr(0.0));
    std::vector<GRBLinExpr> component_position(
        static_cast<std::size_t>(component_count), GRBLinExpr(0.0));
    std::vector<GRBLinExpr> slot_use(
        static_cast<std::size_t>(selected_count) * maximum_horizon,
        GRBLinExpr(0.0));

    for (int component = 0; component < component_count; ++component) {
        const int representative =
            compiled.members[static_cast<std::size_t>(component)][0];
        bool has_pattern = false;
        for (int selected = 0; selected < selected_count; ++selected) {
            const Pattern& pattern = patterns[static_cast<std::size_t>(
                selected_columns[static_cast<std::size_t>(selected)])];
            if (!pattern.contains(representative)) {
                continue;
            }
            has_pattern = true;
            for (int bin = 0; bin < maximum_horizon; ++bin) {
                GRBVar assigned =
                    model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                component_assignment[static_cast<std::size_t>(component)] +=
                    assigned;
                component_position[static_cast<std::size_t>(component)] +=
                    static_cast<double>(bin) * assigned;
                const std::size_t slot =
                    static_cast<std::size_t>(selected) * maximum_horizon +
                    static_cast<std::size_t>(bin);
                slot_use[slot] += assigned;
                model.addConstr(assigned <= placement[slot]);
            }
        }
        if (!has_pattern) {
            throw std::logic_error(
                "integral M cover omitted a together component");
        }
    }
    model.update();

    for (int component = 0; component < component_count; ++component) {
        model.addConstr(
            component_assignment[static_cast<std::size_t>(component)] == 1.0);
    }
    for (int selected = 0; selected < selected_count; ++selected) {
        GRBLinExpr positions = 0.0;
        for (int bin = 0; bin < maximum_horizon; ++bin) {
            const std::size_t slot =
                static_cast<std::size_t>(selected) * maximum_horizon +
                static_cast<std::size_t>(bin);
            positions += placement[slot];
            model.addConstr(slot_use[slot] <=
                            static_cast<double>(component_count) *
                                placement[slot]);
            model.addConstr(placement[slot] <= slot_use[slot]);
            model.addConstr(horizon >=
                            static_cast<double>(bin + 1) * placement[slot]);
        }
        model.addConstr(positions <= 1.0);
    }
    for (int bin = 0; bin < maximum_horizon; ++bin) {
        GRBLinExpr one_pattern = 0.0;
        for (int selected = 0; selected < selected_count; ++selected) {
            one_pattern += placement[
                static_cast<std::size_t>(selected) * maximum_horizon +
                static_cast<std::size_t>(bin)];
        }
        model.addConstr(one_pattern <= 1.0);
    }
    for (const Arc& arc : instance.arcs) {
        const int from = compiled.component_of_item[
            static_cast<std::size_t>(arc.from)];
        const int to = compiled.component_of_item[
            static_cast<std::size_t>(arc.to)];
        if (from == to) {
            if (arc.separation > 0) {
                throw std::logic_error(
                    "positive precedence survived inside a together component");
            }
            continue;
        }
        model.addConstr(
            component_position[static_cast<std::size_t>(to)] -
                component_position[static_cast<std::size_t>(from)] >=
            static_cast<double>(arc.separation));
    }
    model.setObjective(GRBLinExpr(horizon), GRB_MINIMIZE);
    model.optimize();

    const int status = model.get(GRB_IntAttr_Status);
    const int solution_count = model.get(GRB_IntAttr_SolCount);
    if (solution_count > 0) {
        Assignment assignment;
        assignment.bin_of_item.assign(
            static_cast<std::size_t>(instance.size()), -1);
        int last_bin = -1;
        for (int component = 0; component < component_count; ++component) {
            const int bin = static_cast<int>(std::llround(
                component_position[static_cast<std::size_t>(component)]
                    .getValue()));
            for (const int item :
                 compiled.members[static_cast<std::size_t>(component)]) {
                assignment.bin_of_item[static_cast<std::size_t>(item)] = bin;
            }
            last_bin = std::max(last_bin, bin);
        }
        assignment.bin_count = last_bin + 1;
        std::string diagnostic;
        if (!check_assignment(instance, assignment, &diagnostic)) {
            throw std::logic_error(
                "M precedence check returned an invalid assignment: " +
                diagnostic);
        }
        result.assignment = std::move(assignment);
    }
    result.proven = status == GRB_OPTIMAL || status == GRB_INFEASIBLE;
    if (result.proven && cache_key.has_value()) {
        MPrecedenceCheckCacheEntry entry;
        entry.infeasible = status == GRB_INFEASIBLE &&
                           !result.assignment.has_value();
        if (result.assignment.has_value()) {
            entry.optimal_assignment = result.assignment;
        }
        cache.entries.emplace(*cache_key, std::move(entry));
    }
    finish_statistics();
    return result;
}

[[nodiscard]] int open_tree_lower_bound(const std::vector<SearchNode>& stack,
                                        int incumbent,
                                        std::optional<int> current = std::nullopt) {
    int value = incumbent;
    if (current.has_value() && *current < incumbent) {
        value = std::min(value, *current);
    }
    for (const SearchNode& node : stack) {
        if (node.lower_bound < incumbent) {
            value = std::min(value, node.lower_bound);
        }
    }
    return value;
}

}

BranchPriceResult run_branch_price_and_cut(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    if (config.set_covering_master && !config.root_node_only) {
        throw std::invalid_argument(
            "the direct set-covering master is available only as a root bound");
    }
    BranchPriceResult result;
    result.attempted = true;
    result.incumbent = incumbent;
    result.certified_lower_bound = lower_bound;
    result.root_integer_lower_bound = lower_bound;
    result.root_lp_value = std::numeric_limits<double>::quiet_NaN();
    if (deadline.expired() || lower_bound >= incumbent.bin_count) {
        result.optimal = lower_bound >= incumbent.bin_count;
        result.timed_out = deadline.expired() && !result.optimal;
        return result;
    }

    const double tree_budget = config.cg_time_limit_seconds > 0.0
                                   ? std::min(config.cg_time_limit_seconds,
                                              deadline.remaining_seconds())
                                   : deadline.remaining_seconds();
    Deadline tree_deadline(tree_budget);

    const std::vector<std::pair<int, int>> base_conflicts =
        build_base_conflict_edges(instance);
    auto empty_arcs = std::make_shared<const std::vector<int>>();
    auto empty_cuts =
        std::make_shared<const std::vector<std::array<int, 3>>>();
    SearchNode root;
    root.id = 0;
    root.depth = 0;
    root.lower_bound = lower_bound;
    root.pattern_pool = initial_patterns(instance, incumbent);
    root.active_arcs = std::move(empty_arcs);
    root.sr_cuts = std::move(empty_cuts);

    std::vector<SearchNode> stack;
    stack.push_back(std::move(root));
    std::uint64_t next_node_id = 1;
    bool root_recorded = false;
    std::optional<int> interrupted_node_bound;

    while (!stack.empty() && !tree_deadline.expired()) {
        SearchNode node = std::move(stack.back());
        stack.pop_back();
        if (node.lower_bound >= result.incumbent.bin_count) {
            continue;
        }
        ++statistics.explored_nodes;
        const CompiledNode compiled = compile_node(
            instance, result.incumbent.bin_count, base_conflicts, node);
        if (compiled.infeasible) {
            ++statistics.infeasible_nodes;
            continue;
        }

        NodeSolveResult solved = solve_node(
            environment, instance, compiled, node, config, tree_deadline,
            statistics);
        if (node.id == 0) {
            result.root_column_count = solved.column_count;
        }
        if (solved.state == NodeState::kTimedOut) {
            interrupted_node_bound = solved.certified_lower_bound;
            result.timed_out = true;
            break;
        }
        if (solved.state == NodeState::kInfeasible) {
            ++statistics.infeasible_nodes;
            continue;
        }
        node.lower_bound =
            std::max(node.lower_bound, solved.certified_lower_bound);
        if (!root_recorded && node.id == 0) {
            result.root_lp_value = solved.lp_value;
            result.root_integer_lower_bound = node.lower_bound;
            root_recorded = true;
        }

        if (config.set_covering_master) {
            if (!config.root_node_only || node.id != 0) {
                throw std::logic_error(
                    "direct set-covering master escaped its root-only scope");
            }
            result.root_only_completed = true;
            result.certified_lower_bound = std::min(
                result.incumbent.bin_count, node.lower_bound);
            result.optimal =
                result.certified_lower_bound >= result.incumbent.bin_count;
            return result;
        }

        std::optional<BranchDecision> branch;
        const std::optional<Assignment> integral =
            recover_integral_assignment(
                instance, compiled, *solved.pattern_pool,
                solved.variable_values);
        if (integral.has_value()) {
            ++statistics.integral_nodes;
            if (integral->bin_count < result.incumbent.bin_count) {
                result.incumbent = *integral;
            }
            if (config.root_node_only && node.id == 0) {
                result.root_only_completed = true;
                result.certified_lower_bound = std::min(
                    result.incumbent.bin_count, node.lower_bound);
                result.optimal = result.certified_lower_bound >=
                                 result.incumbent.bin_count;
                return result;
            }
            if (node.lower_bound >= result.incumbent.bin_count) {
                continue;
            }
            interrupted_node_bound = node.lower_bound;
            result.timed_out = true;
            break;
        }
        if (config.root_node_only && node.id == 0) {
            result.root_only_completed = true;
            result.certified_lower_bound = std::min(
                result.incumbent.bin_count, node.lower_bound);
            result.optimal =
                result.certified_lower_bound >= result.incumbent.bin_count;
            return result;
        }
        if (node.lower_bound >= result.incumbent.bin_count) {
            continue;
        }
        branch = choose_branch(
            instance, compiled, *solved.pattern_pool,
            solved.variable_values);
        if (!branch.has_value()) {
            throw std::logic_error(
                "fractional master solution has neither a Ryan-Foster nor an "
                "item-position branch candidate");
        }

        SearchNode left = node;
        SearchNode right = node;
        left.id = next_node_id++;
        right.id = next_node_id++;
        left.depth = node.depth + 1;
        right.depth = node.depth + 1;
        left.pattern_pool = solved.pattern_pool;
        right.pattern_pool = solved.pattern_pool;
        left.active_arcs = solved.active_arcs;
        right.active_arcs = solved.active_arcs;
        left.sr_cuts = solved.sr_cuts;
        right.sr_cuts = solved.sr_cuts;
        if (branch->kind == BranchKind::kRyanFoster) {
            left.together.emplace_back(branch->first, branch->second);
            right.separate.emplace_back(branch->first, branch->second);
            ++statistics.rf_branches;
        } else {
            left.position_decisions.push_back(
                PositionDecision{branch->first, branch->second, true});
            right.position_decisions.push_back(
                PositionDecision{branch->first, branch->second, false});
            ++statistics.position_branches;
        }
        stack.push_back(std::move(right));
        stack.push_back(std::move(left));
    }

    if (stack.empty() && !interrupted_node_bound.has_value()) {
        result.optimal = true;
        result.certified_lower_bound = result.incumbent.bin_count;
    } else {
        result.timed_out = true;
        result.certified_lower_bound = open_tree_lower_bound(
            stack, result.incumbent.bin_count, interrupted_node_bound);
    }
    return result;
}

BranchPriceResult run_m_branch_price_and_cut(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    if (config.root_model != RootModelKind::kM) {
        throw std::invalid_argument(
            "the position-free M tree requires the M root model");
    }
    BranchPriceResult result;
    result.attempted = true;
    result.incumbent = incumbent;
    result.certified_lower_bound = lower_bound;
    result.root_integer_lower_bound = lower_bound;
    result.root_lp_value = std::numeric_limits<double>::quiet_NaN();
    if (deadline.expired() || lower_bound >= incumbent.bin_count) {
        result.optimal = lower_bound >= incumbent.bin_count;
        result.timed_out = deadline.expired() && !result.optimal;
        return result;
    }

    const double tree_budget = config.cg_time_limit_seconds > 0.0
        ? std::min(config.cg_time_limit_seconds,
                   deadline.remaining_seconds())
        : deadline.remaining_seconds();
    Deadline tree_deadline(tree_budget);
    const std::vector<std::pair<int, int>> base_conflicts =
        build_base_conflict_edges(instance);
    auto empty_arcs = std::make_shared<const std::vector<int>>();
    auto empty_cuts =
        std::make_shared<const std::vector<std::array<int, 3>>>();

    SearchNode root;
    root.id = 0;
    root.depth = 0;
    root.lower_bound = lower_bound;
    root.pattern_pool = initial_patterns(instance, incumbent);
    root.active_arcs = std::move(empty_arcs);
    root.sr_cuts = std::move(empty_cuts);

    std::vector<SearchNode> stack;
    stack.push_back(std::move(root));
    std::uint64_t next_node_id = 1;
    bool root_recorded = false;
    std::optional<int> interrupted_node_bound;
    MPrecedenceCheckCache precedence_check_cache;

    while (!stack.empty() && !tree_deadline.expired()) {
        SearchNode node = std::move(stack.back());
        stack.pop_back();
        if (node.lower_bound >= result.incumbent.bin_count) {
            continue;
        }
        ++statistics.explored_nodes;
        const CompiledNode compiled = compile_node(
            instance, 1, base_conflicts, node, false);
        if (compiled.infeasible) {
            ++statistics.infeasible_nodes;
            continue;
        }

        MNodeSolveResult solved = solve_m_node(
            environment, instance, compiled, node,
            result.incumbent.bin_count, config, tree_deadline, statistics);
        if (node.id == 0) {
            result.root_column_count = solved.column_count;
        }
        if (solved.state == MNodeState::kTimedOut) {
            interrupted_node_bound = solved.certified_lower_bound;
            result.timed_out = true;
            break;
        }
        if (solved.state == MNodeState::kInfeasible) {
            ++statistics.infeasible_nodes;
            continue;
        }
        node.lower_bound =
            std::max(node.lower_bound, solved.certified_lower_bound);
        if (!root_recorded && node.id == 0) {
            result.root_lp_value = solved.lp_value;
            result.root_integer_lower_bound = node.lower_bound;
            root_recorded = true;
        }
        if (config.root_node_only && node.id == 0) {
            result.root_only_completed = true;
            result.certified_lower_bound = std::min(
                result.incumbent.bin_count, node.lower_bound);
            result.optimal =
                result.certified_lower_bound >= result.incumbent.bin_count;
            return result;
        }
        if (node.lower_bound >= result.incumbent.bin_count) {
            continue;
        }

        std::optional<BranchDecision> branch;
        const bool integral_solution =
            m_solution_is_integral(solved.variable_values);
        bool prefer_together_child = true;
        if (integral_solution) {
            ++statistics.integral_nodes;
            MPrecedenceCheckResult checked = check_m_integral_cover(
                environment, instance, compiled, *solved.pattern_pool,
                solved.variable_values, result.incumbent.bin_count, config,
                tree_deadline.end_time(), statistics,
                precedence_check_cache);
            if (checked.assignment.has_value() &&
                checked.assignment->bin_count < result.incumbent.bin_count) {
                result.incumbent = *checked.assignment;
            }
            if (!checked.proven) {
                interrupted_node_bound = node.lower_bound;
                result.timed_out = true;
                break;
            }
            if (node.lower_bound >= result.incumbent.bin_count) {
                continue;
            }
            const std::vector<std::pair<int, int>> refinement_pairs =
                m_integral_refinement_forest(
                    compiled, *solved.pattern_pool,
                    solved.variable_values);
            if (!refinement_pairs.empty()) {
                std::vector<SearchNode> children;
                children.reserve(refinement_pairs.size() + 1U);
                SearchNode together_prefix = node;
                together_prefix.depth = node.depth + 1;
                together_prefix.pattern_pool = solved.pattern_pool;
                together_prefix.sr_cuts = solved.sr_cuts;
                for (const auto& pair : refinement_pairs) {
                    SearchNode split = together_prefix;
                    split.id = next_node_id++;
                    split.separate.push_back(pair);
                    children.push_back(std::move(split));
                    together_prefix.together.push_back(pair);
                }
                together_prefix.id = next_node_id++;
                children.push_back(std::move(together_prefix));
                statistics.rf_branches += refinement_pairs.size();

                for (auto iterator = children.rbegin();
                     iterator != children.rend(); ++iterator) {
                    stack.push_back(std::move(*iterator));
                }
                continue;
            }
            branch = choose_m_integral_refinement(
                instance, compiled, *solved.pattern_pool,
                solved.variable_values);
            if (!branch.has_value()) {
                continue;
            }
            prefer_together_child = !m_pair_in_positive_pattern(
                branch->first, branch->second, *solved.pattern_pool,
                solved.variable_values);
        } else {
            branch = choose_m_fractional_branch(
                compiled, *solved.pattern_pool, solved.variable_values);
            if (!branch.has_value()) {
                interrupted_node_bound = node.lower_bound;
                result.timed_out = true;
                break;
            }
        }

        SearchNode together = node;
        SearchNode separate = node;
        together.id = next_node_id++;
        separate.id = next_node_id++;
        together.depth = node.depth + 1;
        separate.depth = node.depth + 1;
        together.pattern_pool = solved.pattern_pool;
        separate.pattern_pool = solved.pattern_pool;
        together.sr_cuts = solved.sr_cuts;
        separate.sr_cuts = solved.sr_cuts;
        together.together.emplace_back(branch->first, branch->second);
        separate.separate.emplace_back(branch->first, branch->second);
        ++statistics.rf_branches;
        if (prefer_together_child) {
            stack.push_back(std::move(separate));
            stack.push_back(std::move(together));
        } else {
            stack.push_back(std::move(together));
            stack.push_back(std::move(separate));
        }
    }

    if (stack.empty() && !interrupted_node_bound.has_value()) {
        result.optimal = true;
        result.certified_lower_bound = result.incumbent.bin_count;
    } else {
        result.timed_out = true;
        result.certified_lower_bound = open_tree_lower_bound(
            stack, result.incumbent.bin_count, interrupted_node_bound);
    }
    return result;
}

BppcBoundResult run_bppc_branch_price_bound(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& original_feasible_incumbent,
    int bppc_lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    BppcBoundResult result;
    result.attempted = true;
    result.certified_lower_bound = bppc_lower_bound;
    result.root_integer_lower_bound = bppc_lower_bound;
    result.root_lp_value = std::numeric_limits<double>::quiet_NaN();

    const std::shared_ptr<const std::vector<Pattern>> incumbent_patterns =
        initial_patterns(instance, original_feasible_incumbent);
    result.incumbent_value =
        static_cast<int>(incumbent_patterns->size());
    if (result.incumbent_value <= 0) {
        throw std::invalid_argument(
            "the BPPC relaxation requires a nonempty feasible incumbent");
    }
    const int capacity_lower_bound = static_cast<int>(
        instance.total_weight / instance.capacity +
        (instance.total_weight % instance.capacity != 0 ? 1 : 0));
    bppc_lower_bound = std::max(bppc_lower_bound, capacity_lower_bound);
    if (bppc_lower_bound > result.incumbent_value) {
        throw std::invalid_argument(
            "the supplied BPPC lower bound exceeds its feasible incumbent");
    }
    result.certified_lower_bound = bppc_lower_bound;
    result.root_integer_lower_bound = bppc_lower_bound;
    if (deadline.expired() || bppc_lower_bound >= result.incumbent_value) {
        result.optimal = bppc_lower_bound >= result.incumbent_value;
        result.timed_out = deadline.expired() && !result.optimal;
        return result;
    }

    const double tree_budget = config.cg_time_limit_seconds > 0.0
        ? std::min(config.cg_time_limit_seconds,
                   deadline.remaining_seconds())
        : deadline.remaining_seconds();
    Deadline tree_deadline(tree_budget);
    const std::vector<std::pair<int, int>> base_conflicts =
        build_bppc_relaxation_conflict_edges(instance);
    auto empty_arcs = std::make_shared<const std::vector<int>>();
    auto empty_cuts =
        std::make_shared<const std::vector<std::array<int, 3>>>();

    SearchNode root;
    root.id = 0;
    root.depth = 0;
    root.lower_bound = bppc_lower_bound;
    root.pattern_pool = incumbent_patterns;
    root.active_arcs = std::move(empty_arcs);
    root.sr_cuts = std::move(empty_cuts);

    std::vector<SearchNode> stack;
    stack.push_back(std::move(root));
    std::uint64_t next_node_id = 1;
    bool root_recorded = false;
    std::optional<int> interrupted_node_bound;

    while (!stack.empty() && !tree_deadline.expired()) {
        SearchNode node = std::move(stack.back());
        stack.pop_back();
        if (node.lower_bound >= result.incumbent_value) {
            continue;
        }
        ++statistics.explored_nodes;
        const CompiledNode compiled = compile_node(
            instance, 1, base_conflicts, node, false);
        if (compiled.infeasible) {
            ++statistics.infeasible_nodes;
            continue;
        }

        MNodeSolveResult solved = solve_m_node(
            environment, instance, compiled, node,
            result.incumbent_value, config, tree_deadline, statistics, true);
        if (node.id == 0) {
            result.root_column_count = solved.column_count;
        }
        if (solved.state == MNodeState::kTimedOut) {
            interrupted_node_bound = solved.certified_lower_bound;
            result.timed_out = true;
            break;
        }
        if (solved.state == MNodeState::kInfeasible) {
            ++statistics.infeasible_nodes;
            continue;
        }
        node.lower_bound =
            std::max(node.lower_bound, solved.certified_lower_bound);
        if (!root_recorded && node.id == 0) {
            result.root_lp_value = solved.lp_value;
            result.root_integer_lower_bound = node.lower_bound;
            root_recorded = true;
        }
        if (node.lower_bound >= result.incumbent_value) {
            continue;
        }

        const bool integral_solution =
            m_solution_is_integral(solved.variable_values);
        if (integral_solution) {
            ++statistics.integral_nodes;
            const std::optional<int> integral_value =
                recover_m_integral_bppc_value(
                    instance, compiled, *solved.pattern_pool,
                    solved.variable_values);
            if (!integral_value.has_value()) {
                throw std::logic_error(
                    "integral BPPC set cover could not be converted to a "
                    "partition");
            }
            result.incumbent_value =
                std::min(result.incumbent_value, *integral_value);
            if (node.lower_bound >= result.incumbent_value) {
                continue;
            }
            interrupted_node_bound = node.lower_bound;
            result.timed_out = true;
            break;
        }

        const std::optional<BranchDecision> branch =
            choose_m_fractional_branch(
                compiled, *solved.pattern_pool, solved.variable_values);
        if (!branch.has_value()) {
            interrupted_node_bound = node.lower_bound;
            result.timed_out = true;
            break;
        }

        SearchNode together = node;
        SearchNode separate = node;
        together.id = next_node_id++;
        separate.id = next_node_id++;
        together.depth = node.depth + 1;
        separate.depth = node.depth + 1;
        together.pattern_pool = solved.pattern_pool;
        separate.pattern_pool = solved.pattern_pool;
        together.sr_cuts = solved.sr_cuts;
        separate.sr_cuts = solved.sr_cuts;
        together.together.emplace_back(branch->first, branch->second);
        separate.separate.emplace_back(branch->first, branch->second);
        ++statistics.rf_branches;
        stack.push_back(std::move(separate));
        stack.push_back(std::move(together));
    }

    if (stack.empty() && !interrupted_node_bound.has_value()) {
        result.optimal = true;
        result.certified_lower_bound = result.incumbent_value;
    } else {
        result.timed_out = true;
        result.certified_lower_bound = open_tree_lower_bound(
            stack, result.incumbent_value, interrupted_node_bound);
    }
    return result;
}

}
