#include "bin_indexed_root_bound.hpp"

#include "precpack/exact_arithmetic.hpp"

#include <gurobi_c++.h>

#include <algorithm>
#include <cassert>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace precpack {
namespace {

using Clock = std::chrono::steady_clock;

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

struct RootInput {
    int lower_bound = 0;
    std::shared_ptr<const std::vector<Pattern>> pattern_pool;
};

struct CompiledRoot {
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

[[nodiscard]] CompiledRoot compile_root(
    const Instance& instance,
    int bin_count,
    const std::vector<std::pair<int, int>>& base_conflicts) {
    const int n = instance.size();
    CompiledRoot result;
    result.bin_count = bin_count;
    result.component_of_item.resize(static_cast<std::size_t>(n));
    result.members.resize(static_cast<std::size_t>(n));
    result.weights.resize(static_cast<std::size_t>(n));
    for (int item = 0; item < n; ++item) {
        result.component_of_item[static_cast<std::size_t>(item)] = item;
        result.members[static_cast<std::size_t>(item)].push_back(item);
        result.weights[static_cast<std::size_t>(item)] =
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
        for (int bin = 0; bin < bin_count; ++bin) {
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
                                      const CompiledRoot& compiled) {
    if (pattern.bin < 0 || pattern.bin >= compiled.bin_count) {
        return false;
    }
    std::vector<int> counts(
        static_cast<std::size_t>(compiled.component_count()), 0);
    std::vector<int> selected;
    selected.reserve(pattern.items.size());
    for (const int item : pattern.items) {
        const int component =
            compiled.component_of_item[static_cast<std::size_t>(item)];
        if (counts[static_cast<std::size_t>(component)]++ == 0) {
            selected.push_back(component);
        }
    }
    for (const int component : selected) {
        if (counts[static_cast<std::size_t>(component)] !=
                static_cast<int>(compiled.members[
                    static_cast<std::size_t>(component)].size()) ||
            !compiled.is_eligible(component, pattern.bin)) {
            return false;
        }
    }
    for (std::size_t i = 0; i < selected.size(); ++i) {
        for (std::size_t j = i + 1; j < selected.size(); ++j) {
            if (compiled.conflicts(selected[i], selected[j])) {
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
    };

    RestrictedMaster(GRBEnv& environment,
                     const Instance& instance,
                     const CompiledRoot& compiled,
                     const RootInput& root,
                     const Config& config,
                     Statistics& statistics)
        : instance_(instance),
          compiled_(compiled),
          bin_count_(compiled.bin_count),
          model_(environment),
          active_arc_position_(instance.arcs.size(), -1) {
        model_.set(GRB_IntParam_Threads, 1);
        model_.set(GRB_IntParam_Seed, config.seed);
        model_.set(GRB_IntParam_OutputFlag, 0);
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
                model_.addConstr(GRBLinExpr(0.0) >= 1.0));
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
        const std::size_t expected = root.pattern_pool->size() + 64U;
        patterns_.reserve(expected);
        variables_.reserve(expected);
        variable_values_.reserve(expected);
        pattern_keys_.reserve(expected * 2U + 1U);
        for (const Pattern& pattern : *root.pattern_pool) {
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
        variables_.push_back(
            model_.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, column));
        patterns_.push_back(std::move(pattern));
    }

    [[nodiscard]] Duals duals() {
        Duals result;
        result.item.reserve(item_rows_.size());
        for (GRBConstr& row : item_rows_) {
            const double dual = row.get(GRB_DoubleAttr_Pi);
            result.item.push_back(std::max(0.0, dual));
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
    [[nodiscard]] int bin_count() const noexcept { return bin_count_; }
    [[nodiscard]] std::uint64_t column_count() const noexcept {
        return patterns_.size();
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

    const Instance& instance_;
    const CompiledRoot& compiled_;
    int bin_count_ = 0;
    GRBModel model_;
    std::vector<GRBVar> y_;
    std::vector<GRBConstr> item_rows_;
    std::vector<GRBConstr> usage_rows_;
    std::vector<GRBConstr> prefix_rows_;
    std::vector<int> active_arc_position_;
    std::vector<int> active_arc_indices_;
    std::vector<GRBConstr> active_arc_rows_;
    std::vector<GRBVar> artificial_variables_;
    std::vector<Pattern> patterns_;
    std::vector<GRBVar> variables_;
    std::vector<double> variable_values_;
    std::vector<double> y_values_;
    std::unordered_set<std::string> pattern_keys_;
    bool phase_one_ = true;
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
        const CompiledRoot& compiled,
        const std::unordered_set<std::string>& existing_patterns,
        int bin,
        std::vector<double> item_profit,
        double base_cost,
        double tolerance,
        int maximum_columns,
        Clock::time_point end_time)
        : instance_(instance),
          compiled_(compiled),
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
          diversity_count_(static_cast<std::size_t>(instance.size()), 0) {
        std::vector<Group> easy;
        std::vector<Group> special;
        easy.reserve(static_cast<std::size_t>(compiled.component_count()));
        special.reserve(static_cast<std::size_t>(compiled.component_count()));
        for (int component = 0; component < compiled.component_count(); ++component) {
            if (!compiled.is_eligible(component, bin)) {
                continue;
            }
            double profit = 0.0;
            for (const int item :
                 compiled.members[static_cast<std::size_t>(component)]) {
                profit += item_profit[static_cast<std::size_t>(item)];
            }
            if (profit <= 0.0) {
                continue;
            }
            Group group{component,
                        compiled.weights[static_cast<std::size_t>(component)],
                        profit};
            if (!compiled.conflict_neighbors[static_cast<std::size_t>(component)]
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
    const CompiledRoot& compiled_;
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
                         const CompiledRoot& compiled,
                         int bin,
                         std::vector<std::int64_t> item_profit,
                         Clock::time_point end_time)
        : instance_(instance),
          compiled_(compiled),
          end_time_(end_time),
          blocked_(static_cast<std::size_t>(compiled.component_count()), 0) {
        const int component_count = compiled.component_count();
        std::vector<std::int64_t> component_profit(
            static_cast<std::size_t>(component_count), 0);
        for (int component = 0; component < component_count; ++component) {
            for (const int item :
                 compiled.members[static_cast<std::size_t>(component)]) {
                component_profit[static_cast<std::size_t>(component)] =
                    checked_add(
                        component_profit[static_cast<std::size_t>(component)],
                        item_profit[static_cast<std::size_t>(item)]);
            }
        }

        std::vector<Group> easy;
        std::vector<Group> special;
        easy.reserve(static_cast<std::size_t>(component_count));
        special.reserve(static_cast<std::size_t>(component_count));
        for (int component = component_count - 1; component >= 0; --component) {
            if (!compiled.is_eligible(component, bin) ||
                component_profit[static_cast<std::size_t>(component)] <= 0 ||
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
    const CompiledRoot& compiled_;
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
    std::int64_t best_profit_ = 0;
    std::uint64_t nodes_ = 0;
    bool timed_out_ = false;
};

struct ScaledDuals {
    std::int64_t scale = 1;
    std::vector<std::int64_t> item;
    std::vector<std::int64_t> arc;
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
    const CompiledRoot& compiled,
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
        IntegerPricingSearch pricing(instance, compiled, bin,
                                     std::move(item_profit), end_time);
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
    const CompiledRoot& compiled,
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
            instance, compiled, master.pattern_keys(), bin,
            std::move(item_profit), base_cost, tolerance,
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

enum class RootSolveState {
    kConverged,
    kInfeasible,
    kTimedOut,
};

struct RootSolveResult {
    RootSolveState state = RootSolveState::kTimedOut;
    int certified_lower_bound = 0;
    double lp_value = 0.0;
    std::uint64_t column_count = 0;
};

[[nodiscard]] RootSolveResult solve_root(
    GRBEnv& environment,
    const Instance& instance,
    const CompiledRoot& compiled,
    const RootInput& root,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    RootSolveResult result;
    result.certified_lower_bound = root.lower_bound;
    const auto start = Clock::now();
    RestrictedMaster master(environment, instance, compiled, root, config,
                            statistics);

    for (int iteration = 0; iteration < config.max_cg_iterations; ++iteration) {
        ++statistics.cg_iterations;
        if (deadline.expired()) {
            result.state = RootSolveState::kTimedOut;
            break;
        }
        if (!master.solve(deadline.end_time(), statistics)) {
            result.state = RootSolveState::kTimedOut;
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
                result.state = RootSolveState::kTimedOut;
                break;
            }
            if (pricing.added_columns) {
                continue;
            }
            result.state = RootSolveState::kInfeasible;
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
            result.state = RootSolveState::kTimedOut;
            break;
        }
        if (pricing.added_columns) {
            continue;
        }

        result.lp_value = master.objective_value();
        const SafeBoundResult safe = compute_safe_lagrangian_bound(
            instance, compiled, master, duals, deadline.end_time(), statistics);
        if (!safe.proven) {
            result.state = RootSolveState::kTimedOut;
            break;
        }
        result.certified_lower_bound =
            std::max(root.lower_bound, safe.integer_bound);
        result.state = RootSolveState::kConverged;
        break;
    }

    statistics.cg_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    ++statistics.cg_count;
    result.column_count = master.column_count();
    return result;
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

}

BinIndexedRootBoundResult run_bin_indexed_root_bound(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    BinIndexedRootBoundResult result;
    result.attempted = true;
    result.certified_lower_bound = lower_bound;
    if (lower_bound >= incumbent.bin_count) {
        result.completed = true;
        result.certified_lower_bound = incumbent.bin_count;
        return result;
    }
    if (deadline.expired()) {
        result.timed_out = true;
        return result;
    }

    const std::vector<std::pair<int, int>> base_conflicts =
        build_base_conflict_edges(instance);
    RootInput root;
    root.lower_bound = lower_bound;
    root.pattern_pool = initial_patterns(instance, incumbent);

    ++statistics.explored_nodes;
    const CompiledRoot compiled = compile_root(
        instance, incumbent.bin_count, base_conflicts);
    if (compiled.infeasible) {
        ++statistics.infeasible_nodes;
        result.completed = true;
        result.certified_lower_bound = incumbent.bin_count;
        return result;
    }

    const RootSolveResult solved = solve_root(
        environment, instance, compiled, root, config, deadline,
        statistics);
    result.column_count = solved.column_count;
    result.certified_lower_bound = std::min(
        incumbent.bin_count,
        std::max(lower_bound, solved.certified_lower_bound));
    if (solved.state == RootSolveState::kTimedOut) {
        result.timed_out = true;
        return result;
    }
    if (solved.state == RootSolveState::kInfeasible) {
        ++statistics.infeasible_nodes;
        result.completed = true;
        result.certified_lower_bound = incumbent.bin_count;
        return result;
    }

    result.completed = true;
    result.lp_value = solved.lp_value;
    return result;
}

}
