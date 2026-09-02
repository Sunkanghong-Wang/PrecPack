#include "root_column_generation.hpp"

#include "precpack/exact_arithmetic.hpp"

#include "conflict_bin_packing.hpp"

#include <gurobi_c++.h>

#include "gurobi_compat.hpp"

#include <algorithm>
#include <cassert>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace precpack {
namespace {

using Clock = std::chrono::steady_clock;

inline constexpr double kMaximumObjectiveMultiplier = 400.0;
inline constexpr double kScaledObjectiveTarget = 900000.0;
inline constexpr double kGurobiOptimalityTolerance = 1e-9;
inline constexpr double kPricingResolutionSafety = 1.5;

[[nodiscard]] double rmp_objective_multiplier(int upper_bound) {
    const double range_limited = std::floor(
        kScaledObjectiveTarget /
        static_cast<double>(std::max(1, upper_bound)));
    return std::clamp(range_limited, 1.0, kMaximumObjectiveMultiplier);
}

[[nodiscard]] std::int64_t pricing_scale_limit(
    double objective_multiplier) {
    const long double limit =
        static_cast<long double>(objective_multiplier) /
        (static_cast<long double>(kPricingResolutionSafety) *
         static_cast<long double>(kGurobiOptimalityTolerance));
    if (!std::isfinite(limit) || limit < 1.0L ||
        limit >= static_cast<long double>(
                     std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("invalid position-free pricing scale limit");
    }
    return static_cast<std::int64_t>(std::floor(limit));
}

struct Pattern {
    std::vector<int> items;
    std::vector<std::uint64_t> bits;
};

[[nodiscard]] std::string pattern_key(
    const std::vector<std::uint64_t>& bits) {
    std::string key(bits.size() * sizeof(std::uint64_t), '\0');
    if (!bits.empty()) {
        std::memcpy(key.data(), bits.data(), key.size());
    }
    return key;
}

struct ConflictGraph {
    std::vector<std::vector<int>> neighbors;
    std::vector<int> java_item_order;
};

[[nodiscard]] ConflictGraph build_conflict_graph(const Instance& instance) {
    ConflictGraph graph;
    graph.neighbors.assign(static_cast<std::size_t>(instance.size()), {});
    std::vector<unsigned char> seen(static_cast<std::size_t>(instance.size()),
                                    0U);
    for (int lhs = 0; lhs < instance.size(); ++lhs) {
        for (int rhs = lhs + 1; rhs < instance.size(); ++rhs) {
            if (internal::items_have_same_bin_conflict(instance, lhs, rhs)) {
                graph.neighbors[static_cast<std::size_t>(lhs)].push_back(rhs);
                graph.neighbors[static_cast<std::size_t>(rhs)].push_back(lhs);
                if (seen[static_cast<std::size_t>(lhs)] == 0U) {
                    seen[static_cast<std::size_t>(lhs)] = 1U;
                    graph.java_item_order.push_back(lhs);
                }
                if (seen[static_cast<std::size_t>(rhs)] == 0U) {
                    seen[static_cast<std::size_t>(rhs)] = 1U;
                    graph.java_item_order.push_back(rhs);
                }
            }
        }
    }
    return graph;
}

[[nodiscard]] std::vector<Pattern> make_initial_patterns(
    const Instance& instance, const Assignment& assignment) {
    const std::size_t blocks =
        (static_cast<std::size_t>(instance.size()) + 63U) / 64U;
    std::vector<Pattern> patterns(static_cast<std::size_t>(assignment.bin_count));
    for (Pattern& pattern : patterns) {
        pattern.bits.assign(blocks, 0U);
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
                       [](const Pattern& pattern) {
                           return pattern.items.empty();
                       }),
        patterns.end());
    return patterns;
}

class PositionFreeMaster {
public:
    struct Duals {
        std::vector<double> item;
        double pattern_count = 0.0;
    };

    PositionFreeMaster(GRBEnv& environment,
                       const Instance& instance,
                       int lower_bound,
                       int upper_bound,
        const Config& config)
        : objective_multiplier_(rmp_objective_multiplier(upper_bound)),
          model_(environment),
          horizon_(model_.addVar(
              static_cast<double>(lower_bound),
              static_cast<double>(upper_bound), objective_multiplier_,
              GRB_CONTINUOUS)) {
        model_.set(GRB_IntParam_Threads, 1);
        model_.set(GRB_IntParam_Seed, config.seed);
        model_.set(GRB_IntParam_OutputFlag, 0);
        model_.set(GRB_IntParam_Presolve, 0);
        model_.set(GRB_IntParam_Method, 0);
        model_.set(GRB_IntParam_NumericFocus, 2);
        model_.set(GRB_IntParam_ScaleFlag, 0);
        model_.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
        model_.set(GRB_DoubleParam_OptimalityTol, 1e-9);

        model_.update();

        item_rows_.reserve(static_cast<std::size_t>(instance.size()));
        for (int item = 0; item < instance.size(); ++item) {
            item_rows_.push_back(
                model_.addConstr(GRBLinExpr(0.0) >= 1.0));
        }
        GRBLinExpr expression = horizon_;
        pattern_count_row_ = model_.addConstr(expression >= 0.0);
        model_.update();

        const std::size_t expected =
            static_cast<std::size_t>(std::max(1024, instance.size() * 4));
        patterns_.reserve(expected);
        variables_.reserve(expected);
        pattern_keys_.reserve(expected * 2U);
    }

    void add_initial_patterns(std::vector<Pattern> patterns) {
        for (Pattern& pattern : patterns) {
            static_cast<void>(add_pattern(std::move(pattern)));
        }
        model_.update();
    }

    [[nodiscard]] bool add_pattern(Pattern pattern) {
        const std::string key = pattern_key(pattern.bits);
        if (!pattern_keys_.insert(key).second) {
            return false;
        }
        GRBColumn column;
        for (const int item : pattern.items) {
            column.addTerm(1.0,
                           item_rows_[static_cast<std::size_t>(item)]);
        }
        column.addTerm(-1.0, pattern_count_row_);
        variables_.push_back(model_.addVar(
            0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, column));
        patterns_.push_back(std::move(pattern));
        return true;
    }

    void update() { model_.update(); }

    [[nodiscard]] bool solve(Clock::time_point end_time,
                             Statistics& statistics) {
        ++statistics.rmp_count;
        int status = optimize_once(end_time, statistics);
        if (status != GRB_OPTIMAL) {
            set_timeout_status(status, end_time);
            return false;
        }
        objective_value_ =
            model_.get(GRB_DoubleAttr_ObjVal) / objective_multiplier_;
        return true;
    }

    [[nodiscard]] Duals duals() {
        Duals result;
        result.item.reserve(item_rows_.size());
        for (GRBConstr& row : item_rows_) {
            const double raw = row.get(GRB_DoubleAttr_Pi);
            if (!std::isfinite(raw)) {
                throw std::runtime_error(
                    "Gurobi returned a non-finite item dual");
            }
            result.item.push_back(
                std::max(0.0, raw / objective_multiplier_));
        }
        const double raw = pattern_count_row_.get(GRB_DoubleAttr_Pi);
        if (!std::isfinite(raw)) {
            throw std::runtime_error(
                "Gurobi returned a non-finite pattern-count dual");
        }
        const double normalized = raw / objective_multiplier_;
        result.pattern_count = std::max(0.0, normalized);
        return result;
    }

    [[nodiscard]] double objective_value() const noexcept {
        return objective_value_;
    }
    [[nodiscard]] bool timed_out() const noexcept { return timed_out_; }
    [[nodiscard]] std::size_t pattern_count() const noexcept {
        return patterns_.size();
    }
    [[nodiscard]] double objective_multiplier() const noexcept {
        return objective_multiplier_;
    }
    [[nodiscard]] bool link_objective_multiplier_to_scale(
        std::int64_t certificate_scale) {
        if (certificate_scale <= 0) {
            throw std::invalid_argument(
                "nonpositive position-free certificate scale");
        }
        const long double required = std::ceil(
            static_cast<long double>(kPricingResolutionSafety) *
            static_cast<long double>(kGurobiOptimalityTolerance) *
            static_cast<long double>(certificate_scale));
        if (required >= static_cast<long double>(objective_multiplier_)) {
            return false;
        }
        const double linked = std::max(1.0, static_cast<double>(required));
        objective_multiplier_ = linked;
        horizon_.set(GRB_DoubleAttr_Obj, linked);
        model_.update();
        return true;
    }
    [[nodiscard]] const std::unordered_set<std::string>& pattern_keys() const
        noexcept {
        return pattern_keys_;
    }

private:
    [[nodiscard]] int optimize_once(Clock::time_point end_time,
                                    Statistics& statistics) {
        const double remaining =
            std::chrono::duration<double>(end_time - Clock::now()).count();
        if (remaining <= 0.0) {
            return GRB_TIME_LIMIT;
        }
        model_.set(GRB_DoubleParam_TimeLimit, remaining);
        const auto start = Clock::now();
        model_.optimize();
        statistics.rmp_seconds +=
            std::chrono::duration<double>(Clock::now() - start).count();
        return model_.get(GRB_IntAttr_Status);
    }

    void set_timeout_status(int status, Clock::time_point end_time) noexcept {
        timed_out_ = status == GRB_TIME_LIMIT || status == GRB_INTERRUPTED ||
                     gurobi_compat::is_work_limit_status(status) ||
                     Clock::now() >= end_time;
    }

    double objective_multiplier_ = kMaximumObjectiveMultiplier;
    GRBModel model_;
    GRBVar horizon_;
    std::vector<GRBConstr> item_rows_;
    GRBConstr pattern_count_row_;
    std::vector<Pattern> patterns_;
    std::vector<GRBVar> variables_;
    std::unordered_set<std::string> pattern_keys_;
    bool timed_out_ = false;
    double objective_value_ = 0.0;
};

[[nodiscard]] std::int64_t checked_add(std::int64_t lhs,
                                       std::int64_t rhs) {
    return exact_arithmetic::checked_add(
        lhs, rhs, "position-free pricing overflow");
}

struct LongPricedPattern {
    Pattern pattern;
    std::int64_t profit = 0;
};

struct LongPricingResult {
    bool proven = true;
    std::int64_t maximum_profit = 0;
    std::uint64_t nodes = 0;
    std::vector<LongPricedPattern> patterns;
};

template <typename Profit>
struct PositionFreeCandidate {
    int item = -1;
    int weight = 0;
    Profit profit = 0;
};

inline constexpr int kMaximumDpChoiceSpecialCandidates = 128;

template <typename Candidate>
class FractionalKnapsackBound {
public:
    void initialize(const std::vector<Candidate>& candidates,
                    int first_count,
                    int last_count) {
        const std::size_t candidate_count = candidates.size();
        if (first_count < 0 || last_count < first_count ||
            static_cast<std::size_t>(last_count) > candidate_count) {
            throw std::invalid_argument(
                "invalid fractional-knapsack prefix range");
        }
        first_count_ = first_count;
        last_count_ = last_count;
        for (const Candidate& candidate : candidates) {
            if (candidate.weight < 0 || candidate.profit <= 0) {
                throw std::invalid_argument(
                    "invalid fractional-knapsack candidate");
            }
        }
        std::vector<int> density_order(candidate_count);
        for (std::size_t index = 0; index < candidate_count; ++index) {
            density_order[index] = static_cast<int>(index);
        }
        std::sort(density_order.begin(), density_order.end(),
                  [&candidates](int lhs_index, int rhs_index) {
                      const Candidate& lhs =
                          candidates[static_cast<std::size_t>(lhs_index)];
                      const Candidate& rhs =
                          candidates[static_cast<std::size_t>(rhs_index)];
                      if (lhs.weight == 0 || rhs.weight == 0) {
                          if (lhs.weight != rhs.weight) {
                              return lhs.weight == 0;
                          }
                      } else {
                          const int comparison =
                              exact_arithmetic::compare_nonnegative_fractions(
                                  static_cast<std::uint64_t>(lhs.profit),
                                  static_cast<std::uint64_t>(lhs.weight),
                                  static_cast<std::uint64_t>(rhs.profit),
                                  static_cast<std::uint64_t>(rhs.weight));
                          if (comparison != 0) {
                              return comparison > 0;
                          }
                      }
                      return lhs_index < rhs_index;
                  });

        density_position_.assign(candidate_count, -1);
        density_item_.assign(candidate_count, -1);
        density_weight_.assign(candidate_count, 0);
        density_profit_.assign(candidate_count, 0);
        for (std::size_t density = 0; density < candidate_count; ++density) {
            const int candidate_index = density_order[density];
            const Candidate& candidate =
                candidates[static_cast<std::size_t>(candidate_index)];
            density_position_[static_cast<std::size_t>(candidate_index)] =
                static_cast<int>(density);
            density_item_[density] = candidate.item;
            density_weight_[density] = candidate.weight;
            density_profit_[density] = candidate.profit;
        }

        bit_blocks_ = (candidate_count + 63U) / 64U;
        const std::size_t row_count = static_cast<std::size_t>(
            last_count_ - first_count_ + 1);
        if (bit_blocks_ != 0U &&
            row_count >
                std::numeric_limits<std::size_t>::max() / bit_blocks_) {
            throw std::length_error(
                "fractional-knapsack prefix bitset is too large");
        }
        prefix_density_bits_.assign(
            row_count * bit_blocks_, 0U);
        for (int position = 0; position < first_count_; ++position) {
            const std::size_t density = static_cast<std::size_t>(
                density_position_[static_cast<std::size_t>(position)]);
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
            const std::size_t candidate_index =
                static_cast<std::size_t>(count - 1);
            const std::size_t density = static_cast<std::size_t>(
                density_position_[candidate_index]);
            prefix_density_bits_[current_offset + density / 64U] |=
                std::uint64_t{1} << (density & 63U);
        }
    }

    [[nodiscard]] std::int64_t upper_bound(
        int count,
        int capacity,
        const std::vector<int>& blocked) const {
        assert(count >= first_count_);
        assert(count <= last_count_);
        assert(capacity >= 0);
        std::int64_t remaining = capacity;
        std::int64_t profit = 0;
        const std::size_t row_offset =
            static_cast<std::size_t>(count - first_count_) * bit_blocks_;
        for (std::size_t block = 0; block < bit_blocks_; ++block) {
            std::uint64_t active =
                prefix_density_bits_[row_offset + block];
            while (active != 0U) {
                const unsigned bit = std::countr_zero(active);
                const std::size_t density = block * 64U + bit;
                if (blocked[static_cast<std::size_t>(
                        density_item_[density])] != 0) {
                    active &= active - 1U;
                    continue;
                }
                const std::int64_t weight = density_weight_[density];
                const std::int64_t item_profit = density_profit_[density];
                if (weight <= remaining) {
                    remaining -= weight;
                    profit = checked_add(profit, item_profit);
                } else {
                    assert(weight > 0);
                    const std::int64_t fractional =
                        exact_arithmetic::ceil_nonnegative_product_ratio(
                            item_profit, remaining, weight,
                            "fractional-knapsack bound overflow");
                    return checked_add(
                        profit, fractional);
                }
                active &= active - 1U;
            }
        }
        return profit;
    }

private:
    int first_count_ = 0;
    int last_count_ = 0;
    std::size_t bit_blocks_ = 0U;
    std::vector<int> density_position_;
    std::vector<int> density_item_;
    std::vector<std::int64_t> density_weight_;
    std::vector<std::int64_t> density_profit_;
    std::vector<std::uint64_t> prefix_density_bits_;
};

template <typename Profit>
[[nodiscard]] std::pair<std::vector<PositionFreeCandidate<Profit>>, int>
java_ordered_candidates(const Instance& instance,
                        const ConflictGraph& conflicts,
                        const std::vector<Profit>& item_profit) {
    using Candidate = PositionFreeCandidate<Profit>;
    std::vector<Candidate> easy;
    std::vector<Candidate> conflict;
    easy.reserve(static_cast<std::size_t>(instance.size()));
    conflict.reserve(static_cast<std::size_t>(instance.size()));

    for (int item = instance.size() - 1; item >= 0; --item) {
        const Profit profit = item_profit[static_cast<std::size_t>(item)];
        if (profit > Profit{0} &&
            conflicts.neighbors[static_cast<std::size_t>(item)].empty()) {
            easy.push_back(Candidate{
                item, instance.items[static_cast<std::size_t>(item)].weight,
                profit});
        }
    }
    for (const int item : conflicts.java_item_order) {
        const Profit profit = item_profit[static_cast<std::size_t>(item)];
        if (profit > Profit{0} &&
            !conflicts.neighbors[static_cast<std::size_t>(item)].empty()) {
            conflict.push_back(Candidate{
                item, instance.items[static_cast<std::size_t>(item)].weight,
                profit});
        }
    }
    std::stable_sort(conflict.begin(), conflict.end(),
                     [](const Candidate& lhs, const Candidate& rhs) {
                         return lhs.weight < rhs.weight;
                     });

    const int easy_count = static_cast<int>(easy.size());
    easy.insert(easy.end(), conflict.begin(), conflict.end());
    return {std::move(easy), easy_count};
}

class LongPricingSearch {
public:
    LongPricingSearch(
        const Instance& instance,
        const ConflictGraph& conflicts,
        const std::unordered_set<std::string>& existing_patterns,
        const std::vector<std::int64_t>& item_profit,
        std::int64_t threshold,
        Clock::time_point end_time)
        : instance_(instance),
          conflicts_(conflicts),
          existing_patterns_(existing_patterns),
          threshold_(threshold),
          end_time_(end_time),
          blocked_(static_cast<std::size_t>(instance.size()), 0),
          selected_bits_((static_cast<std::size_t>(instance.size()) + 63U) /
                             64U,
                         0U),
          heuristic_item_use_count_(static_cast<std::size_t>(instance.size()),
                                    0) {
        auto [ordered, easy_count] =
            java_ordered_candidates(instance, conflicts, item_profit);
        candidates_ = std::move(ordered);
        easy_count_ = easy_count;
        use_dp_choice_gate_ =
            static_cast<int>(candidates_.size()) - easy_count_ <=
            kMaximumDpChoiceSpecialCandidates;
        if (use_dp_choice_gate_) {
            fractional_bound_.initialize(
                candidates_, easy_count_,
                static_cast<int>(candidates_.size()));
        }

        const std::size_t width =
            static_cast<std::size_t>(instance.capacity + 1);
        dp_.assign((candidates_.size() + 1U) * width, 0);
        if (use_dp_choice_gate_) {
            choice_blocks_ = (candidates_.size() + 63U) / 64U;
            dp_choice_bits_.assign(width * choice_blocks_, 0U);
        }
        for (std::size_t count = 1; count <= candidates_.size(); ++count) {
            const Candidate& candidate = candidates_[count - 1U];
            const std::int64_t* previous =
                dp_.data() + (count - 1U) * width;
            std::int64_t* current = dp_.data() + count * width;
            for (int capacity = 0; capacity <= instance.capacity; ++capacity) {
                current[static_cast<std::size_t>(capacity)] =
                    previous[static_cast<std::size_t>(capacity)];
                if (candidate.weight <= capacity) {
                    const std::int64_t included = checked_add(
                        previous[static_cast<std::size_t>(
                            capacity - candidate.weight)],
                        candidate.profit);
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
        selected_items_.reserve(candidates_.size());
        heuristic_patterns_.reserve(candidates_.size());
        heuristic_keys_.reserve(candidates_.size());
    }

    [[nodiscard]] LongPricingResult solve() {
        LongPricingResult result;
        if (Clock::now() >= end_time_) {
            result.proven = false;
            return result;
        }

        heuristic_mode_ = true;
        heuristic_limit_ = std::max<std::uint64_t>(
            1U, static_cast<std::uint64_t>(std::max(1, instance_.size())) *
                    static_cast<std::uint64_t>(
                        std::max(1, instance_.capacity)) /
                    10U);
        if (use_dp_choice_gate_) {
            search<true>(static_cast<int>(candidates_.size()),
                         instance_.capacity, 0);
        } else {
            search<false>(static_cast<int>(candidates_.size()),
                          instance_.capacity, 0);
        }

        heuristic_mode_ = false;
        heuristic_stopped_ = false;
        best_new_profit_ = threshold_;
        for (const LongPricedPattern& pattern : heuristic_patterns_) {
            best_new_profit_ = std::max(best_new_profit_, pattern.profit);
        }
        best_new_pattern_.reset();
        const std::int64_t relaxed_upper_bound =
            dp_bound(static_cast<int>(candidates_.size()), instance_.capacity);
        if (relaxed_upper_bound > best_new_profit_) {
            if (use_dp_choice_gate_) {
                search<true>(static_cast<int>(candidates_.size()),
                             instance_.capacity, 0);
            } else {
                search<false>(static_cast<int>(candidates_.size()),
                              instance_.capacity, 0);
            }
        }

        result.proven = !timed_out_;
        result.maximum_profit = maximum_profit_;
        result.nodes = nodes_;
        if (timed_out_) {
            return result;
        }

        if (best_new_pattern_.has_value()) {
            result.patterns.push_back(*best_new_pattern_);
        }
        std::stable_sort(heuristic_patterns_.begin(), heuristic_patterns_.end(),
                         [](const LongPricedPattern& lhs,
                            const LongPricedPattern& rhs) {
                             return lhs.profit > rhs.profit;
                         });
        std::vector<unsigned char> final_item_used(
            static_cast<std::size_t>(instance_.size()), 0U);
        for (const LongPricedPattern& pattern : heuristic_patterns_) {
            bool can_add = true;
            for (const int item : pattern.pattern.items) {
                if (final_item_used[static_cast<std::size_t>(item)] != 0U) {
                    can_add = false;
                    break;
                }
            }
            if (!can_add) {
                continue;
            }
            for (const int item : pattern.pattern.items) {
                final_item_used[static_cast<std::size_t>(item)] = 1U;
            }
            result.patterns.push_back(pattern);
        }
        return result;
    }

private:
    using Candidate = PositionFreeCandidate<std::int64_t>;

    [[nodiscard]] std::int64_t dp_bound(int count,
                                        int capacity) const noexcept {
        const std::size_t width =
            static_cast<std::size_t>(instance_.capacity + 1);
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
            const Candidate& candidate =
                candidates_[static_cast<std::size_t>(position)];
            if (blocked_[static_cast<std::size_t>(candidate.item)] != 0) {
                return true;
            }
            capacity -= candidate.weight;
            count = position;
        }
        return false;
    }

    [[nodiscard]] Pattern make_pattern() const {
        Pattern pattern;
        pattern.items = selected_items_;
        std::sort(pattern.items.begin(), pattern.items.end());
        pattern.bits = selected_bits_;
        return pattern;
    }

    [[nodiscard]] bool selected_pattern_exists(bool include_heuristic) const {
        const std::string key = pattern_key(selected_bits_);
        return existing_patterns_.contains(key) ||
               (include_heuristic && heuristic_keys_.contains(key));
    }

    [[nodiscard]] bool consider(std::int64_t profit) {
        maximum_profit_ = std::max(maximum_profit_, profit);
        if (selected_items_.empty() || profit <= threshold_ ||
            selected_pattern_exists(!heuristic_mode_)) {
            return false;
        }
        if (heuristic_mode_) {
            for (const int item : selected_items_) {
                if (heuristic_item_use_count_[static_cast<std::size_t>(item)] >=
                    2) {
                    return false;
                }
            }
            Pattern pattern = make_pattern();
            const std::string key = pattern_key(pattern.bits);
            if (heuristic_keys_.insert(key).second) {
                for (const int item : pattern.items) {
                    ++heuristic_item_use_count_[static_cast<std::size_t>(item)];
                }
                heuristic_patterns_.push_back(
                    LongPricedPattern{std::move(pattern), profit});
                return true;
            }
        } else if (profit > best_new_profit_) {
            best_new_profit_ = profit;
            best_new_pattern_ =
                LongPricedPattern{make_pattern(), profit};
            return true;
        }
        return false;
    }

    void include(const Candidate& candidate, std::int64_t& profit) {
        selected_items_.push_back(candidate.item);
        selected_bits_[static_cast<std::size_t>(candidate.item) / 64U] |=
            std::uint64_t{1} <<
            (static_cast<unsigned>(candidate.item) & 63U);
        profit = checked_add(profit, candidate.profit);
        for (const int neighbor :
             conflicts_.neighbors[static_cast<std::size_t>(candidate.item)]) {
            ++blocked_[static_cast<std::size_t>(neighbor)];
        }
    }

    void remove(const Candidate& candidate,
                std::int64_t& profit) {
        for (const int neighbor :
             conflicts_.neighbors[static_cast<std::size_t>(candidate.item)]) {
            int& block_count = blocked_[static_cast<std::size_t>(neighbor)];
            assert(block_count > 0);
            --block_count;
        }
        profit = checked_add(profit, -candidate.profit);
        selected_bits_[static_cast<std::size_t>(candidate.item) / 64U] &=
            ~(std::uint64_t{1} <<
              (static_cast<unsigned>(candidate.item) & 63U));
        selected_items_.pop_back();
    }

    [[nodiscard]] bool selected_item_is_saturated() const noexcept {
        for (const int item : selected_items_) {
            if (heuristic_item_use_count_[static_cast<std::size_t>(item)] >= 2) {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] bool reconstruct_easy(std::int64_t target_profit,
                                        int count,
                                        int capacity,
                                        std::int64_t profit) {
        if (timed_out_) {
            return false;
        }
        ++nodes_;
        if ((nodes_ & 1023U) == 0U && Clock::now() >= end_time_) {
            timed_out_ = true;
            return false;
        }
        if (count == 0 || capacity == 0) {
            const bool added = consider(profit);
            return added && profit == target_profit;
        }

        const Candidate& candidate =
            candidates_[static_cast<std::size_t>(count - 1)];
        if (candidate.weight <= capacity) {
            const std::int64_t included_profit =
                checked_add(profit, candidate.profit);
            const std::int64_t included_bound = checked_add(
                included_profit,
                dp_bound(count - 1, capacity - candidate.weight));
            if (included_bound > best_new_profit_) {
                std::int64_t mutable_profit = profit;
                include(candidate, mutable_profit);
                const bool found = reconstruct_easy(
                    target_profit, count - 1, capacity - candidate.weight,
                    mutable_profit);
                remove(candidate, mutable_profit);
                if (found) {
                    return true;
                }
            }
        }
        const std::int64_t excluded_bound =
            checked_add(profit, dp_bound(count - 1, capacity));
        return excluded_bound > best_new_profit_ &&
               reconstruct_easy(target_profit, count - 1, capacity, profit);
    }

    template <bool UseFractionalBound>
    void search(int count, int capacity, std::int64_t profit) {
        if (timed_out_ || heuristic_stopped_) {
            return;
        }
        ++nodes_;
        if ((nodes_ & 1023U) == 0U && Clock::now() >= end_time_) {
            timed_out_ = true;
            return;
        }
        if (count == 0 || capacity == 0) {
            static_cast<void>(consider(profit));
            return;
        }
        if (heuristic_mode_ && ++heuristic_expansions_ > heuristic_limit_) {
            heuristic_stopped_ = true;
            return;
        }
        const std::int64_t target =
            heuristic_mode_ ? threshold_
                            : std::max(threshold_, best_new_profit_);
        const std::int64_t upper_bound =
            checked_add(profit, dp_bound(count, capacity));
        if (upper_bound <= target) {
            return;
        }
        if constexpr (UseFractionalBound) {
            if (count > easy_count_ && profit > 0 &&
                dp_solution_uses_blocked(count, capacity)) {
                const std::int64_t fractional_upper_bound = checked_add(
                    profit,
                    fractional_bound_.upper_bound(count, capacity, blocked_));
                if (fractional_upper_bound <= target) {
                    return;
                }
            }
        }
        if (!heuristic_mode_ && count == easy_count_) {
            static_cast<void>(reconstruct_easy(upper_bound, count, capacity,
                                               profit));
            return;
        }

        const Candidate& candidate =
            candidates_[static_cast<std::size_t>(count - 1)];
        if (candidate.weight <= capacity &&
            blocked_[static_cast<std::size_t>(candidate.item)] == 0) {
            std::int64_t included_profit = profit;
            include(candidate, included_profit);
            search<UseFractionalBound>(
                count - 1, capacity - candidate.weight, included_profit);
            remove(candidate, included_profit);
            if (heuristic_mode_ && selected_item_is_saturated()) {
                return;
            }
        }
        search<UseFractionalBound>(count - 1, capacity, profit);
    }

    const Instance& instance_;
    const ConflictGraph& conflicts_;
    const std::unordered_set<std::string>& existing_patterns_;
    std::int64_t threshold_ = 0;
    Clock::time_point end_time_;
    std::vector<Candidate> candidates_;
    int easy_count_ = 0;
    std::vector<std::int64_t> dp_;
    std::size_t choice_blocks_ = 0U;
    std::vector<std::uint64_t> dp_choice_bits_;
    bool use_dp_choice_gate_ = false;
    FractionalKnapsackBound<Candidate> fractional_bound_;
    std::vector<int> blocked_;
    std::vector<int> selected_items_;
    std::vector<std::uint64_t> selected_bits_;
    std::vector<LongPricedPattern> heuristic_patterns_;
    std::unordered_set<std::string> heuristic_keys_;
    std::vector<int> heuristic_item_use_count_;
    std::optional<LongPricedPattern> best_new_pattern_;
    std::int64_t maximum_profit_ = 0;
    std::int64_t best_new_profit_ = 0;
    std::uint64_t nodes_ = 0;
    std::uint64_t heuristic_limit_ = 0;
    std::uint64_t heuristic_expansions_ = 0;
    bool heuristic_mode_ = false;
    bool heuristic_stopped_ = false;
    bool timed_out_ = false;
};

struct IntegerPricingResult {
    bool proven = true;
    std::int64_t maximum_profit = 0;
    std::uint64_t nodes = 0;
};

class IntegerPricingSearch {
public:
    IntegerPricingSearch(const Instance& instance,
                         const ConflictGraph& conflicts,
                         const std::vector<std::int64_t>& item_profit,
                         Clock::time_point end_time)
        : instance_(instance),
          conflicts_(conflicts),
          end_time_(end_time),
          blocked_(static_cast<std::size_t>(instance.size()), 0) {
        auto [ordered, easy_count] =
            java_ordered_candidates(instance, conflicts, item_profit);
        candidates_ = std::move(ordered);
        easy_count_ = easy_count;
        use_dp_choice_gate_ =
            static_cast<int>(candidates_.size()) - easy_count_ <=
            kMaximumDpChoiceSpecialCandidates;
        if (use_dp_choice_gate_) {
            fractional_bound_.initialize(
                candidates_, easy_count_,
                static_cast<int>(candidates_.size()));
        }

        const std::size_t width =
            static_cast<std::size_t>(instance.capacity + 1);
        dp_.assign((candidates_.size() + 1U) * width, 0);
        if (use_dp_choice_gate_) {
            choice_blocks_ = (candidates_.size() + 63U) / 64U;
            dp_choice_bits_.assign(width * choice_blocks_, 0U);
        }
        for (std::size_t count = 1; count <= candidates_.size(); ++count) {
            const Candidate& candidate = candidates_[count - 1U];
            const std::int64_t* previous =
                dp_.data() + (count - 1U) * width;
            std::int64_t* current = dp_.data() + count * width;
            for (int capacity = 0; capacity <= instance.capacity; ++capacity) {
                current[static_cast<std::size_t>(capacity)] =
                    previous[static_cast<std::size_t>(capacity)];
                if (candidate.weight <= capacity) {
                    const std::int64_t included = checked_add(
                        previous[static_cast<std::size_t>(
                            capacity - candidate.weight)],
                        candidate.profit);
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
            search<true>(static_cast<int>(candidates_.size()),
                         instance_.capacity, 0);
        } else {
            search<false>(static_cast<int>(candidates_.size()),
                          instance_.capacity, 0);
        }
        result.proven = !timed_out_;
        result.maximum_profit = best_profit_;
        result.nodes = nodes_;
        return result;
    }

private:
    using Candidate = PositionFreeCandidate<std::int64_t>;

    [[nodiscard]] std::int64_t dp_bound(int count,
                                        int capacity) const noexcept {
        const std::size_t width =
            static_cast<std::size_t>(instance_.capacity + 1);
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
            const Candidate& candidate =
                candidates_[static_cast<std::size_t>(position)];
            if (blocked_[static_cast<std::size_t>(candidate.item)] != 0) {
                return true;
            }
            capacity -= candidate.weight;
            count = position;
        }
        return false;
    }

    void include(const Candidate& candidate) {
        for (const int neighbor :
             conflicts_.neighbors[static_cast<std::size_t>(candidate.item)]) {
            ++blocked_[static_cast<std::size_t>(neighbor)];
        }
    }

    void remove(const Candidate& candidate) {
        for (const int neighbor :
             conflicts_.neighbors[static_cast<std::size_t>(candidate.item)]) {
            int& block_count = blocked_[static_cast<std::size_t>(neighbor)];
            assert(block_count > 0);
            --block_count;
        }
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
            if (count > easy_count_ && profit > 0 &&
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
            best_profit_ = upper_bound;
            return;
        }
        const Candidate& candidate =
            candidates_[static_cast<std::size_t>(count - 1)];
        if (candidate.weight <= capacity &&
            blocked_[static_cast<std::size_t>(candidate.item)] == 0) {
            include(candidate);
            search<UseFractionalBound>(
                count - 1, capacity - candidate.weight,
                checked_add(profit, candidate.profit));
            remove(candidate);
        }
        search<UseFractionalBound>(count - 1, capacity, profit);
    }

    const Instance& instance_;
    const ConflictGraph& conflicts_;
    Clock::time_point end_time_;
    std::vector<Candidate> candidates_;
    int easy_count_ = 0;
    std::vector<std::int64_t> dp_;
    std::size_t choice_blocks_ = 0U;
    std::vector<std::uint64_t> dp_choice_bits_;
    bool use_dp_choice_gate_ = false;
    FractionalKnapsackBound<Candidate> fractional_bound_;
    std::vector<int> blocked_;
    std::int64_t best_profit_ = 0;
    std::uint64_t nodes_ = 0;
    bool timed_out_ = false;
};

struct ScaledPositionFreeDuals {
    std::int64_t scale = 1;
    std::vector<std::int64_t> item;
    std::int64_t threshold = 0;
    std::int64_t sum = 0;
};

inline constexpr std::int64_t kFixedPointBudget =
    std::int64_t{1} << 62U;

[[nodiscard]] std::int64_t truncate_nonnegative(double value,
                                                std::int64_t scale) {
    if (!std::isfinite(value) || value < 0.0 || scale <= 0) {
        throw std::invalid_argument("invalid position-free dual scaling input");
    }
    const long double scaled =
        static_cast<long double>(value) * static_cast<long double>(scale);
    constexpr long double kInt64ExclusiveUpper =
        9223372036854775808.0L;
    if (!std::isfinite(scaled) || scaled >= kInt64ExclusiveUpper) {
        throw std::overflow_error("position-free dual scaling overflow");
    }
    return static_cast<std::int64_t>(scaled);
}

[[nodiscard]] ScaledPositionFreeDuals scale_position_free_duals(
    const std::vector<double>& duals,
    double threshold,
    bool cap_at_one,
    std::int64_t maximum_scale = kFixedPointBudget) {
    if (!std::isfinite(threshold)) {
        throw std::invalid_argument("non-finite position-free threshold");
    }
    threshold = std::max(0.0, threshold);

    std::vector<double> normalized;
    normalized.reserve(duals.size());
    long double magnitude = static_cast<long double>(threshold);
    for (double value : duals) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("non-finite position-free item dual");
        }
        value = std::max(0.0, value);
        if (cap_at_one) {
            value = std::min(1.0, value);
        }
        normalized.push_back(value);
        magnitude += static_cast<long double>(value);
    }
    if (!std::isfinite(magnitude)) {
        throw std::overflow_error("position-free dual magnitude overflow");
    }

    const long double denominator = std::max(1.0L, magnitude);
    const long double candidate = std::floor(
        static_cast<long double>(kFixedPointBudget) / denominator);
    if (!std::isfinite(candidate) || candidate < 1.0L) {
        throw std::overflow_error(
            "position-free duals are too large for int64 fixed point");
    }
    std::int64_t scale = std::min(
        static_cast<std::int64_t>(candidate), maximum_scale);

    for (;;) {
        ScaledPositionFreeDuals result;
        result.scale = scale;
        result.item.reserve(normalized.size());
        for (const double value : normalized) {
            std::int64_t integer = truncate_nonnegative(value, scale);
            if (cap_at_one) {
                integer = std::min(integer, scale);
            }
            result.item.push_back(integer);
            result.sum = checked_add(result.sum, integer);
        }
        result.threshold = truncate_nonnegative(threshold, scale);
        const std::int64_t total =
            checked_add(result.sum, result.threshold);
        if (total <= kFixedPointBudget) {
            return result;
        }
        if (scale == 1) {
            throw std::overflow_error(
                "position-free fixed-point accumulation overflow");
        }
        const long double ratio =
            static_cast<long double>(kFixedPointBudget) /
            static_cast<long double>(total);
        const long double reduced =
            std::floor(static_cast<long double>(scale) * ratio);
        const std::int64_t next =
            std::max<std::int64_t>(1, static_cast<std::int64_t>(reduced));
        scale = next < scale ? next : scale - 1;
    }
}

[[nodiscard]] int ceil_ratio(std::int64_t numerator,
                             std::int64_t denominator) {
    if (denominator <= 0) {
        throw std::invalid_argument("nonpositive certificate denominator");
    }
    return exact_arithmetic::ceil_ratio_to_int(
        numerator, denominator, "position-free bound does not fit int");
}

[[nodiscard]] LongPricingResult run_long_pricing(
    const Instance& instance,
    const ConflictGraph& conflicts,
    const std::unordered_set<std::string>& existing_patterns,
    const std::vector<std::int64_t>& item_profit,
    std::int64_t threshold,
    Clock::time_point end_time,
    Statistics& statistics) {
    const auto start = Clock::now();
    LongPricingSearch pricing(instance, conflicts, existing_patterns,
                              item_profit, threshold, end_time);
    LongPricingResult result = pricing.solve();
    statistics.pricing_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    ++statistics.pricing_count;
    statistics.pricing_search_nodes += result.nodes;
    return result;
}

[[nodiscard]] IntegerPricingResult run_integer_pricing(
    const Instance& instance,
    const ConflictGraph& conflicts,
    const std::vector<std::int64_t>& item_profit,
    Clock::time_point end_time,
    Statistics& statistics) {
    const auto start = Clock::now();
    IntegerPricingSearch pricing(instance, conflicts, item_profit, end_time);
    IntegerPricingResult result = pricing.solve();
    statistics.pricing_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    ++statistics.pricing_count;
    statistics.pricing_search_nodes += result.nodes;
    return result;
}

void add_priced_patterns(PositionFreeMaster& master,
                         LongPricingResult& pricing,
                         Statistics& statistics) {
    for (LongPricedPattern& priced : pricing.patterns) {
        if (!master.add_pattern(std::move(priced.pattern))) {
            throw std::logic_error(
                "position-free pricing returned a duplicate RLMP column");
        }
        ++statistics.generated_columns;
    }
    master.update();
}

[[nodiscard]] RootStatistics make_root_statistics(
    const RootStatistics& partial,
    const PositionFreeMaster& master,
    const Statistics& statistics,
    const Statistics& initial_statistics,
    Clock::time_point start) {
    RootStatistics result = partial;
    result.column_count = master.pattern_count();
    result.generated_columns =
        statistics.generated_columns - initial_statistics.generated_columns;
    result.iterations =
        statistics.cg_iterations - initial_statistics.cg_iterations;
    result.pricing_count =
        statistics.pricing_count - initial_statistics.pricing_count;
    result.pricing_search_nodes =
        statistics.pricing_search_nodes -
        initial_statistics.pricing_search_nodes;
    result.rmp_count = statistics.rmp_count - initial_statistics.rmp_count;
    result.rmp_objective_multiplier = master.objective_multiplier();
    result.pricing_seconds =
        statistics.pricing_seconds - initial_statistics.pricing_seconds;
    result.rmp_seconds =
        statistics.rmp_seconds - initial_statistics.rmp_seconds;
    result.total_seconds =
        std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

}

RootStatistics run_position_free_root_column_generation(
    GRBEnv& environment,
    const Instance& instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& deadline,
    Statistics& statistics) {
    const auto start = Clock::now();
    const Statistics initial_statistics = statistics;
    RootStatistics partial;
    partial.attempted = true;
    partial.certified_lower_bound = lower_bound;
    if (lower_bound >= incumbent.bin_count) {
        partial.completed = true;
        partial.lp_value = static_cast<double>(lower_bound);
        PositionFreeMaster master(environment, instance, lower_bound,
                                  incumbent.bin_count, config);
        master.add_initial_patterns(make_initial_patterns(instance, incumbent));
        ++statistics.cg_count;
        statistics.cg_seconds +=
            std::chrono::duration<double>(Clock::now() - start).count();
        return make_root_statistics(partial, master, statistics,
                                    initial_statistics, start);
    }

    PositionFreeMaster master(environment, instance, lower_bound,
                              incumbent.bin_count, config);
    master.add_initial_patterns(make_initial_patterns(instance, incumbent));
    const ConflictGraph conflicts = build_conflict_graph(instance);

    for (int iteration = 0; iteration < config.max_cg_iterations; ++iteration) {
        if (deadline.expired()) {
            partial.timed_out = true;
            break;
        }
        ++statistics.cg_iterations;
        if (!master.solve(deadline.end_time(), statistics)) {
            partial.timed_out = master.timed_out() || deadline.expired();
            partial.numerical_failure = !partial.timed_out;
            break;
        }

        const PositionFreeMaster::Duals duals = master.duals();
        const double threshold = duals.pattern_count;
        const ScaledPositionFreeDuals certificate_scaled =
            scale_position_free_duals(duals.item, threshold, false);
        if (master.link_objective_multiplier_to_scale(
                certificate_scaled.scale)) {
            continue;
        }
        const ScaledPositionFreeDuals scaled = scale_position_free_duals(
            duals.item, threshold, false,
            std::min(certificate_scaled.scale,
                     pricing_scale_limit(master.objective_multiplier())));
        const std::uint64_t scale =
            static_cast<std::uint64_t>(scaled.scale);
        if (partial.fixed_point_scale_min == 0U) {
            partial.fixed_point_scale_min = scale;
        } else {
            partial.fixed_point_scale_min =
                std::min(partial.fixed_point_scale_min, scale);
        }
        partial.fixed_point_scale_max =
            std::max(partial.fixed_point_scale_max, scale);
        const std::uint64_t certificate_scale =
            static_cast<std::uint64_t>(certificate_scaled.scale);
        if (partial.certificate_scale_min == 0U) {
            partial.certificate_scale_min = certificate_scale;
        } else {
            partial.certificate_scale_min =
                std::min(partial.certificate_scale_min, certificate_scale);
        }
        LongPricingResult priced = run_long_pricing(
            instance, conflicts, master.pattern_keys(), scaled.item,
            scaled.threshold, deadline.end_time(), statistics);
        if (!priced.proven) {
            partial.timed_out = true;
            break;
        }
        if (!priced.patterns.empty()) {
            add_priced_patterns(master, priced, statistics);
            continue;
        }

        const IntegerPricingResult exact = run_integer_pricing(
            instance, conflicts, certificate_scaled.item, deadline.end_time(),
            statistics);
        if (!exact.proven) {
            partial.timed_out = true;
            break;
        }

        const std::int64_t denominator =
            std::max(certificate_scaled.scale, exact.maximum_profit);
        const int certified =
            ceil_ratio(certificate_scaled.sum, denominator);
        partial.certified_lower_bound = std::min(
            incumbent.bin_count,
            std::max(lower_bound, certified));
        partial.lp_value = master.objective_value();
        partial.completed = true;
        break;
    }

    if (!partial.completed && !partial.timed_out &&
        !partial.numerical_failure) {
        partial.timed_out = deadline.expired();
        partial.numerical_failure = !partial.timed_out;
    }
    ++statistics.cg_count;
    statistics.cg_seconds +=
        std::chrono::duration<double>(Clock::now() - start).count();
    return make_root_statistics(partial, master, statistics,
                                initial_statistics, start);
}

}
