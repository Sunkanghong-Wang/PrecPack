#include "precpack/column_generation.hpp"

#include <gurobi_c++.h>

#include <algorithm>
#include <array>
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
using Int128 = __int128_t;

struct Pattern {
    int bin = -1;
    std::vector<int> items;
    std::vector<std::uint64_t> bits;

    [[nodiscard]] bool contains(int item) const noexcept {
        return ((bits[static_cast<std::size_t>(item) / 64U] >>
                 (static_cast<unsigned>(item) & 63U)) &
                1U) != 0U;
    }
};

struct SrCut {
    std::array<int, 3> items{};
    double dual = 0.0;
};

[[nodiscard]] std::string pattern_key(int bin,
                                      const std::vector<std::uint64_t>& bits) {
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

[[nodiscard]] int sr_coefficient(const Pattern& pattern,
                                 const std::array<int, 3>& triple) noexcept {
    const int count = static_cast<int>(pattern.contains(triple[0])) +
                      static_cast<int>(pattern.contains(triple[1])) +
                      static_cast<int>(pattern.contains(triple[2]));
    return count >= 2 ? 1 : 0;
}

class InitialBppMaster {
public:
    InitialBppMaster(GRBEnv& environment,
                     const Instance& instance,
                     int lower_bound,
                     const Config& config)
        : instance_(instance), model_(environment) {
        model_.set(GRB_IntParam_Threads, 1);
        model_.set(GRB_IntParam_Seed, static_cast<int>(config.seed));
        model_.set(GRB_IntParam_OutputFlag, config.gurobi_log ? 1 : 0);
        model_.set(GRB_IntParam_Method, 1);
        model_.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
        model_.set(GRB_DoubleParam_OptimalityTol, 1e-9);
        model_.set(GRB_IntParam_NumericFocus, 2);

        const std::size_t initial_capacity =
            static_cast<std::size_t>(std::max(1024, instance.size() * 4));
        patterns_.reserve(initial_capacity);
        variables_.reserve(initial_capacity);
        variable_values_.reserve(initial_capacity);
        pattern_keys_.reserve(initial_capacity * 2U);
        item_rows_.reserve(static_cast<std::size_t>(instance.size()));
        for (int item = 0; item < instance.size(); ++item) {
            item_rows_.push_back(model_.addConstr(GRBLinExpr(0.0) >= 1.0));
        }
        valid_bound_row_ =
            model_.addConstr(GRBLinExpr(0.0) >= static_cast<double>(lower_bound));
        sr_cuts_.reserve(1000U);
        sr_rows_.reserve(1000U);
    }

    void add_initial_patterns(const Assignment& assignment) {
        const std::size_t blocks =
            (static_cast<std::size_t>(instance_.size()) + 63U) / 64U;
        std::vector<Pattern> patterns(static_cast<std::size_t>(assignment.bin_count));
        for (Pattern& pattern : patterns) {
            pattern.bin = 0;
            pattern.bits.assign(blocks, 0U);
        }
        for (int item = 0; item < instance_.size(); ++item) {
            Pattern& pattern = patterns[static_cast<std::size_t>(
                assignment.bin_of_item[static_cast<std::size_t>(item)])];
            pattern.items.push_back(item);
            pattern.bits[static_cast<std::size_t>(item) / 64U] |=
                std::uint64_t{1} << (static_cast<unsigned>(item) & 63U);
        }
        for (Pattern& pattern : patterns) {
            if (!pattern.items.empty()) {
                add_pattern(std::move(pattern));
            }
        }
        model_.update();
    }

    void add_pattern(Pattern pattern) {
        const std::string key = pattern_key(0, pattern.bits);
        if (!pattern_keys_.insert(key).second) {
            return;
        }
        GRBColumn column;
        for (const int item : pattern.items) {
            column.addTerm(1.0, item_rows_[static_cast<std::size_t>(item)]);
        }
        column.addTerm(1.0, valid_bound_row_);
        for (std::size_t cut = 0; cut < sr_cuts_.size(); ++cut) {
            if (sr_coefficient(pattern, sr_cuts_[cut].items) != 0) {
                column.addTerm(1.0, sr_rows_[cut]);
            }
        }
        variables_.push_back(model_.addVar(0.0, GRB_INFINITY, 1.0,
                                           GRB_CONTINUOUS, column));
        patterns_.push_back(std::move(pattern));
    }

    void add_sr_cut(const std::array<int, 3>& triple) {
        GRBLinExpr expression = 0.0;
        for (std::size_t pattern = 0; pattern < patterns_.size(); ++pattern) {
            if (sr_coefficient(patterns_[pattern], triple) != 0) {
                expression += variables_[pattern];
            }
        }
        sr_cuts_.push_back(SrCut{triple, 0.0});
        sr_rows_.push_back(model_.addConstr(expression <= 1.0));
        sr_keys_.insert(triple_key(triple[0], triple[1], triple[2]));
        model_.update();
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
        if (model_.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
            return false;
        }
        objective_value_ = model_.get(GRB_DoubleAttr_ObjVal);
        variable_values_.resize(variables_.size());
        for (std::size_t variable = 0; variable < variables_.size(); ++variable) {
            variable_values_[variable] = variables_[variable].get(GRB_DoubleAttr_X);
        }
        return true;
    }

    struct Duals {
        std::vector<double> item;
        std::vector<double> sr;
        double valid_bound = 0.0;
    };

    [[nodiscard]] Duals duals() {
        Duals result;
        result.item.reserve(item_rows_.size());
        for (GRBConstr& row : item_rows_) {
            result.item.push_back(std::max(0.0, row.get(GRB_DoubleAttr_Pi)));
        }
        result.sr.reserve(sr_rows_.size());
        for (std::size_t cut = 0; cut < sr_rows_.size(); ++cut) {
            const double dual = std::min(0.0, sr_rows_[cut].get(GRB_DoubleAttr_Pi));
            sr_cuts_[cut].dual = dual;
            result.sr.push_back(dual);
        }
        result.valid_bound = std::max(0.0, valid_bound_row_.get(GRB_DoubleAttr_Pi));
        return result;
    }

    [[nodiscard]] std::vector<std::array<int, 3>> violated_sr_cuts(
        int maximum_new_cuts,
        double tolerance) const {
        struct Candidate {
            double violation = 0.0;
            std::array<int, 3> triple{};
        };
        std::vector<int> positive_patterns;
        positive_patterns.reserve(patterns_.size());
        for (std::size_t pattern = 0; pattern < patterns_.size(); ++pattern) {
            if (variable_values_[pattern] >= tolerance) {
                positive_patterns.push_back(static_cast<int>(pattern));
            }
        }
        std::vector<Candidate> candidates;
        const int n = instance_.size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (instance_.items[static_cast<std::size_t>(i)].weight +
                        instance_.items[static_cast<std::size_t>(j)].weight >
                    instance_.capacity) {
                    continue;
                }
                double pair_value = 0.0;
                for (const int pattern : positive_patterns) {
                    if (patterns_[static_cast<std::size_t>(pattern)].contains(i) &&
                        patterns_[static_cast<std::size_t>(pattern)].contains(j)) {
                        pair_value += variable_values_[static_cast<std::size_t>(pattern)];
                    }
                }
                if (pair_value < tolerance) {
                    continue;
                }
                for (int k = j + 1; k < n; ++k) {
                    if (instance_.items[static_cast<std::size_t>(i)].weight +
                                instance_.items[static_cast<std::size_t>(k)].weight >
                            instance_.capacity ||
                        instance_.items[static_cast<std::size_t>(j)].weight +
                                instance_.items[static_cast<std::size_t>(k)].weight >
                            instance_.capacity ||
                        sr_keys_.contains(triple_key(i, j, k))) {
                        continue;
                    }
                    double lhs = pair_value;
                    for (const int pattern : positive_patterns) {
                        const Pattern& column = patterns_[static_cast<std::size_t>(pattern)];
                        if (!column.contains(k)) {
                            continue;
                        }
                        const int pair_count = static_cast<int>(column.contains(i)) +
                                               static_cast<int>(column.contains(j));
                        if (pair_count == 1) {
                            lhs += variable_values_[static_cast<std::size_t>(pattern)];
                        }
                    }
                    if (lhs > 1.0 + tolerance) {
                        candidates.push_back(Candidate{lhs - 1.0, {i, j, k}});
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
        std::vector<std::array<int, 3>> result;
        result.reserve(candidates.size());
        for (const Candidate& candidate : candidates) {
            result.push_back(candidate.triple);
        }
        return result;
    }

    [[nodiscard]] double objective_value() const noexcept { return objective_value_; }
    [[nodiscard]] int sr_cut_count() const noexcept {
        return static_cast<int>(sr_cuts_.size());
    }
    [[nodiscard]] const std::vector<SrCut>& sr_cuts() const noexcept { return sr_cuts_; }
    [[nodiscard]] const std::unordered_set<std::string>& pattern_keys() const noexcept {
        return pattern_keys_;
    }

private:
    const Instance& instance_;
    GRBModel model_;
    std::vector<GRBConstr> item_rows_;
    GRBConstr valid_bound_row_;
    std::vector<SrCut> sr_cuts_;
    std::vector<GRBConstr> sr_rows_;
    std::unordered_set<std::uint64_t> sr_keys_;
    std::vector<Pattern> patterns_;
    std::vector<GRBVar> variables_;
    std::vector<double> variable_values_;
    std::unordered_set<std::string> pattern_keys_;
    double objective_value_ = 0.0;
};

[[nodiscard]] int certified_initial_bpp_dual_bound(
    const Instance& instance,
    const InitialBppMaster::Duals& duals,
    int valid_bound_rhs) {
    long double magnitude = 1.0L;
    for (const double value : duals.item) {
        magnitude += std::abs(static_cast<long double>(value));
    }
    magnitude += static_cast<long double>(valid_bound_rhs) *
                 std::abs(static_cast<long double>(duals.valid_bound));
    for (const double value : duals.sr) {
        magnitude += std::abs(static_cast<long double>(value));
    }
    magnitude = std::max(1.0L, magnitude);

    std::int64_t scale = std::int64_t{1} << 40U;
    constexpr long double kLimit =
        static_cast<long double>(std::int64_t{1} << 60U);
    while (scale > 1 && magnitude * static_cast<long double>(scale) > kLimit) {
        scale /= 2;
    }
    const auto scaled = [scale](double value) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("non-finite initial BPP dual");
        }
        const long double product =
            static_cast<long double>(value) * static_cast<long double>(scale);
        if (product >
                static_cast<long double>(std::numeric_limits<std::int64_t>::max()) ||
            product <
                static_cast<long double>(std::numeric_limits<std::int64_t>::min())) {
            throw std::overflow_error("initial BPP dual scaling overflow");
        }
        return static_cast<std::int64_t>(std::llround(product));
    };

    std::vector<std::int64_t> item_profit;
    item_profit.reserve(duals.item.size());
    Int128 numerator = 0;
    for (const double value : duals.item) {
        const std::int64_t profit = std::max<std::int64_t>(0, scaled(value));
        item_profit.push_back(profit);
        numerator += profit;
    }
    const std::int64_t valid_dual =
        std::max<std::int64_t>(0, scaled(duals.valid_bound));
    numerator += static_cast<Int128>(valid_bound_rhs) * valid_dual;
    for (const double value : duals.sr) {
        numerator += std::min<std::int64_t>(0, scaled(value));
    }

    std::vector<std::int64_t> dp(
        static_cast<std::size_t>(instance.capacity + 1), 0);
    for (int item = 0; item < instance.size(); ++item) {
        const int weight =
            instance.items[static_cast<std::size_t>(item)].weight;
        const std::int64_t profit = item_profit[static_cast<std::size_t>(item)];
        for (int capacity = instance.capacity; capacity >= weight; --capacity) {
            const Int128 candidate =
                static_cast<Int128>(dp[static_cast<std::size_t>(capacity - weight)]) +
                profit;
            if (candidate > std::numeric_limits<std::int64_t>::max()) {
                throw std::overflow_error("initial BPP knapsack certificate overflow");
            }
            dp[static_cast<std::size_t>(capacity)] = std::max(
                dp[static_cast<std::size_t>(capacity)],
                static_cast<std::int64_t>(candidate));
        }
    }
    const Int128 maximum_column_score_wide =
        static_cast<Int128>(valid_dual) +
        *std::max_element(dp.begin(), dp.end());
    if (maximum_column_score_wide >
        std::numeric_limits<std::int64_t>::max()) {
        throw std::overflow_error("initial BPP column-score overflow");
    }
    const std::int64_t maximum_column_score =
        static_cast<std::int64_t>(maximum_column_score_wide);
    const std::int64_t denominator =
        std::max(scale, maximum_column_score);
    Int128 quotient = 0;
    if (numerator >= 0) {
        quotient = (numerator + denominator - 1) / denominator;
    } else {
        quotient = numerator / denominator;
    }
    if (quotient > std::numeric_limits<int>::max() ||
        quotient < std::numeric_limits<int>::min()) {
        throw std::overflow_error("initial BPP certified bound does not fit int");
    }
    return static_cast<int>(quotient);
}

struct ConflictGraph {
    std::vector<std::vector<int>> neighbors;
};

[[nodiscard]] ConflictGraph build_direct_conflict_graph(const Instance& instance) {
    ConflictGraph graph;
    graph.neighbors.assign(static_cast<std::size_t>(instance.size()), {});
    for (const Arc& arc : instance.arcs) {
        if (arc.separation <= 0) {
            continue;
        }
        graph.neighbors[static_cast<std::size_t>(arc.from)].push_back(arc.to);
        graph.neighbors[static_cast<std::size_t>(arc.to)].push_back(arc.from);
    }
    return graph;
}

class PricingSearch {
public:
    struct Result {
        std::optional<Pattern> pattern;
        double reduced_cost = 0.0;
        bool proven = true;
        std::uint64_t nodes = 0;
    };

    PricingSearch(const Instance& instance,
                  const ConflictGraph& conflicts,
                  const std::vector<SrCut>& sr_cuts,
                  const std::vector<double>& sr_duals,
                  const std::unordered_set<std::string>& existing_patterns,
                  int bin,
                  int bin_count,
                  std::vector<double> item_profit,
                  double base_cost,
                  double tolerance,
                  Clock::time_point end_time,
                  bool respect_positions = true,
                  int minimum_load = 0)
        : instance_(instance),
          conflicts_(conflicts),
          sr_cuts_(sr_cuts),
          sr_duals_(sr_duals),
          existing_patterns_(existing_patterns),
          bin_(bin),
          base_cost_(base_cost),
          tolerance_(tolerance),
          end_time_(end_time),
          minimum_load_(minimum_load),
          blocked_(static_cast<std::size_t>(instance.size()), 0),
          cut_count_(sr_cuts.size(), 0),
          selected_bits_((static_cast<std::size_t>(instance.size()) + 63U) / 64U,
                         0U),
          item_cut_indices_(static_cast<std::size_t>(instance.size())) {
        candidates_.reserve(static_cast<std::size_t>(instance.size()));
        selected_items_.reserve(static_cast<std::size_t>(instance.size()));
        best_items_.reserve(static_cast<std::size_t>(instance.size()));
        lookup_key_.resize(sizeof(bin_) +
                           selected_bits_.size() * sizeof(std::uint64_t));
        std::memcpy(lookup_key_.data(), &bin_, sizeof(bin_));
        for (std::size_t cut = 0; cut < sr_cuts_.size(); ++cut) {
            for (const int item : sr_cuts_[cut].items) {
                item_cut_indices_[static_cast<std::size_t>(item)].push_back(
                    static_cast<int>(cut));
            }
        }
        for (int item = 0; item < instance.size(); ++item) {
            bool participates_in_active_cut = false;
            for (const int cut : item_cut_indices_[static_cast<std::size_t>(item)]) {
                if (std::abs(sr_duals_[static_cast<std::size_t>(cut)]) > kEpsilon) {
                    participates_in_active_cut = true;
                    break;
                }
            }
            if ((respect_positions &&
                 (bin < instance.front[static_cast<std::size_t>(item)] ||
                  bin + instance.back[static_cast<std::size_t>(item)] >= bin_count)) ||
                (item_profit[static_cast<std::size_t>(item)] <= 0.0 &&
                 !participates_in_active_cut)) {
                continue;
            }
            candidates_.push_back(Candidate{
                item,
                instance.items[static_cast<std::size_t>(item)].weight,
                item_profit[static_cast<std::size_t>(item)]});
        }
        std::sort(candidates_.begin(), candidates_.end(), [](const Candidate& lhs,
                                                             const Candidate& rhs) {
            const long double lhs_cross =
                static_cast<long double>(lhs.profit) * rhs.weight;
            const long double rhs_cross =
                static_cast<long double>(rhs.profit) * lhs.weight;
            if (lhs_cross != rhs_cross) {
                return lhs_cross > rhs_cross;
            }
            if (lhs.profit != rhs.profit) {
                return lhs.profit > rhs.profit;
            }
            return lhs.item < rhs.item;
        });
        best_profit_ = base_cost_ + tolerance_;
    }

    [[nodiscard]] Result solve() {
        if (Clock::now() >= end_time_) {
            timed_out_ = true;
        }
        consider(0.0);
        if (!timed_out_) {
            search(0, instance_.capacity, 0.0);
        }
        Result result;
        result.proven = !timed_out_;
        result.nodes = nodes_;
        if (have_best_) {
            Pattern pattern;
            pattern.bin = bin_;
            pattern.items = best_items_;
            std::sort(pattern.items.begin(), pattern.items.end());
            pattern.bits = best_bits_;
            result.reduced_cost = base_cost_ - best_profit_;
            result.pattern = std::move(pattern);
        }
        return result;
    }

private:
    struct Candidate {
        int item = -1;
        int weight = 0;
        double profit = 0.0;
    };

    [[nodiscard]] double fractional_upper_bound(std::size_t position,
                                                int capacity,
                                                double profit) const {
        double bound = profit;
        for (std::size_t i = position; i < candidates_.size() && capacity > 0; ++i) {
            const Candidate& candidate = candidates_[i];
            if (blocked_[static_cast<std::size_t>(candidate.item)] != 0) {
                continue;
            }
            if (candidate.weight <= capacity) {
                capacity -= candidate.weight;
                bound += candidate.profit;
            } else {
                bound += candidate.profit *
                         (static_cast<double>(capacity) / candidate.weight);
                break;
            }
        }
        return bound;
    }

    void consider(double profit) {
        if (profit <= best_profit_ || selected_weight_ < minimum_load_) {
            return;
        }
        if (!selected_bits_.empty()) {
            std::memcpy(lookup_key_.data() + sizeof(bin_), selected_bits_.data(),
                        selected_bits_.size() * sizeof(std::uint64_t));
        }
        if (existing_patterns_.contains(lookup_key_)) {
            return;
        }
        best_profit_ = profit;
        best_items_ = selected_items_;
        best_bits_ = selected_bits_;
        have_best_ = true;
    }

    void include_item(const Candidate& candidate, double& profit) {
        selected_items_.push_back(candidate.item);
        selected_weight_ += candidate.weight;
        selected_bits_[static_cast<std::size_t>(candidate.item) / 64U] |=
            std::uint64_t{1} << (static_cast<unsigned>(candidate.item) & 63U);
        profit += candidate.profit;
        for (const int cut : item_cut_indices_[static_cast<std::size_t>(candidate.item)]) {
            if (cut_count_[static_cast<std::size_t>(cut)] == 1) {
                profit += sr_duals_[static_cast<std::size_t>(cut)];
            }
            ++cut_count_[static_cast<std::size_t>(cut)];
        }
        for (const int neighbor :
             conflicts_.neighbors[static_cast<std::size_t>(candidate.item)]) {
            ++blocked_[static_cast<std::size_t>(neighbor)];
        }
    }

    void remove_item(const Candidate& candidate, double& profit) {
        for (const int neighbor :
             conflicts_.neighbors[static_cast<std::size_t>(candidate.item)]) {
            --blocked_[static_cast<std::size_t>(neighbor)];
        }
        for (const int cut : item_cut_indices_[static_cast<std::size_t>(candidate.item)]) {
            --cut_count_[static_cast<std::size_t>(cut)];
            if (cut_count_[static_cast<std::size_t>(cut)] == 1) {
                profit -= sr_duals_[static_cast<std::size_t>(cut)];
            }
        }
        profit -= candidate.profit;
        selected_weight_ -= candidate.weight;
        selected_bits_[static_cast<std::size_t>(candidate.item) / 64U] &=
            ~(std::uint64_t{1} << (static_cast<unsigned>(candidate.item) & 63U));
        selected_items_.pop_back();
    }

    void search(std::size_t position, int capacity, double profit) {
        if (timed_out_ || position >= candidates_.size()) {
            return;
        }
        ++nodes_;
        if ((nodes_ & 2047U) == 0U && Clock::now() >= end_time_) {
            timed_out_ = true;
            return;
        }
        const double upper_bound =
            fractional_upper_bound(position, capacity, profit);
        const double safety =
            1e-10 * std::max(1.0, std::abs(upper_bound));
        if (upper_bound + safety <= best_profit_) {
            return;
        }

        const Candidate& candidate = candidates_[position];
        if (candidate.weight <= capacity &&
            blocked_[static_cast<std::size_t>(candidate.item)] == 0) {
            double included_profit = profit;
            include_item(candidate, included_profit);
            consider(included_profit);
            search(position + 1, capacity - candidate.weight, included_profit);
            remove_item(candidate, included_profit);
        }
        search(position + 1, capacity, profit);
    }

    const Instance& instance_;
    const ConflictGraph& conflicts_;
    const std::vector<SrCut>& sr_cuts_;
    const std::vector<double>& sr_duals_;
    const std::unordered_set<std::string>& existing_patterns_;
    int bin_;
    double base_cost_;
    double tolerance_;
    Clock::time_point end_time_;
    int minimum_load_ = 0;
    std::vector<Candidate> candidates_;
    std::vector<int> blocked_;
    std::vector<int> cut_count_;
    std::vector<std::uint64_t> selected_bits_;
    std::vector<int> selected_items_;
    std::vector<std::vector<int>> item_cut_indices_;
    double best_profit_ = 0.0;
    bool have_best_ = false;
    bool timed_out_ = false;
    std::uint64_t nodes_ = 0;
    int selected_weight_ = 0;
    std::vector<int> best_items_;
    std::vector<std::uint64_t> best_bits_;
    std::string lookup_key_;
};

}

ColumnGenerationResult run_initial_bpp_column_generation(
    GRBEnv& environment,
    const Instance& strengthened_instance,
    const Assignment& incumbent,
    int lower_bound,
    const Config& config,
    Deadline& global_deadline,
    Statistics& statistics) {
    ColumnGenerationResult result;
    result.attempted = true;
    result.integer_lower_bound = lower_bound;
    result.lp_value = static_cast<double>(lower_bound);
    if (global_deadline.expired() || incumbent.bin_count <= lower_bound) {
        result.converged = incumbent.bin_count <= lower_bound;
        result.pricing_proven = result.converged;
        return result;
    }

    const auto cg_start = Clock::now();
    const Clock::time_point cg_end = global_deadline.end_time();
    InitialBppMaster master(environment, strengthened_instance, lower_bound, config);
    master.add_initial_patterns(incumbent);
    const ConflictGraph conflicts =
        build_direct_conflict_graph(strengthened_instance);

    const std::int64_t incumbent_capacity =
        static_cast<std::int64_t>(incumbent.bin_count) *
        strengthened_instance.capacity;
    const int maximum_waste = static_cast<int>(std::clamp<std::int64_t>(
        incumbent_capacity - strengthened_instance.total_weight, 0,
        strengthened_instance.capacity));
    const int minimum_load = strengthened_instance.capacity - maximum_waste;

    constexpr int kMaximumInitialSrCuts = 1000;
    int valid_lower_bound = lower_bound;
    bool terminated = false;
    for (int iteration = 0;
         iteration < config.max_cg_iterations && !terminated; ++iteration) {
        if (Clock::now() >= cg_end) {
            result.timed_out = true;
            break;
        }

        bool exact_pricing = false;
        bool have_solved_master = false;
        bool stopped_by_valid_bound = false;
        while (Clock::now() < cg_end) {
            if (!master.solve(cg_end, statistics)) {
                result.timed_out = Clock::now() >= cg_end;
                terminated = true;
                break;
            }
            have_solved_master = true;
            result.lp_value = master.objective_value();
            InitialBppMaster::Duals duals = master.duals();
            result.integer_lower_bound = std::max(
                result.integer_lower_bound,
                certified_initial_bpp_dual_bound(
                    strengthened_instance, duals, lower_bound));
            result.integer_lower_bound =
                std::min(result.integer_lower_bound, incumbent.bin_count);
            if (result.integer_lower_bound >= incumbent.bin_count) {
                result.pricing_proven = true;
                terminated = true;
                break;
            }

            const auto pricing_start = Clock::now();
            PricingSearch pricing(
                strengthened_instance, conflicts, master.sr_cuts(), duals.sr,
                master.pattern_keys(), 0, 1, std::move(duals.item),
                1.0 - duals.valid_bound, config.reduced_cost_tolerance, cg_end,
                false, minimum_load);
            PricingSearch::Result priced = pricing.solve();
            statistics.pricing_seconds +=
                std::chrono::duration<double>(Clock::now() - pricing_start).count();
            ++statistics.pricing_count;
            statistics.pricing_search_nodes += priced.nodes;

            if (!priced.proven) {
                result.timed_out = true;
                terminated = true;
                break;
            }
            if (!priced.pattern.has_value() ||
                priced.reduced_cost >= -config.reduced_cost_tolerance) {
                exact_pricing = true;
                const double margin =
                    1e-5 + static_cast<double>(master.sr_cut_count() + 2) *
                                 config.reduced_cost_tolerance;
                valid_lower_bound = std::max(
                    valid_lower_bound,
                    static_cast<int>(std::ceil(result.lp_value - margin)));
                break;
            }

            const double denominator = 1.0 - priced.reduced_cost;
            if (denominator > 0.0) {
                const double valid_fractional_bound = result.lp_value / denominator;
                const double margin =
                    1e-5 + static_cast<double>(master.sr_cut_count() + 2) *
                                 config.reduced_cost_tolerance;
                valid_lower_bound = std::max(
                    valid_lower_bound,
                    static_cast<int>(std::ceil(valid_fractional_bound - margin)));
            }
            valid_lower_bound =
                std::min(valid_lower_bound, incumbent.bin_count);

            const double objective_margin =
                1e-5 + static_cast<double>(master.sr_cut_count() + 2) *
                             config.reduced_cost_tolerance;
            const int rounded_rmp =
                static_cast<int>(std::ceil(result.lp_value - objective_margin));
            if (valid_lower_bound >= incumbent.bin_count ||
                valid_lower_bound >= rounded_rmp) {
                stopped_by_valid_bound = true;
                break;
            }

            master.add_pattern(std::move(*priced.pattern));
            ++statistics.generated_columns;
        }

        if (terminated || result.integer_lower_bound >= incumbent.bin_count) {
            break;
        }
        if (!have_solved_master) {
            result.timed_out = Clock::now() >= cg_end;
            break;
        }

        if (config.enable_sr_cuts &&
            master.sr_cut_count() < kMaximumInitialSrCuts) {
            auto cuts = master.violated_sr_cuts(
                kMaximumInitialSrCuts - master.sr_cut_count(),
                config.row_violation_tolerance);
            if (!cuts.empty()) {
                for (const auto& cut : cuts) {
                    master.add_sr_cut(cut);
                }
                statistics.generated_sr_rows += cuts.size();
                continue;
            }
        }

        result.converged = exact_pricing;
        result.pricing_proven = exact_pricing;
        static_cast<void>(stopped_by_valid_bound);
        break;
    }

    if (!result.converged && !result.timed_out &&
        result.integer_lower_bound >= incumbent.bin_count) {
        result.pricing_proven = true;
    }
    statistics.cg_seconds +=
        std::chrono::duration<double>(Clock::now() - cg_start).count();
    ++statistics.cg_count;
    return result;
}

}
