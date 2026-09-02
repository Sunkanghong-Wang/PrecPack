#include "precpack/dff.hpp"

#include "precpack/exact_arithmetic.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace precpack {
namespace {

struct Candidate {
    std::int64_t capacity = 1;
    std::vector<std::int64_t> values;
    std::int64_t total = 0;
    int root_bound = 0;
};

struct Threshold {
    std::int64_t numerator = 0;
    std::int64_t denominator = 1;
};

void validate_input(const std::vector<int>& weights, int capacity) {
    if (capacity <= 0 || weights.empty()) {
        throw std::invalid_argument("invalid DFF input");
    }
    for (const int weight : weights) {
        if (weight <= 0 || weight > capacity) {
            throw std::invalid_argument("DFF item weight is outside capacity");
        }
    }
}

void canonicalize(Candidate& candidate) {
    std::int64_t divisor = candidate.capacity;
    for (const std::int64_t value : candidate.values) {
        divisor = std::gcd(divisor, value);
    }
    if (divisor > 1) {
        candidate.capacity /= divisor;
        for (std::int64_t& value : candidate.values) {
            value /= divisor;
        }
    }
    candidate.total = 0;
    for (const std::int64_t value : candidate.values) {
        candidate.total = exact_arithmetic::checked_add(
            candidate.total, value, "DFF contribution total overflow");
    }
    candidate.root_bound = exact_arithmetic::ceil_ratio_to_int(
        candidate.total, candidate.capacity, "DFF root bound overflow");
}

[[nodiscard]] bool candidate_less(const Candidate& lhs,
                                  const Candidate& rhs) {
    if (lhs.capacity != rhs.capacity) {
        return lhs.capacity < rhs.capacity;
    }
    return lhs.values < rhs.values;
}

[[nodiscard]] bool candidate_equal(const Candidate& lhs,
                                   const Candidate& rhs) {
    return lhs.capacity == rhs.capacity && lhs.values == rhs.values;
}

[[nodiscard]] bool candidate_better(const Candidate& lhs,
                                    const Candidate& rhs) {
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
}

[[nodiscard]] bool ranked_candidate_less(const Candidate& lhs,
                                         const Candidate& rhs) {
    if (candidate_better(lhs, rhs)) {
        return true;
    }
    if (candidate_better(rhs, lhs)) {
        return false;
    }
    return candidate_less(lhs, rhs);
}

[[nodiscard]] std::vector<Threshold> build_thresholds(
    const std::vector<int>& weights,
    int capacity) {
    std::vector<Threshold> thresholds{{0, 1}, {1, 2}};
    thresholds.reserve(weights.size() + 2U);
    for (const int weight : weights) {
        if (2LL * weight <= capacity) {
            const std::int64_t divisor = std::gcd(weight, capacity);
            thresholds.push_back(
                Threshold{weight / divisor, capacity / divisor});
        }
    }
    const auto threshold_less = [](const Threshold& lhs,
                                   const Threshold& rhs) {
        return exact_arithmetic::compare_nonnegative_fractions(
                   static_cast<std::uint64_t>(lhs.numerator),
                   static_cast<std::uint64_t>(lhs.denominator),
                   static_cast<std::uint64_t>(rhs.numerator),
                   static_cast<std::uint64_t>(rhs.denominator)) < 0;
    };
    const auto threshold_equal = [](const Threshold& lhs,
                                    const Threshold& rhs) {
        return exact_arithmetic::compare_nonnegative_fractions(
                   static_cast<std::uint64_t>(lhs.numerator),
                   static_cast<std::uint64_t>(lhs.denominator),
                   static_cast<std::uint64_t>(rhs.numerator),
                   static_cast<std::uint64_t>(rhs.denominator)) == 0;
    };
    std::sort(thresholds.begin(), thresholds.end(), threshold_less);
    thresholds.erase(
        std::unique(thresholds.begin(), thresholds.end(), threshold_equal),
        thresholds.end());
    return thresholds;
}

template <typename Emit>
void generate_complete_dff_candidates(const std::vector<int>& weights,
                                      int capacity,
                                      bool include_dff3_family,
                                      Emit&& emit) {
    const std::vector<Threshold> thresholds =
        build_thresholds(weights, capacity);
    const auto emit_compositions = [&](std::int64_t base_capacity,
                                       const std::vector<std::int64_t>& base) {
        Candidate candidate;
        candidate.values.resize(weights.size(), 0);
        for (const Threshold threshold : thresholds) {
            candidate.capacity = base_capacity;
            for (std::size_t item = 0; item < weights.size(); ++item) {
                const std::int64_t value = base[item];
                const int upper_comparison =
                    exact_arithmetic::compare_nonnegative_fractions(
                        static_cast<std::uint64_t>(value),
                        static_cast<std::uint64_t>(base_capacity),
                        static_cast<std::uint64_t>(
                            threshold.denominator - threshold.numerator),
                        static_cast<std::uint64_t>(threshold.denominator));
                if (upper_comparison > 0) {
                    candidate.values[item] = base_capacity;
                } else if (exact_arithmetic::compare_nonnegative_fractions(
                               static_cast<std::uint64_t>(value),
                               static_cast<std::uint64_t>(base_capacity),
                               static_cast<std::uint64_t>(threshold.numerator),
                               static_cast<std::uint64_t>(
                                   threshold.denominator)) >= 0) {
                    candidate.values[item] = value;
                } else {
                    candidate.values[item] = 0;
                }
            }
            canonicalize(candidate);
            emit(candidate);
        }
    };

    std::vector<std::int64_t> base(weights.begin(), weights.end());
    emit_compositions(capacity, base);
    for (int parameter = 1; parameter <= 100; ++parameter) {
        const std::int64_t transformed_capacity =
            static_cast<std::int64_t>(parameter) * capacity;
        for (std::size_t item = 0; item < weights.size(); ++item) {
            const std::int64_t weight = weights[item];
            const std::int64_t numerator =
                static_cast<std::int64_t>(parameter + 1) * weight;
            base[item] = numerator % capacity == 0
                             ? static_cast<std::int64_t>(parameter) * weight
                             : (numerator / capacity) * capacity;
        }
        emit_compositions(transformed_capacity, base);
    }

    if (!include_dff3_family) {
        return;
    }
    Candidate candidate;
    candidate.values.resize(weights.size(), 0);
    for (int numerator = 1; numerator <= 500; ++numerator) {
        candidate.capacity = 1000 / numerator;
        std::fill(candidate.values.begin(), candidate.values.end(), 0);
        const std::int64_t scaled_parameter =
            static_cast<std::int64_t>(capacity) * numerator;
        for (std::size_t item = 0; item < weights.size(); ++item) {
            const std::int64_t weight = weights[item];
            if (2LL * weight > capacity) {
                candidate.values[item] = candidate.capacity -
                    (static_cast<std::int64_t>(capacity - weight) * 1000) /
                        scaled_parameter;
            } else if (static_cast<std::int64_t>(weight) * 1000 >=
                       scaled_parameter) {
                candidate.values[item] =
                    (static_cast<std::int64_t>(weight) * 1000) /
                    scaled_parameter;
            }
        }
        canonicalize(candidate);
        emit(candidate);
    }
}

[[nodiscard]] DffTransformSet make_result(
    const std::vector<Candidate>& candidates,
    std::size_t item_count) {
    DffTransformSet result;
    result.item_count = static_cast<int>(item_count);
    result.capacities.reserve(candidates.size());
    result.contributions.reserve(candidates.size() * item_count);
    for (const Candidate& candidate : candidates) {
        result.capacities.push_back(candidate.capacity);
        result.contributions.insert(result.contributions.end(),
                                    candidate.values.begin(),
                                    candidate.values.end());
    }
    return result;
}

}

DffTransformSet build_complete_dff_transforms(
    const std::vector<int>& weights,
    int capacity,
    bool include_dff3_family) {
    validate_input(weights, capacity);
    std::vector<Candidate> candidates;
    generate_complete_dff_candidates(
        weights, capacity, include_dff3_family,
        [&](const Candidate& candidate) {
            candidates.push_back(candidate);
        });
    std::sort(candidates.begin(), candidates.end(), candidate_less);
    candidates.erase(
        std::unique(candidates.begin(), candidates.end(), candidate_equal),
        candidates.end());
    return make_result(candidates, weights.size());
}

DffTransformSet select_ranked_complete_dff_transforms(
    const std::vector<int>& weights,
    int capacity,
    bool include_dff3_family,
    std::size_t maximum_transforms) {
    validate_input(weights, capacity);
    if (maximum_transforms == 0U) {
        DffTransformSet result;
        result.item_count = static_cast<int>(weights.size());
        return result;
    }
    std::vector<Candidate> selected;
    selected.reserve(maximum_transforms);
    generate_complete_dff_candidates(
        weights, capacity, include_dff3_family,
        [&](const Candidate& candidate) {
            if (std::any_of(selected.begin(), selected.end(),
                            [&](const Candidate& previous) {
                                return candidate_equal(previous, candidate);
                            })) {
                return;
            }
            const auto position = std::lower_bound(
                selected.begin(), selected.end(), candidate,
                ranked_candidate_less);
            if (selected.size() == maximum_transforms &&
                position == selected.end()) {
                return;
            }
            selected.insert(position, candidate);
            if (selected.size() > maximum_transforms) {
                selected.pop_back();
            }
        });
    return make_result(selected, weights.size());
}

}
