#include "precpack/dff.hpp"

#include "precpack/exact_arithmetic.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace precpack {
namespace {

struct Candidate {
    std::int64_t capacity = 1;
    std::vector<std::int64_t> values;
};

struct Threshold {
    std::int64_t numerator = 0;
    std::int64_t denominator = 1;
};

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

}

DffTransformSet build_complete_dff_transforms(
    const std::vector<int>& weights,
    int capacity,
    bool include_dff3_family) {
    if (capacity <= 0 || weights.empty()) {
        throw std::invalid_argument("invalid DFF input");
    }
    for (const int weight : weights) {
        if (weight <= 0 || weight > capacity) {
            throw std::invalid_argument("DFF item weight is outside capacity");
        }
    }

    std::vector<Threshold> thresholds{{0, 1}, {1, 2}};
    thresholds.reserve(weights.size() + 2U);
    for (const int weight : weights) {
        if (2LL * weight <= capacity) {
            const std::int64_t divisor = std::gcd(weight, capacity);
            thresholds.push_back(Threshold{weight / divisor,
                                           capacity / divisor});
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

    std::vector<Candidate> candidates;
    candidates.reserve(101U * thresholds.size() +
                       (include_dff3_family ? 500U : 0U));
    const auto append_compositions = [&](std::int64_t base_capacity,
                                         const std::vector<std::int64_t>& base) {
        for (const Threshold threshold : thresholds) {
            Candidate candidate;
            candidate.capacity = base_capacity;
            candidate.values.resize(weights.size(), 0);
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
                }
            }
            canonicalize(candidate);
            candidates.push_back(std::move(candidate));
        }
    };

    std::vector<std::int64_t> base(weights.begin(), weights.end());
    append_compositions(capacity, base);
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
        append_compositions(transformed_capacity, base);
    }

    if (include_dff3_family) {
        for (int numerator = 1; numerator <= 500; ++numerator) {
            Candidate candidate;
            candidate.capacity = 1000 / numerator;
            candidate.values.resize(weights.size(), 0);
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
            candidates.push_back(std::move(candidate));
        }
    }

    std::sort(candidates.begin(), candidates.end(), candidate_less);
    candidates.erase(
        std::unique(candidates.begin(), candidates.end(), candidate_equal),
        candidates.end());

    DffTransformSet result;
    result.item_count = static_cast<int>(weights.size());
    result.capacities.reserve(candidates.size());
    result.contributions.reserve(candidates.size() * weights.size());
    for (const Candidate& candidate : candidates) {
        result.capacities.push_back(candidate.capacity);
        result.contributions.insert(result.contributions.end(),
                                    candidate.values.begin(),
                                    candidate.values.end());
    }
    return result;
}

}
