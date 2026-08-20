#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace precpack::exact_arithmetic {

[[nodiscard]] inline std::int64_t checked_add(
    std::int64_t lhs,
    std::int64_t rhs,
    const char* message) {
    if ((rhs > 0 && lhs > std::numeric_limits<std::int64_t>::max() - rhs) ||
        (rhs < 0 && lhs < std::numeric_limits<std::int64_t>::min() - rhs)) {
        throw std::overflow_error(message);
    }
    return lhs + rhs;
}

[[nodiscard]] inline std::int64_t checked_subtract(
    std::int64_t lhs,
    std::int64_t rhs,
    const char* message) {
    if ((rhs > 0 && lhs < std::numeric_limits<std::int64_t>::min() + rhs) ||
        (rhs < 0 && lhs > std::numeric_limits<std::int64_t>::max() + rhs)) {
        throw std::overflow_error(message);
    }
    return lhs - rhs;
}

[[nodiscard]] inline std::int64_t checked_multiply(
    std::int64_t lhs,
    std::int64_t rhs,
    const char* message) {
    if (lhs == 0 || rhs == 0) {
        return 0;
    }
    if ((lhs == -1 && rhs == std::numeric_limits<std::int64_t>::min()) ||
        (rhs == -1 && lhs == std::numeric_limits<std::int64_t>::min())) {
        throw std::overflow_error(message);
    }
    if (lhs > 0) {
        if ((rhs > 0 &&
             lhs > std::numeric_limits<std::int64_t>::max() / rhs) ||
            (rhs < 0 &&
             rhs < std::numeric_limits<std::int64_t>::min() / lhs)) {
            throw std::overflow_error(message);
        }
    } else if ((rhs > 0 &&
                lhs < std::numeric_limits<std::int64_t>::min() / rhs) ||
               (rhs < 0 &&
                lhs < std::numeric_limits<std::int64_t>::max() / rhs)) {
        throw std::overflow_error(message);
    }
    return lhs * rhs;
}

[[nodiscard]] inline int ceil_ratio_to_int(std::int64_t numerator,
                                           std::int64_t denominator,
                                           const char* message) {
    if (denominator <= 0) {
        throw std::invalid_argument("ratio denominator must be positive");
    }
    std::int64_t quotient = numerator / denominator;
    if (numerator % denominator > 0) {
        quotient = checked_add(quotient, 1, message);
    }
    if (quotient > std::numeric_limits<int>::max() ||
        quotient < std::numeric_limits<int>::min()) {
        throw std::overflow_error(message);
    }
    return static_cast<int>(quotient);
}

[[nodiscard]] inline std::int64_t ceil_nonnegative_product_ratio(
    std::int64_t value,
    std::int64_t multiplier,
    std::int64_t denominator,
    const char* message) {
    if (value < 0 || multiplier < 0 || denominator <= 0) {
        throw std::invalid_argument(
            "product ratio requires nonnegative terms");
    }
    const std::int64_t whole = checked_multiply(
        value / denominator, multiplier, message);
    const std::int64_t remainder_product = checked_multiply(
        value % denominator, multiplier, message);
    std::int64_t fractional = remainder_product / denominator;
    if (remainder_product % denominator != 0) {
        fractional = checked_add(fractional, 1, message);
    }
    return checked_add(whole, fractional, message);
}

[[nodiscard]] inline int compare_nonnegative_fractions(
    std::uint64_t lhs_numerator,
    std::uint64_t lhs_denominator,
    std::uint64_t rhs_numerator,
    std::uint64_t rhs_denominator) {
    if (lhs_denominator == 0U || rhs_denominator == 0U) {
        throw std::invalid_argument("fraction denominator must be positive");
    }

    constexpr std::uint64_t kMaximum =
        std::numeric_limits<std::uint64_t>::max();
    if (lhs_numerator <= kMaximum / rhs_denominator &&
        rhs_numerator <= kMaximum / lhs_denominator) {
        const std::uint64_t lhs = lhs_numerator * rhs_denominator;
        const std::uint64_t rhs = rhs_numerator * lhs_denominator;
        return lhs < rhs ? -1 : (lhs > rhs ? 1 : 0);
    }

    bool reversed = false;
    for (;;) {
        const std::uint64_t lhs_quotient =
            lhs_numerator / lhs_denominator;
        const std::uint64_t rhs_quotient =
            rhs_numerator / rhs_denominator;
        if (lhs_quotient != rhs_quotient) {
            const int comparison = lhs_quotient < rhs_quotient ? -1 : 1;
            return reversed ? -comparison : comparison;
        }

        const std::uint64_t lhs_remainder =
            lhs_numerator % lhs_denominator;
        const std::uint64_t rhs_remainder =
            rhs_numerator % rhs_denominator;
        if (lhs_remainder == 0U || rhs_remainder == 0U) {
            const int comparison =
                lhs_remainder == rhs_remainder
                    ? 0
                    : (lhs_remainder == 0U ? -1 : 1);
            return reversed ? -comparison : comparison;
        }

        lhs_numerator = lhs_denominator;
        lhs_denominator = lhs_remainder;
        rhs_numerator = rhs_denominator;
        rhs_denominator = rhs_remainder;
        reversed = !reversed;
    }
}

class NonnegativeRatioSum {
public:
    explicit NonnegativeRatioSum(std::int64_t denominator)
        : denominator_(denominator) {
        if (denominator_ <= 0) {
            throw std::invalid_argument("ratio denominator must be positive");
        }
    }

    void add(std::int64_t numerator) {
        if (numerator < 0) {
            throw std::invalid_argument("ratio numerator must be nonnegative");
        }
        increase_quotient(numerator / denominator_);
        const std::int64_t added_remainder = numerator % denominator_;
        if (added_remainder == 0) {
            return;
        }
        if (remainder_ >= denominator_ - added_remainder) {
            increase_quotient(1);
            remainder_ -= denominator_ - added_remainder;
        } else {
            remainder_ += added_remainder;
        }
    }

    [[nodiscard]] int ceil_to_int() const {
        std::int64_t result = quotient_;
        if (remainder_ != 0) {
            if (result == std::numeric_limits<int>::max()) {
                throw std::overflow_error("ratio result does not fit int");
            }
            ++result;
        }
        return static_cast<int>(result);
    }

private:
    void increase_quotient(std::int64_t value) {
        if (value > std::numeric_limits<int>::max() - quotient_) {
            throw std::overflow_error("ratio result does not fit int");
        }
        quotient_ += value;
    }

    std::int64_t denominator_ = 1;
    std::int64_t quotient_ = 0;
    std::int64_t remainder_ = 0;
};

}
