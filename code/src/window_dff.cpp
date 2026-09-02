#include "window_dff.hpp"

#include "precpack/exact_arithmetic.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace precpack::internal {
namespace {

class WindowDffEvaluator {
public:
    WindowDffEvaluator(const std::vector<int>& weights,
                       const std::vector<int>& front,
                       const std::vector<int>& back,
                       int capacity)
        : weights_(weights),
          front_(front),
          back_(back),
          capacity_(capacity) {
        if (capacity_ <= 0 || weights_.size() != front_.size() ||
            weights_.size() != back_.size()) {
            throw std::invalid_argument("invalid window-DFF input");
        }
        for (std::size_t item = 0; item < weights_.size(); ++item) {
            if (weights_[item] <= 0 || weights_[item] > capacity_ ||
                front_[item] < 0 || back_[item] < 0) {
                throw std::invalid_argument(
                    "window-DFF item data is outside the supported range");
            }
            longest_ = std::max(longest_, front_[item]);
        }
        if (longest_ < 2) {
            return;
        }
        stride_ = static_cast<std::size_t>(longest_) + 2U;
        if (stride_ > std::numeric_limits<std::size_t>::max() / stride_) {
            throw std::length_error("window-DFF table is too large");
        }
        sums_.assign(stride_ * stride_, 0);
        parameter_max_back_.assign(stride_, -1);
        lower_bound_ = longest_ + 1;
    }

    [[nodiscard]] int solve() {
        if (longest_ < 2) {
            return lower_bound_;
        }

        for (int parameter = 1; parameter <= 100; ++parameter) {
            const std::int64_t denominator =
                static_cast<std::int64_t>(capacity_) * parameter;
            evaluate(denominator, [this, parameter](int weight) {
                const std::int64_t scaled =
                    static_cast<std::int64_t>(parameter + 1) * weight;
                return scaled % capacity_ == 0
                           ? static_cast<std::int64_t>(weight) * parameter
                           : (scaled / capacity_) * capacity_;
            });
        }

        evaluate(capacity_, [this](int weight) -> std::int64_t {
            if (2LL * weight > capacity_) {
                return capacity_;
            }
            return 2LL * weight == capacity_ ? weight : 0;
        });

        std::vector<int> thresholds = weights_;
        std::sort(thresholds.begin(), thresholds.end());
        thresholds.erase(std::unique(thresholds.begin(), thresholds.end()),
                         thresholds.end());
        for (const int threshold : thresholds) {
            if (2LL * threshold >= capacity_) {
                continue;
            }
            evaluate(capacity_, [this, threshold](int weight) -> std::int64_t {
                if (weight > capacity_ - threshold) {
                    return capacity_;
                }
                return weight >= threshold ? weight : 0;
            },
                     threshold);
        }

        for (int numerator = 1; numerator <= 500; ++numerator) {
            const std::int64_t transformed_capacity = 1000 / numerator;
            const std::int64_t ratio_denominator =
                static_cast<std::int64_t>(capacity_) * numerator;
            evaluate(transformed_capacity,
                     [this, ratio_denominator,
                      transformed_capacity](int weight) -> std::int64_t {
                         if (2LL * weight > capacity_) {
                             const std::int64_t complement_units =
                                 (static_cast<std::int64_t>(capacity_ - weight) *
                                  1000) /
                                 ratio_denominator;
                             return transformed_capacity - complement_units;
                         }
                         if (static_cast<std::int64_t>(weight) * 1000 >=
                             ratio_denominator) {
                             return (static_cast<std::int64_t>(weight) * 1000) /
                                    ratio_denominator;
                         }
                         return 0;
                     });
        }
        return lower_bound_;
    }

private:
    template <typename Contribution>
    void evaluate(std::int64_t denominator,
                  Contribution&& contribution,
                  int required_parameter_weight = -1) {
        std::fill(sums_.begin(), sums_.end(), 0);
        if (required_parameter_weight >= 0) {
            std::fill(parameter_max_back_.begin(),
                      parameter_max_back_.end(), -1);
        }

        for (std::size_t item = 0; item < weights_.size(); ++item) {
            const int front = std::min(front_[item], longest_);
            const int back = std::min(back_[item], longest_);
            const std::size_t index =
                static_cast<std::size_t>(front) * stride_ +
                static_cast<std::size_t>(back);
            sums_[index] = exact_arithmetic::checked_add(
                sums_[index], contribution(weights_[item]),
                "window-DFF contribution sum overflow");
            if (weights_[item] == required_parameter_weight) {
                parameter_max_back_[static_cast<std::size_t>(front)] =
                    std::max(parameter_max_back_[static_cast<std::size_t>(front)],
                             back);
            }
        }

        for (int front = longest_; front >= 0; --front) {
            for (int back = longest_; back >= 0; --back) {
                const std::size_t index =
                    static_cast<std::size_t>(front) * stride_ +
                    static_cast<std::size_t>(back);
                const std::int64_t down = sums_[index + stride_];
                const std::int64_t diagonal = sums_[index + stride_ + 1U];
                const std::int64_t right = sums_[index + 1U];
                const std::int64_t down_only =
                    exact_arithmetic::checked_subtract(
                        down, diagonal,
                        "window-DFF suffix sum became negative");
                sums_[index] = exact_arithmetic::checked_add(
                    sums_[index],
                    exact_arithmetic::checked_add(
                        down_only, right,
                        "window-DFF suffix sum overflow"),
                    "window-DFF suffix sum overflow");
            }
        }

        if (required_parameter_weight >= 0) {
            for (int front = longest_ - 1; front >= 0; --front) {
                parameter_max_back_[static_cast<std::size_t>(front)] = std::max(
                    parameter_max_back_[static_cast<std::size_t>(front)],
                    parameter_max_back_[static_cast<std::size_t>(front + 1)]);
            }
        }

        for (int front = 1; front <= longest_; ++front) {
            const int last_back = longest_ - front;
            for (int back = 1; back <= last_back; ++back) {
                if (required_parameter_weight >= 0 &&
                    parameter_max_back_[static_cast<std::size_t>(front)] < back) {
                    continue;
                }
                const std::int64_t sum =
                    sums_[static_cast<std::size_t>(front) * stride_ +
                          static_cast<std::size_t>(back)];
                const int dff_bound = exact_arithmetic::ceil_ratio_to_int(
                    sum, denominator, "window-DFF lower bound overflow");
                const std::int64_t candidate =
                    static_cast<std::int64_t>(front) + back + dff_bound;
                if (candidate > std::numeric_limits<int>::max()) {
                    throw std::overflow_error(
                        "window-DFF lower bound exceeds the supported range");
                }
                lower_bound_ =
                    std::max(lower_bound_, static_cast<int>(candidate));
            }
        }
    }

    const std::vector<int>& weights_;
    const std::vector<int>& front_;
    const std::vector<int>& back_;
    int capacity_ = 0;
    int longest_ = 0;
    int lower_bound_ = 1;
    std::size_t stride_ = 0;
    std::vector<std::int64_t> sums_;
    std::vector<int> parameter_max_back_;
};

}

int compute_window_dff_lower_bound(const std::vector<int>& weights,
                                   const std::vector<int>& front,
                                   const std::vector<int>& back,
                                   int capacity) {
    return WindowDffEvaluator(weights, front, back, capacity).solve();
}

}
