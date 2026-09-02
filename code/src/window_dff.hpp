#pragma once

#include <vector>

namespace precpack::internal {

[[nodiscard]] int compute_window_dff_lower_bound(
    const std::vector<int>& weights,
    const std::vector<int>& front,
    const std::vector<int>& back,
    int capacity);

}
