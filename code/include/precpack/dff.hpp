#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace precpack {

struct DffTransformSet {
    int item_count = 0;
    std::vector<std::int64_t> capacities;
    std::vector<std::int64_t> contributions;

    [[nodiscard]] std::size_t size() const noexcept {
        return capacities.size();
    }

    [[nodiscard]] const std::int64_t* row(std::size_t transform) const noexcept {
        return contributions.data() + transform *
               static_cast<std::size_t>(item_count);
    }
};

[[nodiscard]] DffTransformSet build_complete_dff_transforms(
    const std::vector<int>& weights,
    int capacity,
    bool include_dff3_family = true);

[[nodiscard]] DffTransformSet select_ranked_complete_dff_transforms(
    const std::vector<int>& weights,
    int capacity,
    bool include_dff3_family,
    std::size_t maximum_transforms);

}
