#pragma once

#include "precpack/types.hpp"

#include <utility>
#include <vector>

namespace precpack {

[[nodiscard]] std::vector<std::pair<int, int>>
build_bppc_relaxation_conflict_edges(const Instance& instance);

}
