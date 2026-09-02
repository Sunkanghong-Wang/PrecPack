#pragma once

#include "precpack/types.hpp"

namespace precpack::internal {

[[nodiscard]] int compute_bppp_flow_lower_bound(const Instance& instance);

}
