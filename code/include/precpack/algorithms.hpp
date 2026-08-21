#pragma once

#include "precpack/types.hpp"

#include <string>

namespace precpack {

[[nodiscard]] bool check_assignment(
    const Instance& instance,
    const Assignment& assignment,
    std::string* diagnostic = nullptr);

}
