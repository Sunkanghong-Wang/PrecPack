#pragma once

#ifndef PRECPACK_HAS_GUROBI
#define PRECPACK_HAS_GUROBI 0
#endif

namespace precpack {

inline constexpr bool kHasGurobiSupport = PRECPACK_HAS_GUROBI != 0;

}
