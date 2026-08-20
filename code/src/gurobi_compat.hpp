#pragma once

#include <gurobi_c++.h>

namespace precpack::gurobi_compat {

[[nodiscard]] inline bool is_work_limit_status(int status) noexcept {
#if GRB_VERSION_MAJOR > 9 || \
    (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 5)
    return status == GRB_WORK_LIMIT;
#else
    static_cast<void>(status);
    return false;
#endif
}

}
