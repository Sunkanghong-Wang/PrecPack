#include "precpack/algorithms.hpp"

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <utility>
#include <vector>

namespace precpack {

bool check_assignment(const Instance& instance,
                      const Assignment& assignment,
                      std::string* diagnostic) {
    const auto fail = [&](std::string message) {
        if (diagnostic != nullptr) {
            *diagnostic = std::move(message);
        }
        return false;
    };
    if (assignment.bin_of_item.size() != instance.items.size()) {
        return fail("assignment length differs from item count");
    }
    if (assignment.bin_count <= 0) {
        return fail("bin count is not positive");
    }

    std::vector<std::int64_t> loads(
        static_cast<std::size_t>(assignment.bin_count), 0);
    int maximum_bin = -1;
    for (int item = 0; item < instance.size(); ++item) {
        const int bin = assignment.bin_of_item[static_cast<std::size_t>(item)];
        if (bin < 0 || bin >= assignment.bin_count) {
            return fail("item " + std::to_string(item + 1) +
                        " has an invalid bin index");
        }
        maximum_bin = std::max(maximum_bin, bin);
        loads[static_cast<std::size_t>(bin)] +=
            instance.items[static_cast<std::size_t>(item)].weight;
    }
    if (maximum_bin + 1 != assignment.bin_count) {
        return fail("bin count is not the last occupied position plus one");
    }
    for (int bin = 0; bin < assignment.bin_count; ++bin) {
        if (loads[static_cast<std::size_t>(bin)] > instance.capacity) {
            return fail("capacity is exceeded in bin " +
                        std::to_string(bin));
        }
    }
    for (const Arc& arc : instance.arcs) {
        if (assignment.bin_of_item[static_cast<std::size_t>(arc.to)] -
                assignment.bin_of_item[static_cast<std::size_t>(arc.from)] <
            arc.separation) {
            std::ostringstream message;
            message << "precedence arc (" << arc.from + 1 << ','
                    << arc.to + 1 << ',' << arc.separation
                    << ") is violated";
            return fail(message.str());
        }
    }
    if (diagnostic != nullptr) {
        diagnostic->clear();
    }
    return true;
}

}
