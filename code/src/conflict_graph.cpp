#include "precpack/conflict_graph.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace precpack {

std::vector<std::pair<int, int>> build_bppc_relaxation_conflict_edges(
    const Instance& instance) {
    const int n = instance.size();
    const std::size_t blocks =
        (static_cast<std::size_t>(n) + 63U) / 64U;
    std::vector<std::uint64_t> reachable(
        static_cast<std::size_t>(n) * blocks, 0U);
    const auto set_reachable = [&](int from, int to) {
        reachable[static_cast<std::size_t>(from) * blocks +
                  static_cast<std::size_t>(to) / 64U] |=
            std::uint64_t{1} << (static_cast<unsigned>(to) & 63U);
    };
    for (int item = 0; item < n; ++item) {
        set_reachable(item, item);
    }
    for (auto iterator = instance.topological_order.rbegin();
         iterator != instance.topological_order.rend(); ++iterator) {
        const int from = *iterator;
        std::uint64_t* row = reachable.data() +
                             static_cast<std::size_t>(from) * blocks;
        for (const int to :
             instance.successors[static_cast<std::size_t>(from)]) {
            const std::uint64_t* successor_row =
                reachable.data() + static_cast<std::size_t>(to) * blocks;
            for (std::size_t block = 0; block < blocks; ++block) {
                row[block] |= successor_row[block];
            }
        }
    }

    std::vector<std::uint64_t> reaching(
        static_cast<std::size_t>(n) * blocks, 0U);
    for (int from = 0; from < n; ++from) {
        const std::uint64_t* row = reachable.data() +
                                   static_cast<std::size_t>(from) * blocks;
        for (std::size_t block = 0; block < blocks; ++block) {
            std::uint64_t bits = row[block];
            while (bits != 0U) {
                const unsigned bit = std::countr_zero(bits);
                const int to = static_cast<int>(block * 64U + bit);
                if (to < n) {
                    reaching[static_cast<std::size_t>(to) * blocks +
                             static_cast<std::size_t>(from) / 64U] |=
                        std::uint64_t{1} <<
                        (static_cast<unsigned>(from) & 63U);
                }
                bits &= bits - 1U;
            }
        }
    }

    const auto is_reachable = [&](int from, int to) noexcept {
        return ((reachable[static_cast<std::size_t>(from) * blocks +
                           static_cast<std::size_t>(to) / 64U] >>
                 (static_cast<unsigned>(to) & 63U)) &
                std::uint64_t{1}) != 0U;
    };
    const auto interval_overloads = [&](int from, int to) {
        std::int64_t load = 0;
        const std::uint64_t* descendants =
            reachable.data() + static_cast<std::size_t>(from) * blocks;
        const std::uint64_t* ancestors =
            reaching.data() + static_cast<std::size_t>(to) * blocks;
        for (std::size_t block = 0; block < blocks; ++block) {
            std::uint64_t bits = descendants[block] & ancestors[block];
            while (bits != 0U) {
                const unsigned bit = std::countr_zero(bits);
                const int item = static_cast<int>(block * 64U + bit);
                if (item < n) {
                    load += instance.items[static_cast<std::size_t>(item)].weight;
                    if (load > instance.capacity) {
                        return true;
                    }
                }
                bits &= bits - 1U;
            }
        }
        return false;
    };

    std::vector<std::pair<int, int>> edges;
    for (int lhs = 0; lhs < n; ++lhs) {
        for (int rhs = lhs + 1; rhs < n; ++rhs) {
            if (instance.separation(lhs, rhs) > 0 ||
                instance.separation(rhs, lhs) > 0) {
                edges.emplace_back(lhs, rhs);
                continue;
            }
            if ((is_reachable(lhs, rhs) && interval_overloads(lhs, rhs)) ||
                (is_reachable(rhs, lhs) && interval_overloads(rhs, lhs))) {
                edges.emplace_back(lhs, rhs);
            }
        }
    }
    return edges;
}

}
