#include "precpack/algorithms.hpp"
#include "precpack/bbr.hpp"
#include "precpack/initial_bounds.hpp"
#include "precpack/types.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

[[nodiscard]] precpack::Instance make_instance(
    const std::vector<int>& weights,
    int capacity,
    std::vector<precpack::Arc> arcs = {}) {
    precpack::Instance instance;
    instance.problem_type = "BPP-GP";
    instance.capacity = capacity;
    for (int item = 0; item < static_cast<int>(weights.size()); ++item) {
        instance.items.push_back({item, weights[static_cast<std::size_t>(item)]});
    }
    instance.arcs = std::move(arcs);
    instance.initialize();
    return instance;
}

[[nodiscard]] precpack::PreparedInstance make_prepared(
    const precpack::Instance& instance,
    precpack::Assignment incumbent) {
    precpack::PreparedInstance prepared;
    prepared.search_instance = instance;
    prepared.search_incumbent = std::move(incumbent);
    prepared.search_to_original.resize(static_cast<std::size_t>(instance.size()));
    std::iota(prepared.search_to_original.begin(),
              prepared.search_to_original.end(), 0);
    prepared.original_to_search = prepared.search_to_original;
    return prepared;
}

[[nodiscard]] precpack::Config base_config() {
    precpack::Config config;
    config.time_limit_seconds = 10.0;
    config.bbr_memory_limit_mb = 64;
    config.bbr_state_limit = 100'000;
    config.bbr_enable_complete_dff = false;
    config.bbr_enable_generalized_item_dominance = false;
    config.bbr_enable_jackson = false;
    config.bbr_enable_no_successor = false;
    config.bbr_enable_superset_memory = false;
    config.bbr_enable_profile_dominance = false;
    return config;
}

void test_parallel_exactness() {
    const precpack::Instance instance = make_instance(
        std::vector<int>(12U, 1), 3);
    precpack::Assignment incumbent;
    incumbent.bin_of_item.resize(12U);
    std::iota(incumbent.bin_of_item.begin(), incumbent.bin_of_item.end(), 0);
    incumbent.bin_count = 12;
    const precpack::PreparedInstance prepared =
        make_prepared(instance, incumbent);

    precpack::Config serial_config = base_config();
    serial_config.threads = 1;
    precpack::Deadline serial_deadline(10.0);
    const precpack::BbrResult serial = precpack::run_branch_bound_remember(
        prepared, 4, serial_config, serial_deadline);

    std::string diagnostic;
    for (const int threads : {2, 4, 8}) {
        precpack::Config parallel_config = base_config();
        parallel_config.threads = threads;
        precpack::Deadline parallel_deadline(10.0);
        const precpack::BbrResult parallel =
            precpack::run_branch_bound_remember(
                prepared, 4, parallel_config, parallel_deadline);
        require(serial.optimal && parallel.optimal &&
                    serial.incumbent.bin_count == 4 &&
                    parallel.incumbent.bin_count ==
                        serial.incumbent.bin_count &&
                    precpack::check_assignment(instance, parallel.incumbent,
                                                &diagnostic),
                "parallel BBR disagreed with the serial result at " +
                    std::to_string(threads) + " threads: " + diagnostic);
        require(
            parallel.statistics.parallel &&
                parallel.statistics.threads == threads &&
                parallel.statistics.parallel_tasks_generated > 1U &&
                parallel.statistics.parallel_tasks_completed ==
                    parallel.statistics.parallel_tasks_generated &&
                parallel.statistics.peak_memory_bytes <=
                    64U * 1024U * 1024U,
            "parallel task, termination, or memory contract failed at " +
                std::to_string(threads) + " threads");
    }
}

void test_generalized_parallel_exactness() {
    const precpack::Instance instance = make_instance(
        {142, 34, 140, 214, 121, 279, 50, 282, 129, 175,
         97, 132, 107, 132, 69, 169, 73, 231, 120, 186},
        1000,
        {{0, 5, 1},   {1, 6, 0},   {3, 7, 0},   {4, 8, 1},
         {5, 9, 0},   {6, 10, 0},  {7, 11, 0},  {9, 12, 0},
         {10, 12, 1}, {11, 13, 1}, {11, 14, 0}, {12, 15, 1},
         {12, 16, 0}, {12, 17, 0}, {13, 19, 1}, {14, 18, 1}});
    precpack::Assignment incumbent;
    incumbent.bin_of_item.resize(20U);
    std::iota(incumbent.bin_of_item.begin(), incumbent.bin_of_item.end(), 0);
    incumbent.bin_count = 20;
    const precpack::PreparedInstance prepared =
        make_prepared(instance, incumbent);

    precpack::Config serial_config = base_config();
    serial_config.threads = 1;
    precpack::Deadline serial_deadline(10.0);
    const precpack::BbrResult serial = precpack::run_branch_bound_remember(
        prepared, 3, serial_config, serial_deadline);

    std::string diagnostic;
    for (const int threads : {2, 4, 8}) {
        precpack::Config parallel_config = base_config();
        parallel_config.threads = threads;
        precpack::Deadline parallel_deadline(10.0);
        const precpack::BbrResult parallel =
            precpack::run_branch_bound_remember(
                prepared, 3, parallel_config, parallel_deadline);
        require(serial.optimal && parallel.optimal &&
                    parallel.incumbent.bin_count ==
                        serial.incumbent.bin_count &&
                    precpack::check_assignment(instance, parallel.incumbent,
                                                &diagnostic),
                "generalized parallel BBR disagreed with the serial result at " +
                    std::to_string(threads) + " threads: " + diagnostic);
        require(parallel.statistics.parallel &&
                    parallel.statistics.threads == threads &&
                    parallel.statistics.parallel_tasks_completed ==
                        parallel.statistics.parallel_tasks_generated &&
                    parallel.statistics.states_created > 1U,
                "generalized parallel task accounting failed at " +
                    std::to_string(threads) + " threads: generated=" +
                    std::to_string(
                        parallel.statistics.parallel_tasks_generated) +
                    ", completed=" +
                    std::to_string(
                        parallel.statistics.parallel_tasks_completed) +
                    ", states=" +
                    std::to_string(parallel.statistics.states_created));
    }
}

void test_global_limits() {
    const precpack::Instance instance = make_instance(
        std::vector<int>(12U, 1), 3);
    precpack::Assignment incumbent;
    incumbent.bin_of_item.resize(12U);
    std::iota(incumbent.bin_of_item.begin(), incumbent.bin_of_item.end(), 0);
    incumbent.bin_count = 12;
    const precpack::PreparedInstance prepared =
        make_prepared(instance, incumbent);
    precpack::Config config = base_config();
    config.threads = 2;

    precpack::Config state_config = config;
    state_config.bbr_state_limit = 1U;
    precpack::Deadline state_deadline(10.0);
    const precpack::BbrResult state_limited =
        precpack::run_branch_bound_remember(
            prepared, 4, state_config, state_deadline);
    require(!state_limited.optimal && state_limited.state_limited &&
                !state_limited.memory_limited &&
                state_limited.statistics.states_created <= 1U,
            "the global state cap was treated as a per-worker budget");

    precpack::Config memory_config = config;
    memory_config.bbr_memory_limit_mb = 4U;
    precpack::Deadline memory_deadline(10.0);
    const precpack::BbrResult memory_limited =
        precpack::run_branch_bound_remember(
            prepared, 4, memory_config, memory_deadline);
    require(!memory_limited.optimal && memory_limited.memory_limited &&
                !memory_limited.state_limited &&
                memory_limited.statistics.peak_memory_bytes <=
                    4U * 1024U * 1024U,
            "the global memory cap was exceeded or misclassified");

    precpack::Deadline timeout_deadline(1e-12);
    const precpack::BbrResult timed_out =
        precpack::run_branch_bound_remember(
            prepared, 4, config, timeout_deadline);
    require(!timed_out.optimal && timed_out.timed_out &&
                !timed_out.state_limited && !timed_out.memory_limited,
            "the global deadline was misclassified");
}

}

int main() {
    try {
        test_parallel_exactness();
        test_generalized_parallel_exactness();
        test_global_limits();
        std::cout << "Parallel BBR tests passed (threads=1,2,4,8).\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Parallel BBR test failed: " << error.what() << '\n';
        return 1;
    }
}
