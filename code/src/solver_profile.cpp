#include "precpack/solver_profile.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <thread>

namespace precpack {

ProblemKind parse_problem_kind(std::string_view value) {
    if (value == "salbp-i") {
        return ProblemKind::kSalbpI;
    }
    if (value == "bpp-p") {
        return ProblemKind::kBppP;
    }
    if (value == "bpp-gp") {
        return ProblemKind::kBppGp;
    }
    throw std::invalid_argument(
        "problem must be salbp-i, bpp-p, or bpp-gp: " +
        std::string(value));
}

const char* to_string(ProblemKind problem) noexcept {
    switch (problem) {
        case ProblemKind::kSalbpI:
            return "SALBP-I";
        case ProblemKind::kBppP:
            return "BPP-P";
        case ProblemKind::kBppGp:
            return "BPP-GP";
    }
    return "UNKNOWN";
}

const char* to_slug(ProblemKind problem) noexcept {
    switch (problem) {
        case ProblemKind::kSalbpI:
            return "salbp-i";
        case ProblemKind::kBppP:
            return "bpp-p";
        case ProblemKind::kBppGp:
            return "bpp-gp";
    }
    return "unknown";
}

Config make_solver_config(ProblemKind problem,
                          double time_limit_seconds,
                          std::uint64_t memory_limit_mb,
                          int threads) {
    if (!std::isfinite(time_limit_seconds) || time_limit_seconds <= 0.0) {
        throw std::invalid_argument("time limit must be positive and finite");
    }
    if (memory_limit_mb == 0U) {
        throw std::invalid_argument("memory limit must be positive");
    }
    if (threads == 0 || threads < -1) {
        throw std::invalid_argument("threads must be -1 or a positive integer");
    }

    Config config;
    config.time_limit_seconds = time_limit_seconds;
    config.seed = 1;
    config.threads = threads;
    config.bbr_memory_limit_mb = memory_limit_mb;
    config.bbr_enable_early_exact_probe = true;
    config.bbr_enable_jackson = true;
    config.bbr_enable_no_successor = true;
    config.bbr_enable_superset_memory = true;
    config.bbr_enable_profile_dominance = true;
    config.bbr_enable_structured_preprocessing = true;
    config.bbr_root_cg_time_limit_seconds = 5.0;

    if (problem == ProblemKind::kSalbpI) {
        config.bbr_enable_initial_bdp = false;
        config.bbr_enable_bbr12_mhh = true;
        config.bbr_enable_bbr12_mhh_portfolio = true;
        config.bbr12_mhh_portfolio_max_items = 200;
        config.bbr12_mhh_full_load_limit = 1'000;
        config.bbr_enable_generalized_item_dominance = false;
        config.bbr_enable_paper_queue_order = false;
        config.bbr_enable_bbr12_jackson = true;
        config.bbr_enable_complete_dff = true;
        config.bbr_dff_transform_limit = 15;
        config.bbr_enable_closure_bound = false;
        config.bbr_enable_binlb = true;
        config.bbr_binlb_call_time_limit_seconds = 1.0;
        config.bbr_binlb_total_time_limit_seconds = 0.1;
        config.bbr_binlb_node_limit = 1'000'000U;
        config.bbr_binlb_load_limit = 50U;
        config.bbr_binlb_memo_limit = 200'000U;
        config.bbr_binlb_max_items = 400;
        config.bbr_enable_root_strengthening = false;
    } else {
        config.bbr_enable_initial_bdp = true;
        config.bbr_enable_bbr12_mhh = false;
        config.bbr_enable_bbr12_mhh_portfolio = false;
        config.bbr_enable_generalized_item_dominance = true;
        config.bbr_enable_paper_queue_order = true;
        config.bbr_enable_bbr12_jackson = false;
        config.bbr_enable_complete_dff = true;
        config.bbr_dff_transform_limit = 0;
        config.bbr_enable_closure_bound = true;
        config.bbr_enable_binlb = false;
        config.bbr_enable_root_strengthening = true;
    }
    return config;
}

int resolve_thread_count(int requested_threads) noexcept {
    if (requested_threads != -1) {
        return std::max(1, requested_threads);
    }
    const unsigned available = std::thread::hardware_concurrency();
    return available == 0U ? 1 : static_cast<int>(available);
}

}
