#include "gurobi_oracle.hpp"

#include <gurobi_c++.h>

#include "gurobi_compat.hpp"
#include "precpack/algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace precpack::test {

MipResult solve_compact_mip(GRBEnv& environment,
                            const Instance& instance,
                            const Assignment& incumbent,
                            int lower_bound,
                            const Config& config,
                            Deadline& deadline) {
    if (!incumbent.complete()) {
        throw std::invalid_argument("compact MIP requires a complete incumbent");
    }
    std::string diagnostic;
    if (!check_assignment(instance, incumbent, &diagnostic)) {
        throw std::invalid_argument("invalid MIP start: " + diagnostic);
    }

    MipResult result;
    result.assignment = incumbent;
    result.has_incumbent = true;
    result.certified_lower_bound = lower_bound;
    if (lower_bound >= incumbent.bin_count) {
        result.status = OracleStatus::kOptimal;
        result.optimal = true;
        result.best_bound = static_cast<double>(lower_bound);
        return result;
    }
    if (deadline.expired()) {
        result.status = OracleStatus::kTimeLimit;
        result.best_bound = static_cast<double>(lower_bound);
        return result;
    }

    const int n = instance.size();
    const int m = incumbent.bin_count;
    GRBModel model(environment);
    model.set(GRB_IntParam_Threads, 1);
    model.set(GRB_IntParam_Seed, static_cast<int>(config.seed));
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_DoubleParam_MIPGap, 0.0);
    model.set(GRB_DoubleParam_MIPGapAbs, 0.0);
    model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
    model.set(GRB_DoubleParam_OptimalityTol, 1e-9);
    model.set(GRB_DoubleParam_IntFeasTol, 1e-9);
    model.set(GRB_IntParam_NumericFocus, 2);
    model.set(GRB_IntParam_IntegralityFocus, 1);

    std::vector<GRBVar> y;
    y.reserve(static_cast<std::size_t>(m));
    for (int bin = 0; bin < m; ++bin) {
        y.push_back(model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
    }

    std::vector<int> variable_index(static_cast<std::size_t>(n) * m, -1);
    std::vector<GRBVar> x;
    x.reserve(static_cast<std::size_t>(n) * m);
    for (int item = 0; item < n; ++item) {
        const int first = instance.front[static_cast<std::size_t>(item)];
        const int last = m - 1 - instance.back[static_cast<std::size_t>(item)];
        if (first > last) {
            throw std::logic_error("incumbent contradicts precedence position bounds");
        }
        for (int bin = first; bin <= last; ++bin) {
            variable_index[static_cast<std::size_t>(item) * m + bin] =
                static_cast<int>(x.size());
            x.push_back(model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
        }
    }

    GRBLinExpr objective = 0.0;
    for (const GRBVar& variable : y) {
        objective += variable;
    }
    model.setObjective(objective, GRB_MINIMIZE);
    model.addConstr(objective >= lower_bound);

    for (int item = 0; item < n; ++item) {
        GRBLinExpr assigned_once = 0.0;
        for (int bin = 0; bin < m; ++bin) {
            const int index =
                variable_index[static_cast<std::size_t>(item) * m + bin];
            if (index >= 0) {
                assigned_once += x[static_cast<std::size_t>(index)];
            }
        }
        model.addConstr(assigned_once == 1.0);
    }

    for (int bin = 0; bin < m; ++bin) {
        GRBLinExpr load = 0.0;
        for (int item = 0; item < n; ++item) {
            const int index =
                variable_index[static_cast<std::size_t>(item) * m + bin];
            if (index >= 0) {
                load += instance.items[static_cast<std::size_t>(item)].weight *
                        x[static_cast<std::size_t>(index)];
            }
        }
        model.addConstr(load <= instance.capacity * y[static_cast<std::size_t>(bin)]);
    }

    for (const Arc& arc : instance.arcs) {
        GRBLinExpr distance = 0.0;
        for (int bin = 0; bin < m; ++bin) {
            const int from_index =
                variable_index[static_cast<std::size_t>(arc.from) * m + bin];
            const int to_index =
                variable_index[static_cast<std::size_t>(arc.to) * m + bin];
            if (to_index >= 0) {
                distance += bin * x[static_cast<std::size_t>(to_index)];
            }
            if (from_index >= 0) {
                distance -= bin * x[static_cast<std::size_t>(from_index)];
            }
        }
        model.addConstr(distance >= arc.separation);
    }

    for (int bin = 0; bin + 1 < m; ++bin) {
        model.addConstr(y[static_cast<std::size_t>(bin)] >=
                        y[static_cast<std::size_t>(bin + 1)]);
    }

    for (int bin = 0; bin < m; ++bin) {
        y[static_cast<std::size_t>(bin)].set(GRB_DoubleAttr_Start, 1.0);
    }
    for (int item = 0; item < n; ++item) {
        for (int bin = 0; bin < m; ++bin) {
            const int index =
                variable_index[static_cast<std::size_t>(item) * m + bin];
            if (index >= 0) {
                x[static_cast<std::size_t>(index)].set(
                    GRB_DoubleAttr_Start,
                    incumbent.bin_of_item[static_cast<std::size_t>(item)] == bin ? 1.0
                                                                                 : 0.0);
            }
        }
    }

    if (deadline.expired()) {
        result.status = OracleStatus::kTimeLimit;
        result.best_bound = static_cast<double>(lower_bound);
        return result;
    }
    model.set(GRB_DoubleParam_TimeLimit, deadline.remaining_seconds());
    model.optimize();
    const int status = model.get(GRB_IntAttr_Status);
    const int solution_count = model.get(GRB_IntAttr_SolCount);
    result.explored_nodes = static_cast<std::uint64_t>(
        std::max(0.0, model.get(GRB_DoubleAttr_NodeCount)));

    double best_bound = static_cast<double>(lower_bound);
    int certified_lower_bound = lower_bound;
    try {
        const double gurobi_bound = model.get(GRB_DoubleAttr_ObjBound);
        if (std::isfinite(gurobi_bound)) {
            best_bound = std::max(best_bound, gurobi_bound);
            const double integral_bound = std::floor(gurobi_bound + 0.5);
            if (integral_bound >= static_cast<double>(lower_bound) &&
                integral_bound <= static_cast<double>(incumbent.bin_count) &&
                integral_bound <=
                    static_cast<double>(std::numeric_limits<int>::max())) {
                certified_lower_bound = std::max(
                    certified_lower_bound, static_cast<int>(integral_bound));
            }
        }
    } catch (const GRBException&) {
    }
    result.best_bound = best_bound;
    result.certified_lower_bound = certified_lower_bound;

    if (solution_count > 0) {
        Assignment candidate;
        candidate.bin_of_item.assign(static_cast<std::size_t>(n), -1);
        int last_bin = -1;
        for (int item = 0; item < n; ++item) {
            for (int bin = 0; bin < m; ++bin) {
                const int index =
                    variable_index[static_cast<std::size_t>(item) * m + bin];
                if (index >= 0 &&
                    x[static_cast<std::size_t>(index)].get(GRB_DoubleAttr_X) > 0.5) {
                    candidate.bin_of_item[static_cast<std::size_t>(item)] = bin;
                    last_bin = std::max(last_bin, bin);
                    break;
                }
            }
        }
        candidate.bin_count = last_bin + 1;
        if (!check_assignment(instance, candidate, &diagnostic)) {
            throw std::logic_error("Gurobi returned an invalid incumbent: " + diagnostic);
        }
        if (candidate.bin_count <= result.assignment.bin_count) {
            result.assignment = std::move(candidate);
        }
        result.has_incumbent = true;
    }

    if (status == GRB_OPTIMAL) {
        result.status = OracleStatus::kOptimal;
        result.optimal = true;
        result.certified_lower_bound = result.assignment.bin_count;
        result.best_bound = static_cast<double>(result.assignment.bin_count);
    } else if (status == GRB_TIME_LIMIT || status == GRB_INTERRUPTED ||
               gurobi_compat::is_work_limit_status(status) ||
               status == GRB_NODE_LIMIT) {
        result.status = OracleStatus::kTimeLimit;
    } else if (status == GRB_INFEASIBLE) {
        result.status = OracleStatus::kError;
        throw std::logic_error("compact MIP rejected a validated incumbent");
    } else if (result.has_incumbent) {
        result.status = OracleStatus::kFeasible;
    } else {
        result.status = OracleStatus::kError;
    }
    return result;
}

}
