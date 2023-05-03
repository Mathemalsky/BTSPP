#pragma once

#include <vector>

#include <solve/definitions.hpp>

// graph library
#include "graph.hpp"

namespace exactsolver {

struct Result {
  std::vector<unsigned int> tour;
  double opt;
  graph::Edge bottleneckEdge;
};

/*!
 * @brief solves an instance of BTSP, BTSPP or TSP
 * @param problemType type of instance
 * @param noCrossing if BTSP the solution can be forced to have no crossings
 */
Result solve(const graph::Euclidean& euclidean, const ProblemType problemType, const bool noCrossing = false);

}  // namespace exactsolver
