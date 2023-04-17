#pragma once

#include <vector>

#include <solve/definitions.hpp>

#include "graph/graph.hpp"

namespace exactsolver {

struct Result {
  std::vector<unsigned int> tour;
  double opt;
  graph::Edge bottleneckEdge;
};

Result solve(const graph::Euclidean& euclidean, const ProblemType problemType, const bool noCrossing = false);

}  // namespace exactsolver
