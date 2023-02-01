#pragma once

#include <vector>

#include <solve/definitions.hpp>

#include "graph/graph.hpp"

namespace exactsolver {

struct Result {
  std::vector<unsigned int> tour;
  double opt;
  Edge bottleneckEdge;
};

Result solve(const Euclidean& euclidean, const ProblemType problemType);

}  // namespace exactsolver
