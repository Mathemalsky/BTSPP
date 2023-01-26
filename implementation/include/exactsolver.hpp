#pragma once

#include <vector>

#include <draw/definitions.hpp>

#include "graph/graph.hpp"

namespace exactsolver {

struct Result {
  const std::vector<unsigned int> tour;
  const double opt;
  const Edge bottleneckEdge;
};

Result solve(const Euclidean& euclidean, const ProblemType problemType);

}  // namespace exactsolver
