#pragma once

#include <vector>

#include <draw/definitions.hpp>

#include "graph/graph.hpp"

namespace exactsolver {

struct Result {
  const std::vector<unsigned int> tour;
  const Edge bottleneck;
  const double opt;
};

Result solve(const Euclidean& euclidean, const ProblemType problemType);

}  // namespace exactsolver
