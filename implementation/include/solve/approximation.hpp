#pragma once

#include <vector>

#include "definitions.hpp"

#include "graph/graph.hpp"
#include "graph/algorithm.hpp"

namespace approximation {
struct Result {
  AdjacencyMatrixGraph biconnectedGraph;
  EarDecomposition openEarDecomposition;
  std::vector<unsigned int> tour;
  // const Edge bottleneckEdge;
};

Result approximate(const Euclidean& euclidean, const ProblemType problemType);
}  // namespace approximation
