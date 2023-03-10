#pragma once

#include <vector>

#include "definitions.hpp"

#include "graph/algorithm.hpp"
#include "graph/graph.hpp"

namespace approximation {
struct Result {
  graph::AdjacencyMatrixGraph biconnectedGraph;
  graph::EarDecomposition openEarDecomposition;
  std::vector<unsigned int> tour;
  double opt;
  graph::Edge bottleneckEdge;
};

Result approximate(const graph::Euclidean& euclidean, const ProblemType problemType, const bool printInfo = true);
}  // namespace approximation
