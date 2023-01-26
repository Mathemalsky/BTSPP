#pragma once

#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"

namespace approximation {
struct Result {
  AdjacencyMatrixGraph<Directionality::Undirected> biconnectedGraph;
  OpenEarDecomposition openEarDecomposition;
  // const std::vector<unsigned int> tour;
  // const Edge bottleneckEdge;
};

Result approximate(const Euclidean& euclidean, const ProblemType problemType);
}  // namespace approximation
