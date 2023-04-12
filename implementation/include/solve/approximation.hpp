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

/*!
 * \brief approximates the BTSP
 * \param euclidean complete graph, providing distances between nodes
 * \param printInfo controls if objective, lower bound on OPT and a fortiori guarantee are printed to console
 * \return Result
 */
Result approximateBTSP(const graph::Euclidean& euclidean, const bool printInfo = true);

/*!
 * @brief approximates the BTSPP
 * @param euclidean complete graph, providing distances between nodes
 * @param s start node
 * @param t end node
 * @param printInfo controls if objective, lower bound on OPT and a fortiori guarantee are printed to console
 * @return Result
 */
Result approximateBTSPP(const graph::Euclidean& euclidean, const size_t s = 0, const size_t t = 1, const bool printInfo = true);
}  // namespace approximation
