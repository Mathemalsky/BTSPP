/*
 * pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
 * BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <vector>

// graph library
#include "algorithm.hpp"
#include "graph.hpp"

#include "solve/commonfunctions.hpp"
#include "solve/definitions.hpp"

namespace approximation {
/*!
 * @brief Result bundles all important measures from the approximation
 */
struct Result {
  graph::AdjacencyListGraph biconnectedGraph;     /**< bottleneck optimal biconnected subgraph */
  graph::EarDecomposition openEarDecomposition;   /**< open ear decomposition */
  std::vector<size_t> tour;                       /**< hamilton cycle in square of original graph */
  graph::Edge bottleneckEdge;                     /**< a longest edge in the tour */
  double objective;                               /**< length of the longest edge */
  double lowerBoundOnOPT;                         /**< lower bound on opt */
  size_t numberOfEdgesInMinimallyBiconectedGraph; /**< number of edges in th minimally biconnected subgraph */
};

/*!
 * @brief lists all important information from solving in the terminal
 * @param res result
 * @param problemType type of instance
 * @param runtime elapsed time during approximation, -1.0 is default as invalid value and indicates the function to not print the runtime
 */
void printInfo(const approximation::Result& res, const ProblemType problemType, const double runtime = -1.0);

/*!
 * @brief remove all edges whcih are not 2-essential
 * @details The ear decomposition is computed to cheaply get rid of many edges at once. The removal has roughly the same computaional costs
 * as every check (using schmidt algorithm) for biconnectivity after removal of a single edge.
 * @tparam G type of graph
 * @param biconnectedGraph a biconnected graph as input
 * @return AdjacencyListGraph: minimally biconnected graph
 */
template <typename G>
  requires(std::is_base_of_v<graph::Graph, G>)
static graph::AdjacencyListGraph makeMinimallyBiconnected(const G& biconnectedGraph) {
  const graph::EarDecomposition ears       = schmidt(biconnectedGraph);
  const graph::AdjacencyListGraph fromEars = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());
  return minimallyBiconnectedSubgraph(fromEars);
}

/*!
 * @brief finds a hamilton cycle in the given open ear decomposition
 * @param openEars
 * @param numberOfNodes
 * @return hamilton cycle, first node is not repeated as last
 */
std::vector<size_t> findHamiltonCycleInOpenEarDecomposition(const graph::EarDecomposition& openEars, const size_t numberOfNodes);

/*!
 * @brief approximates a BTSP
 * @tparam G type of graph, must be complete and weighted
 * @param completeGraph complete weighted graph, providing distances between nodes
 * @return Result
 */
template <typename G>
  requires(std::is_base_of_v<graph::CompleteGraph, G> && std::is_base_of_v<graph::WeightedGraph, G>)
Result approximateBTSP(const G& completeGraph) {
  double maxEdgeWeight;
  const graph::AdjacencyListGraph biconnectedGraph = bottleneckOptimalBiconnectedSubgraph(completeGraph, maxEdgeWeight);
  const graph::AdjacencyListGraph minimal          = makeMinimallyBiconnected(biconnectedGraph);
  const graph::EarDecomposition openEars           = schmidt(minimal);  // calculate proper ear decomposition
  const std::vector<size_t> tour                   = findHamiltonCycleInOpenEarDecomposition(openEars, completeGraph.numberOfNodes());
  const graph::Edge bottleneckEdge                 = findBottleneck(completeGraph, tour, true);
  const double objective                           = completeGraph.weight(bottleneckEdge);

  assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
  return Result{biconnectedGraph, openEars, tour, bottleneckEdge, objective, maxEdgeWeight, minimal.numberOfEdges()};
}

/*!
 * @brief extracts hamiton path from hamilton cycle in fivefold graph
 * @param wholeTour hamilton cycle in fivefold graph
 * @param s start node
 * @param t end node
 * @return std::vector<size_t> hamiltonian s-t path in original graph
 */
std::vector<size_t> extractHamiltonPath(const std::vector<size_t>& wholeTour, const size_t s, const size_t t);

/*!
 * @brief removes edges that are not 2-essential if the egde (s,t) is added to the graph
 * @details The graph is minimally biconnected when adding the edge (s,t).
 * @param biconnectedGraph
 * @param s start node
 * @param t end node
 * @return graph that is biconnected when (s,t) is added
 */
template <typename G>
  requires(std::is_base_of_v<graph::Graph, G>)
graph::AdjacencyListGraph makeEdgeAugmentedMinimallyBiconnected(const G& biconnectedGraph, const size_t s, const size_t t) {
  const graph::Edge st_Edge{s, t};
  const graph::EarDecomposition ears = schmidt(biconnectedGraph);
  graph::AdjacencyListGraph fromEars = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());
  if (!fromEars.adjacent(s, t)) {  // if the s-t edge is one of the removed ones,
    fromEars.addEdge(st_Edge);     // add it again.
  }
  graph::AdjacencyListGraph minimal = edgeKeepingMinimallyBiconectedSubgraph(fromEars, st_Edge);
  minimal.removeEdge(st_Edge);

  return minimal;
}

/*!
 * @brief creates graph consisting of 5 copies of the original graph
 * @details Copies original graph, adds 2 nodes x and y and connects x to all copies of s and y to all copies of t
 * @param minimallyBiconnected
 * @param s start node
 * @param t end node
 * @return graph consisting of 5 copies of minimallyBiconnected plus nodes x and y
 */
graph::AdjacencyListGraph createFiveFoldGraph(const graph::AdjacencyListGraph& minimallyBiconnected, const size_t s, const size_t t);

/*!
 * @brief approximates the BTSPP
 * @param euclidean complete graph, providing distances between nodes
 * @param s start node
 * @param t end node
 * @param printInfo controls if objective, lower bound on OPT and a fortiori guarantee are printed to console
 * @return Result
 */
template <typename G>
  requires(std::is_base_of_v<graph::CompleteGraph, G> && std::is_base_of_v<graph::WeightedGraph, G>)
Result approximateBTSPP(const G& completeGraph, const size_t s = 0, const size_t t = 1) {
  double maxEdgeWeight;

  // find graph s.t. G = (V,E) + (s,t) is biconnected
  const graph::AdjacencyListGraph biconnectedGraph = edgeAugmentedBiconnectedSubgraph(completeGraph, graph::Edge{s, t}, maxEdgeWeight);
  const graph::AdjacencyListGraph minimal          = makeEdgeAugmentedMinimallyBiconnected(biconnectedGraph, s, t);
  graph::AdjacencyListGraph fiveFoldGraph          = createFiveFoldGraph(minimal, s, t);
  const size_t numberOfNodes5FoldGraph             = fiveFoldGraph.numberOfNodes();
  const graph::EarDecomposition openEars           = schmidt(fiveFoldGraph);                // calculate open ear decomposition
  std::vector<size_t> wholeTour                    = findHamiltonCycleInOpenEarDecomposition(openEars, numberOfNodes5FoldGraph);
  const std::vector<size_t> tour                   = extractHamiltonPath(wholeTour, s, t);  // extract s-t-path from solution
  const graph::Edge bottleneckEdge                 = findBottleneck(completeGraph, tour, false);
  const double objective                           = completeGraph.weight(bottleneckEdge);

  assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
  return Result{biconnectedGraph, openEars, tour, bottleneckEdge, objective, maxEdgeWeight, minimal.numberOfEdges()};
}
}  // namespace approximation
