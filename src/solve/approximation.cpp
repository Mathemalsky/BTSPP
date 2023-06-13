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
#include "solve/approximation.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "draw/definitions.hpp"

#include "algorithm.hpp"
#include "graph.hpp"
#include "utils.hpp"

#include "solve/commonfunctions.hpp"

namespace approximation {

/***********************************************************************************************************************
 *                                                  output function
 **********************************************************************************************************************/

void printInfo(const approximation::Result& res, const ProblemType problemType, const double runtime) {
  std::cout << "-------------------------------------------------------\n";
  std::cout << "Approximated an instance of " << problemType << "." << std::endl;
  std::cout << "objective                            : " << res.objective << std::endl;
  std::cout << "lower bound on OPT                   : " << res.lowerBoundOnOPT << std::endl;
  std::cout << "a fortiori guarantee                 : " << res.objective / res.lowerBoundOnOPT << std::endl;
  std::cout << "edges in biconnected graph           : " << res.biconnectedGraph.numberOfEdges() << std::endl;
  std::cout << "edges in minimally biconnected graph : " << res.numberOfEdgesInMinimallyBiconectedGraph << std::endl;
  if (runtime != -1.0) {
    std::cout << "elapsed time                         : " << runtime << " ms\n";
  }
}

/***********************************************************************************************************************
 *                                               algorithms for BTSP
 **********************************************************************************************************************/

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
 * @brief doubles and deletes edges to make open ear decomposition eularian
 * @param ears open ear decomposition
 * @param numberofNodes tells the function the ranges of indices ocurring in ears
 */
static graph::AdjacencyListDigraph constructDigraph(const graph::EarDecomposition& ears, const size_t numberOfNodes) {
  graph::AdjacencyListGraph graph = earDecompToAdjacencyListGraph(ears, numberOfNodes);
  graph::AdjacencyListDigraph digraph(numberOfNodes);

  // in the first (or last) ear there is never an edge to be deleted
  const std::vector<size_t>& firstEar = ears.ears.back();
  for (size_t i = 0; i < firstEar.size() - 2; ++i) {
    digraph.addEdge(firstEar[i], firstEar[i + 1]);
  }
  digraph.addEdge(firstEar.back(), firstEar[firstEar.size() - 2]);

  // the other ears deletion of at most one edge can occur
  for (long j = ears.ears.size() - 2; j >= 0; --j) {
    const std::vector<size_t>& ear = ears.ears[static_cast<size_t>(j)];
    const size_t pos_y =
        std::distance(ear.begin(), std::find_if(ear.begin() + 1, ear.end() - 1, [&](size_t u) { return graph.degree(u) == 2; }));

    digraph.addEdge(ear[0], ear[1]);  // add the first edge directing into the ear
    size_t earPosOfLastDoubledEdge = 0;

    // struct to bundle edge and index
    struct EdgeIndex {
      graph::Edge e;
      size_t index;
    };
    std::vector<EdgeIndex> edgesToBeDirected;
    for (size_t i = 1; i < ear.size() - 1; ++i) {
      const size_t u = ear[i];
      const size_t v = ear[i + 1];

      if (graph.degree(u) % 2 == 1) {
        graph.addEdge(u, v);
        digraph.addEdge(u, v);
        digraph.addEdge(v, u);
        earPosOfLastDoubledEdge = i;
      }
      else {
        edgesToBeDirected.emplace_back(graph::Edge{u, v}, i);
      }
    }

    if (earPosOfLastDoubledEdge == ear.size() - 2) {  // if the last edge was doubled
      for (const EdgeIndex& edgeIndex : edgesToBeDirected) {
        digraph.addEdge(edgeIndex.e);
      }

      const graph::Edge lastDoubledEdge{ear[earPosOfLastDoubledEdge], ear[earPosOfLastDoubledEdge + 1]};
      graph.removeEdge(lastDoubledEdge);
      graph.removeEdge(lastDoubledEdge);  // the edge was doubled so it needs to be removed twice

      digraph.removeEdge(lastDoubledEdge);
      digraph.removeEdge(lastDoubledEdge.reverse());
    }
    else {
      const size_t y           = ear[pos_y];
      const size_t y_successor = ear[pos_y + 1];
      if (digraph.adjacent(y, y_successor) && digraph.adjacent(y_successor, y)) {  // edges adjcent to y are doubled
        digraph.removeEdge(y, y_successor);
        digraph.removeEdge(y_successor, y);
      }
      for (const EdgeIndex& edge : edgesToBeDirected) {
        if (edge.index < pos_y) {
          digraph.addEdge(edge.e);
        }
        else {
          digraph.addEdge(edge.e.reverse());
        }
      }
    }
  }
  return digraph;
}

/*!
 * @brief prepairs the graphs for finding the euler tour
 * @param graphPair graph and digraph
 */
static void prepareForEulertour(graph::AdjacencyListDigraph& digraph) {
  const size_t numberOfNodes = digraph.numberOfNodes();

  graph::AdjacencyListDigraph reverseDigraph(numberOfNodes);
  for (const graph::Edge& e : digraph.edges()) {
    reverseDigraph.addEdge(e.reverse());
  }

  digraph.addIsolatedNodes(numberOfNodes);
  for (const size_t u : reverseDigraph.nodes()) {
    if (reverseDigraph.degree(u) == 2 && digraph.degree(u) > 0) {
      const size_t v = reverseDigraph.neighbours(u)[0];
      const size_t w = reverseDigraph.neighbours(u)[1];
      digraph.removeEdge(v, u);
      digraph.removeEdge(w, u);
      digraph.addEdge(v, u + numberOfNodes);
      digraph.addEdge(w, u + numberOfNodes);
    }
  }
}

/*!
 * @brief finds an euler tour in graph
 * @param graph
 * @param digraph holds information needed to construct an euler tour that can be shortcutted to hamiltonian cycle
 * @return std::vector<size_t> of node indices; first is not repeated as last
 */
static std::vector<size_t> findEulertour(graph::AdjacencyListDigraph& digraph) {
  prepareForEulertour(digraph);
  std::vector<size_t> eulertourInGMinus = hierholzer(digraph.undirected());
  return eulertourInGMinus;
}

/*!
 * @brief shortcuts euler tour to hamiltonian cycle
 * @param longEulertour euler tour in multigraph
 * @param digraph holds information for shortcutting
 * @return std::vector<size_t> of node indices, first is not repeated as last
 */
static std::vector<size_t> shortcutToHamiltoncycle(const std::vector<size_t>& longEulertour,
                                                   graph::AdjacencyListDigraph& digraph,
                                                   const size_t numberOfNodes) {
  std::vector<size_t> hamiltoncycle;
  hamiltoncycle.reserve(numberOfNodes);
  hamiltoncycle.push_back(longEulertour[0]);

  for (size_t i = 1; i < longEulertour.size() - 1; ++i) {
    const size_t u = longEulertour[i - 1];
    const size_t v = longEulertour[i];
    const size_t w = longEulertour[i + 1];
    if (digraph.adjacent(v, u) && digraph.adjacent(v, w)) {
      hamiltoncycle.push_back(w);
      ++i;  // do not consider next node for skipping. It's alredy added.
    }
    else {
      hamiltoncycle.push_back(v);
    }
  }

  for (size_t& node : hamiltoncycle) {
    if (node >= numberOfNodes) {
      node -= numberOfNodes;
    }
  }

  return hamiltoncycle;
}

static std::vector<size_t> findHamiltonCycleInOpenEarDecomposition(const graph::EarDecomposition& openEars, const size_t numberOfNodes) {
  std::vector<size_t> tour;
  if (openEars.ears.size() == 1) {
    tour = std::vector<size_t>(openEars.ears[0].begin(), openEars.ears[0].end() - 1);  // do not repeat first node
  }
  else {
    graph::AdjacencyListDigraph digraph = constructDigraph(openEars, numberOfNodes);
    std::vector<size_t> tmp             = findEulertour(digraph);
    tour                                = shortcutToHamiltoncycle(tmp, digraph, numberOfNodes);
  }
  assert(tour.size() == numberOfNodes && "Missmatching number of nodes in hamilton cycle!");
  return tour;
}

Result approximateBTSP(const graph::Euclidean& euclidean) {
  double maxEdgeWeight;
  const graph::AdjacencyListGraph biconnectedGraph = bottleneckOptimalBiconnectedSubgraph(euclidean, maxEdgeWeight);
  const graph::AdjacencyListGraph minimal          = makeMinimallyBiconnected(biconnectedGraph);
  const graph::EarDecomposition openEars           = schmidt(minimal);  // calculate proper ear decomposition
  const std::vector<size_t> tour                   = findHamiltonCycleInOpenEarDecomposition(openEars, euclidean.numberOfNodes());
  const graph::Edge bottleneckEdge                 = findBottleneck(euclidean, tour, true);
  const double objective                           = euclidean.weight(bottleneckEdge);

  assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
  return Result{biconnectedGraph, openEars, tour, bottleneckEdge, objective, maxEdgeWeight, minimal.numberOfEdges()};
}

/***********************************************************************************************************************
 *                                               algorithms for BTSPP
 **********************************************************************************************************************/

static size_t increaseModulo(const size_t number, const size_t modulus) {
  return number == modulus - 1 ? 0 : +1;
}

static size_t decreaseModulo(const size_t number, const size_t modulus) {
  return number == 0 ? modulus - 1 : number - 1;
}

template <typename Type>
static void reverse(std::vector<Type>& vec) {
  const size_t size       = vec.size();
  const size_t sizeMinus1 = size - 1;
  for (size_t i = 0; i < size / 2; ++i) {
    std::swap(vec[i], vec[sizeMinus1 - i]);
  }
}

/*!
 * @brief creates graph consisting of 5 copies of the original graph
 * @details Copies original graph, adds 2 nodes x and y and connects x to all copies of s and y to all copies of t
 */
static graph::AdjacencyListGraph createFiveFoldGraph(const graph::Euclidean& euclidean,
                                                     const graph::AdjacencyListGraph& minimalBiconnected,
                                                     const size_t s,
                                                     const size_t t) {
  const size_t numberOfNodes           = euclidean.numberOfNodes();
  const size_t numberOfNodes5FoldGraph = 5 * numberOfNodes + 2;
  std::vector<std::vector<size_t>> adjacencyList;
  adjacencyList.reserve(numberOfNodes5FoldGraph);

  for (size_t i = 0; i < 5; ++i) {
    adjacencyList.insert(adjacencyList.end(), minimalBiconnected.adjacencyList().begin(), minimalBiconnected.adjacencyList().end());
  }

  for (size_t i = 1; i < 5; ++i) {
    for (size_t j = 0; j < numberOfNodes; ++j) {
      for (size_t& k : adjacencyList[i * numberOfNodes + j]) {
        k += i * numberOfNodes;
      }
    }
  }

  const size_t x = 5 * numberOfNodes;
  const size_t y = 5 * numberOfNodes + 1;
  adjacencyList.resize(numberOfNodes5FoldGraph);  // add empty adjacency vector for x and y
  graph::AdjacencyListGraph fiveFoldGraph(adjacencyList);

  for (size_t i = 0; i < 5; i++) {
    fiveFoldGraph.addEdge(i * numberOfNodes + s, x);  // connect all copies of s to x
    fiveFoldGraph.addEdge(i * numberOfNodes + t, y);  // connect all copies of t to y
  }

  return fiveFoldGraph;
}

/*!
 * @brief extracts hamiton path from hamilton cycle in fivefold graph
 * @param wholeTour hamilton cycle in fivefold graph
 * @param s start node
 * @param t end node
 * @return std::vector<size_t> hamiltonian s-t path in original graph
 */
static std::vector<size_t> extractHamiltonPath(const std::vector<size_t>& wholeTour, const size_t s, const size_t t) {
  const size_t numberOfNodes5FoldGraph = wholeTour.size();
  const size_t numberOfNodes           = (numberOfNodes5FoldGraph - 2) / 5;
  const size_t x                       = numberOfNodes5FoldGraph - 2;
  const size_t y                       = numberOfNodes5FoldGraph - 1;
  const size_t pos_x                   = graph::findPosition(wholeTour, x);
  const size_t pos_y                   = graph::findPosition(wholeTour, y);

  std::array<bool, 5> graphCopyIsSolution{true, true, true, true, true};
  graphCopyIsSolution[wholeTour[increaseModulo(pos_x, numberOfNodes5FoldGraph)] / numberOfNodes] = false;
  graphCopyIsSolution[wholeTour[decreaseModulo(pos_x, numberOfNodes5FoldGraph)] / numberOfNodes] = false;
  graphCopyIsSolution[wholeTour[increaseModulo(pos_y, numberOfNodes5FoldGraph)] / numberOfNodes] = false;
  graphCopyIsSolution[wholeTour[decreaseModulo(pos_y, numberOfNodes5FoldGraph)] / numberOfNodes] = false;

  const size_t solutionIndex = graph::findPosition(graphCopyIsSolution, true);
  const size_t pos_s         = graph::findPosition(wholeTour, solutionIndex * numberOfNodes + s);
  const size_t pos_t         = graph::findPosition(wholeTour, solutionIndex * numberOfNodes + t);
  const size_t maxPos        = std::max(pos_s, pos_t);
  const size_t minPos        = std::min(pos_s, pos_t);

  std::vector<size_t> tour;
  tour.reserve(numberOfNodes);
  if (maxPos - minPos == numberOfNodes - 1) {
    for (size_t i = minPos; i <= maxPos; ++i) {
      tour.push_back(wholeTour[i]);
    }
  }
  else {
    for (size_t i = maxPos; i != minPos; i = increaseModulo(i, numberOfNodes5FoldGraph)) {
      tour.push_back(wholeTour[i]);
    }
    tour.push_back(wholeTour[minPos]);
  }

  for (size_t& node : tour) {
    node -= solutionIndex * numberOfNodes;
  }

  if (tour.front() != s) {
    reverse(tour);
  }
  assert(tour.front() == s && "Hamilton path must start with correct node!");
  assert(tour.back() == t && "Hamilton path must end with correct node!");

  return tour;
}

/*!
 * @brief removes all edges, not 2-essential when adding (s,t)
 * @details The graph is minimally biconnected when adding the edge (s,t).
 * @param biconnectedGraph
 * @param s start node
 * @param t end node
 * @return AdjacencyListGraph minimal with desired property
 */
static graph::AdjacencyListGraph makeEdgeAugmentedMinimallyBiconnected(const graph::AdjacencyListGraph& biconnectedGraph,
                                                                       const size_t s,
                                                                       const size_t t) {
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

Result approximateBTSPP(const graph::Euclidean& euclidean, const size_t s, const size_t t) {
  double maxEdgeWeight;

  // find graph s.t. G = (V,E) + (s,t) is biconnected
  const graph::AdjacencyListGraph biconnectedGraph = edgeAugmentedBiconnectedSubgraph(euclidean, graph::Edge{s, t}, maxEdgeWeight);
  const graph::AdjacencyListGraph minimal          = makeEdgeAugmentedMinimallyBiconnected(biconnectedGraph, s, t);
  graph::AdjacencyListGraph fiveFoldGraph          = createFiveFoldGraph(euclidean, minimal, s, t);
  const size_t numberOfNodes5FoldGraph             = fiveFoldGraph.numberOfNodes();
  const graph::EarDecomposition openEars           = schmidt(fiveFoldGraph);                // calculate open ear decomposition
  std::vector<size_t> wholeTour                    = findHamiltonCycleInOpenEarDecomposition(openEars, numberOfNodes5FoldGraph);
  const std::vector<size_t> tour                   = extractHamiltonPath(wholeTour, s, t);  // extract s-t-path from solution
  const graph::Edge bottleneckEdge                 = findBottleneck(euclidean, tour, false);
  const double objective                           = euclidean.weight(bottleneckEdge);

  assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
  return Result{biconnectedGraph, openEars, tour, bottleneckEdge, objective, maxEdgeWeight, minimal.numberOfEdges()};
}
}  // namespace approximation
