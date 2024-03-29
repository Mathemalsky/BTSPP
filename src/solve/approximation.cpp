/*
 * BTSPP is a tool to solve, approximate and draw instances of BTSVPP,
 * BTSPP, BTSP and TSP. Drawing is limited to euclidean graphs.
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
 * @brief doubles and deletes edges to make open ear decomposition eularian
 * @param ears open ear decomposition
 * @param numberOfNodes tells the function the ranges of indices ocurring in ears
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
 * @brief prepairs the graphs for finding the Euler tour
 * @param digraph directed graph where all graphs have even degree
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
 * @brief finds an Euler tour in graph
 * @param digraph holds information needed to construct an Euler tour that can be shortcutted to hamiltonian cycle
 * @return std::vector<size_t> of node indices; first is not repeated as last
 */
static std::vector<size_t> findEulertour(graph::AdjacencyListDigraph& digraph) {
  prepareForEulertour(digraph);
  std::vector<size_t> eulertourInGMinus = hierholzer(digraph.undirected());
  return eulertourInGMinus;
}

/*!
 * @brief shortcuts Euler tour to hamiltonian cycle
 * @param longEulertour Euler tour in multigraph
 * @param digraph holds information for shortcutting
 * @param numberOfNodes number of nodes in the oroginal graph
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

std::vector<size_t> findHamiltonCycleInOpenEarDecomposition(const graph::EarDecomposition& openEars, const size_t numberOfNodes) {
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

/*
 * creates graph consisting of 5 copies of the original graph
 * copies original graph, adds 2 nodes x and y and connects x to all copies of s and y to all copies of t
 */
graph::AdjacencyListGraph createFiveFoldGraph(const graph::AdjacencyListGraph& minimallyBiconnected, const size_t s, const size_t t) {
  const size_t numberOfNodes           = minimallyBiconnected.numberOfNodes();
  const size_t numberOfNodes5FoldGraph = 5 * numberOfNodes + 2;
  std::vector<std::vector<size_t>> adjacencyList;
  adjacencyList.reserve(numberOfNodes5FoldGraph);

  for (size_t i = 0; i < 5; ++i) {
    adjacencyList.insert(adjacencyList.end(), minimallyBiconnected.adjacencyList().begin(), minimallyBiconnected.adjacencyList().end());
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

std::vector<size_t> extractHamiltonPath(const std::vector<size_t>& wholeTour, const size_t s, const size_t t) {
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
}  // namespace approximation
