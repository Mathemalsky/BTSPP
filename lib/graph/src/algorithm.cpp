/*
 * GRAPH is a library to store and manipulate graphs as adjacency list or
 * as sparse eigen matrix. Different specialized types of graphs are
 * supported.
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
#include "algorithm.hpp"

#include <cmath>
#include <cstddef>

#include <Eigen/SparseCore>

namespace graph {

AdjacencyListGraph earDecompToAdjacencyListGraph(const EarDecomposition& earDecomposition, const size_t numberOfNodes) {
  std::vector<std::vector<size_t>> adjacencyList(numberOfNodes);
  for (const std::vector<size_t>& ear : earDecomposition.ears) {
    for (size_t i = 0; i < ear.size() - 1; ++i) {
      adjacencyList[ear[i]].push_back(ear[i + 1]);
      adjacencyList[ear[i + 1]].push_back(ear[i]);
    }
  }
  return AdjacencyListGraph(adjacencyList);
}

AdjacencyListGraph minimallyBiconnectedSubgraph(const AdjacencyListGraph& graph) {
  AdjacencyListGraph saveCopy         = graph;
  AdjacencyListGraph experimentalCopy = graph;

  for (const Edge& e : graph.edges()) {
    if (experimentalCopy.degree(e.u) > 2 && experimentalCopy.degree(e.v) > 2) {
      experimentalCopy.removeEdge(e);
      if (experimentalCopy.biconnected()) {
        saveCopy = experimentalCopy;
      }
      else {
        experimentalCopy = saveCopy;
      }
    }
  }
  return saveCopy;
}

AdjacencyListGraph edgeKeepingMinimallyBiconectedSubgraph(const AdjacencyListGraph& graph, const Edge& keepEdge) {
  AdjacencyListGraph saveCopy         = graph;
  AdjacencyListGraph experimentalCopy = graph;

  for (const Edge& e : graph.edges()) {
    if (experimentalCopy.degree(e.u) > 2 && experimentalCopy.degree(e.v) > 2 && e != keepEdge && e.reverse() != keepEdge) {
      experimentalCopy.removeEdge(e);
      if (experimentalCopy.biconnected()) {
        saveCopy = experimentalCopy;
      }
      else {
        experimentalCopy = saveCopy;
      }
    }
  }
  return saveCopy;
}

}  // namespace graph
