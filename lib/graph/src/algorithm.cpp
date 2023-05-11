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

/*!
 * @brief precompute the values to improve performance
 */
class ExplicitEdges {
public:
  struct EdgeInfo {
    Edge edge;
    double weightSquared;
  };

  ExplicitEdges(const ExplicitEdges&) = delete;  // copying would destroy the pointers

  ExplicitEdges(const Euclidean& euclidean) {
    pEdges.reserve(euclidean.numberOfEdges());
    pEdgePointer.reserve(euclidean.numberOfEdges());
    for (const Edge& edge : euclidean.edges()) {
      pEdges.emplace_back(edge, euclidean.weightSquared(edge));
      pEdgePointer.push_back(&pEdges.back());
    }
  }

  ExplicitEdges(const Euclidean& euclidean, const Edge& augmentationEdge) {
    pEdges.reserve(euclidean.numberOfEdges());
    pEdgePointer.reserve(euclidean.numberOfEdges());
    pEdges.emplace_back(augmentationEdge, euclidean.weightSquared(augmentationEdge));
    pEdgePointer.push_back(&pEdges.back());

    for (const Edge& edge : euclidean.edges()) {
      if (edge != augmentationEdge) {
        pEdges.emplace_back(edge, euclidean.weightSquared(edge));
        pEdgePointer.push_back(&pEdges.back());
      }
    }
  }

  std::vector<const EdgeInfo*>& edgePointers() { return pEdgePointer; }

private:
  std::vector<EdgeInfo> pEdges;
  std::vector<const EdgeInfo*> pEdgePointer;
};

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

static AdjacencyListGraph addEdgesUntilBiconnected(const size_t numberOfNodes,
                                                   double& maxEdgeWeight,
                                                   const std::vector<const ExplicitEdges::EdgeInfo*>& edges) {
  const size_t numberOfEdges = edges.size();

  // add the first numberOfNodes many edges
  AdjacencyListGraph graph(numberOfNodes);
  for (size_t i = 0; i < numberOfNodes; ++i) {
    graph.addEdge(edges[i]->edge);
  }
  AdjacencyListGraph graphCopy = graph;  // and a copy

  // use bisection search to find bottleneck optimal biconnected subgraph
  size_t upperbound = numberOfEdges;
  size_t lowerbound = numberOfNodes;
  while (upperbound != lowerbound) {
    size_t middle = (lowerbound + upperbound) / 2;
    for (size_t i = lowerbound; i < middle; ++i) {
      graphCopy.addEdge(edges[i]->edge);
    }

    if (graphCopy.biconnected()) {
      upperbound = middle;
      graphCopy  = graph;
    }
    else {
      lowerbound = middle + 1;
      graphCopy.addEdge(edges[middle]->edge);
      graph = graphCopy;
    }
  }
  maxEdgeWeight = std::sqrt(edges[lowerbound - 1]->weightSquared);  // for lower bound on opt
  return graph;
}

AdjacencyListGraph biconnectedSubgraph(const Euclidean& euclidean, double& maxEdgeWeight) {
  // precompute the squared edge weights
  ExplicitEdges explicitEdges(euclidean);
  std::vector<const ExplicitEdges::EdgeInfo*> edges = explicitEdges.edgePointers();

  // sort the edges
  std::sort(edges.begin(), edges.end(), [&](const ExplicitEdges::EdgeInfo* a, const ExplicitEdges::EdgeInfo* b) {
    return a->weightSquared < b->weightSquared;
  });

  return addEdgesUntilBiconnected(euclidean.numberOfNodes(), maxEdgeWeight, edges);
}

AdjacencyListGraph edgeAugmentedBiconnectedSubgraph(const Euclidean& euclidean, Edge augmentationEdge, double& maxEdgeWeight) {
  assert(augmentationEdge.u != augmentationEdge.v && "Start node and end node must be different!");
  if (augmentationEdge.u < augmentationEdge.v) {
    augmentationEdge.invert();
  }

  ExplicitEdges explicitEdges(euclidean, augmentationEdge);
  std::vector<const ExplicitEdges::EdgeInfo*> edges = explicitEdges.edgePointers();
  std::sort(edges.begin() + 1, edges.end(), [&](const ExplicitEdges::EdgeInfo* a, const ExplicitEdges::EdgeInfo* b) {
    return a->weightSquared < b->weightSquared;
  });

  return addEdgesUntilBiconnected(euclidean.numberOfNodes(), maxEdgeWeight, edges);
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
