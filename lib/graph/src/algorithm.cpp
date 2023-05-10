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
using Entry = Eigen::Triplet<EdgeWeight>;

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

class Index {
public:
  Index(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

  size_t edgeIndex(const Edge& edge) const { return edge.u * pNumberOfNodes + edge.v; }

  Edge edge(const size_t index) const {
    const size_t u = index / pNumberOfNodes;
    const size_t v = index - u * pNumberOfNodes;  // avoiding modulo operation appeared to be faster
    return Edge{u, v};
  }

private:
  const size_t pNumberOfNodes;
};

static std::vector<size_t> createEdgeIndeces(const Euclidean& euclidean, const Index& index) {
  std::vector<size_t> edgeIndices;
  edgeIndices.reserve(euclidean.numberOfEdges());
  for (const Edge& edge : euclidean.edges()) {
    edgeIndices.push_back(index.edgeIndex(edge));
  }
  return edgeIndices;
}

static AdjacencyListGraph addEdgesUntilBiconnected(const Euclidean& euclidean,
                                                   const Index index,
                                                   const std::vector<size_t>& edgeIndices,
                                                   double& maxEdgeWeight) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  const size_t numberOfEdges = euclidean.numberOfEdges();

  // add the first numberOfNodes many edges
  AdjacencyListGraph graph(numberOfNodes);
  for (size_t i = 0; i < numberOfNodes; ++i) {
    const Edge e = index.edge(edgeIndices[i]);
    graph.addEdge(e);
  }
  AdjacencyListGraph graphCopy = graph;  // and a copy

  // use bisection search to find bottleneck optimal biconnected subgraph
  size_t upperbound = numberOfEdges;
  size_t lowerbound = numberOfNodes;
  while (upperbound != lowerbound) {
    size_t middle = (lowerbound + upperbound) / 2;
    for (size_t i = lowerbound; i < middle; ++i) {
      const Edge e = index.edge(edgeIndices[i]);
      graphCopy.addEdge(e, euclidean.weight(e));
    }

    if (graphCopy.biconnected()) {
      upperbound = middle;
      graphCopy  = graph;
    }
    else {
      lowerbound   = middle + 1;
      const Edge e = index.edge(edgeIndices[middle]);
      graphCopy.addEdge(e, euclidean.weight(e));
      graph = graphCopy;
    }
  }
  maxEdgeWeight = euclidean.weight(index.edge(edgeIndices[lowerbound - 1]));  // for lower bound on opt
  return graph;
}

AdjacencyListGraph biconnectedSubgraph(const Euclidean& euclidean, double& maxEdgeWeight) {
  // sort the edges
  const Index index(euclidean.numberOfNodes());
  std::vector<size_t> edgeIndices = createEdgeIndeces(euclidean, index);
  std::sort(edgeIndices.begin(), edgeIndices.end(), [&euclidean, &index](const size_t a, const size_t b) {
    return euclidean.weightSquared(index.edge(a)) < euclidean.weightSquared(index.edge(b));
  });

  return addEdgesUntilBiconnected(euclidean, index, edgeIndices, maxEdgeWeight);
}

AdjacencyListGraph edgeAugmentedBiconnectedSubgraph(const Euclidean& euclidean, Edge augmentationEdge, double& maxEdgeWeight) {
  assert(augmentationEdge.u != augmentationEdge.v && "Start node and end node must be different!");
  if (augmentationEdge.u < augmentationEdge.v) {
    augmentationEdge.invert();
  }
  // sort the edges, put the augmentation edge in first position
  const Index index(euclidean.numberOfNodes());
  std::vector<size_t> edgeIndices            = createEdgeIndeces(euclidean, index);
  const size_t augmentationEdgeIndexPosition = findPosition(edgeIndices, index.edgeIndex(augmentationEdge));

  assert(augmentationEdgeIndexPosition < edgeIndices.size() && "Augmentation edge must exist in euclidean graph!");

  std::swap(edgeIndices.front(), edgeIndices[augmentationEdgeIndexPosition]);
  std::sort(edgeIndices.begin() + 1, edgeIndices.end(), [&euclidean, &index, &augmentationEdge](const size_t a, const size_t b) {
    return euclidean.weight(index.edge(a)) < euclidean.weight(index.edge(b));
  });

  return addEdgesUntilBiconnected(euclidean, index, edgeIndices, maxEdgeWeight);
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
