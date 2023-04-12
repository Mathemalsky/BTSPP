#include "graph/algorithm.hpp"

#include <cmath>
#include <cstddef>

#include <Eigen/SparseCore>

// DEBUG
#include <iostream>
#include "graph/ostream.hpp"

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

  size_t edgeIndex(const size_t i, const size_t j) const { return i > j ? 0.5 * (i * i - i) + j : 0.5 * (j * j - j) + i; }
  size_t edgeIndex(const Edge& edge) const { return edgeIndex(edge.u, edge.v); }

  Edge edge(const size_t k) const {
    const size_t i = std::floor(std::sqrt(0.25 + 2 * k) + 0.5);
    return Edge{i, k - (i * i - i) / 2};
  }

private:
  const size_t pNumberOfNodes;
};

static std::vector<size_t> createEdgeIndeces(const Euclidean& euclidean) {
  std::vector<size_t> edgeIndices(euclidean.numberOfEdges());
  std::iota(edgeIndices.begin(), edgeIndices.end(), 0);
  return edgeIndices;
}

static AdjacencyMatrixGraph addEdgesUntilBiconnected(const Euclidean& euclidean, const Index index, const std::vector<size_t>& edgeIndices,
                                                     double& maxEdgeWeight) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  const size_t numberOfEdges = euclidean.numberOfEdges();

  // add the first numberOfNodes many edges
  std::vector<Entry> entries;
  entries.reserve(numberOfNodes);
  for (size_t i = 0; i < numberOfNodes; ++i) {
    const Edge e = index.edge(edgeIndices[i]);
    entries.push_back(Entry(e.u, e.v, euclidean.weight(e)));
    entries.push_back(Entry(e.v, e.u, euclidean.weight(e)));
  }
  AdjacencyMatrixGraph graph(numberOfNodes, entries);  // create an undirected graph from that
  AdjacencyMatrixGraph graphCopy = graph;              // and a copy

  // use bisection search to find bottleneck optimal biconnected subgraph
  size_t upperbound = numberOfEdges;
  size_t lowerbound = numberOfNodes;
  while (upperbound != lowerbound) {
    size_t middle = (lowerbound + upperbound) / 2;
    for (size_t i = lowerbound; i < middle; ++i) {
      const Edge e = index.edge(edgeIndices[i]);
      graphCopy.addEdge(e, euclidean.weight(e));
    }
    graphCopy.compressMatrix();

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

  // DEBUG
  std::cerr << "edge for lower bound: " << index.edge(edgeIndices[lowerbound - 1]) << std::endl;

  maxEdgeWeight = euclidean.weight(index.edge(edgeIndices[lowerbound - 1]));  // for lower bound on opt
  return graph;
}

AdjacencyMatrixGraph biconnectedSubgraph(const Euclidean& euclidean, double& maxEdgeWeight) {
  // sort the edges
  const Index index(euclidean.numberOfNodes());
  std::vector<size_t> edgeIndices = createEdgeIndeces(euclidean);
  std::sort(edgeIndices.begin(), edgeIndices.end(), [euclidean, index](const size_t a, const size_t b) {
    return euclidean.weight(index.edge(a)) < euclidean.weight(index.edge(b));
  });

  AdjacencyMatrixGraph graph = addEdgesUntilBiconnected(euclidean, index, edgeIndices, maxEdgeWeight);
  graph.compressMatrix();  // matrix became uncommpressed when adding edges
  return graph;
}

AdjacencyMatrixGraph edgeAugmentedBiconnectedSubgraph(const Euclidean& euclidean, const Edge augmentationEdge, double& maxEdgeWeight) {
  // sort the edges, put the augmentation edge in first position
  const Index index(euclidean.numberOfNodes());
  std::vector<size_t> edgeIndices    = createEdgeIndeces(euclidean);
  const size_t augmentationEdgeIndex = findPosition(edgeIndices, index.edgeIndex(augmentationEdge));
  std::swap(edgeIndices.front(), edgeIndices[augmentationEdgeIndex]);
  std::sort(edgeIndices.begin() + 1, edgeIndices.end(), [euclidean, index, augmentationEdge](const size_t a, const size_t b) {
    return euclidean.weight(index.edge(a)) < euclidean.weight(index.edge(b));
  });

  // DEBUG
  std::cerr << "edgPos of augmentation edge: "
            << std::distance(edgeIndices.begin(),
                             std::find(edgeIndices.begin(), edgeIndices.end(), index.edgeIndex(augmentationEdge.u, augmentationEdge.v)))
            << std::endl;

  AdjacencyMatrixGraph graph = addEdgesUntilBiconnected(euclidean, index, edgeIndices, maxEdgeWeight);
  // graph.removeEdge(augmentationEdge);
  // graph.pruneMatrix();
  graph.compressMatrix();  // matrix became uncommpressed when adding edges
  return graph;
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
