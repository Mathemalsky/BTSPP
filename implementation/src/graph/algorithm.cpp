#include "graph/algorithm.hpp"

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

  size_t numberOfEdges() const { return pNumberOfNodes * (pNumberOfNodes - 1) / 2; }
  size_t edgeIndex(const size_t i, const size_t j) const { return i > j ? 0.5 * (i * i - i) + j : 0.5 * (j * j - j) + i; }

  Edge edge(const unsigned int k) const {
    const unsigned int i = std::floor(std::sqrt(0.25 + 2 * k) + 0.5);
    return Edge{i, k - (i * i - i) / 2};
  }

private:
  const size_t pNumberOfNodes;
};

AdjacencyMatrixGraph biconnectedSpanningGraph(const Euclidean& euclidean, double& maxEdgeWeight) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  const size_t numberOfEdges = euclidean.numberOfEdges();

  // sort the edges
  const Index index(numberOfNodes);
  std::vector<unsigned int> edgeIndices(index.numberOfEdges());
  std::iota(edgeIndices.begin(), edgeIndices.end(), 0);
  std::sort(edgeIndices.begin(), edgeIndices.end(), [euclidean, index](const unsigned int a, const unsigned int b) {
    return euclidean.weight(index.edge(a)) < euclidean.weight(index.edge(b));
  });

  // add the first numberOfNodes many edges
  std::vector<Entry> entries;
  entries.reserve(numberOfNodes);
  for (unsigned int i = 0; i < numberOfNodes; ++i) {
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
  maxEdgeWeight = euclidean.weight(index.edge(edgeIndices[lowerbound - 1]));  // for lower bound on opt
  graph.compressMatrix();                                                     // matrix became uncommpressed when adding edges
  return graph;
}
}  // namespace graph
