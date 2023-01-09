#include "approximation.hpp"

#include <cmath>
#include <numeric>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"

using Entry = Eigen::Triplet<EdgeCost>;

class Index {
public:
  Index(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

  size_t numberOfEdges() const { return pNumberOfNodes * (pNumberOfNodes - 1) / 2; }
  size_t edgeIndex(const size_t i, const size_t j) const {
    return i > j ? 0.5 * (i * i + i) + j : 0.5 * (j * j + j) + i;
  }

  Edge edge(const unsigned int k) const {
    const unsigned int i = std::floor(std::sqrt(0.25 + 2 * k) - 0.5);
    return Edge{i, k - (i * i + i) / 2};
  }

private:
  const size_t pNumberOfNodes;
};

std::vector<unsigned int> approximate(const Euclidean& euclidean, const ProblemType problemType) {
  const Index index(euclidean.numberOfNodes());
  const size_t numberOfNodes = euclidean.numberOfNodes();
  std::vector<unsigned int> edgeIndeces(index.numberOfEdges());
  std::iota(edgeIndeces.begin(), edgeIndeces.end(), 0);
  std::sort(edgeIndeces.begin(), edgeIndeces.end(), [euclidean, index](const unsigned int a, const unsigned int b) {
    return euclidean.distance(index.edge(a)) < euclidean.distance(index.edge(b));
  });

  // add the first numberOfNodes many edges
  std::vector<Entry> entries;
  entries.reserve(numberOfNodes);
  for (unsigned int i = 0; i < numberOfNodes; ++i) {
    const Edge e = index.edge(edgeIndeces[i]);
    entries.push_back(Entry(e.u, e.v, euclidean.distance(e)));
  }

  Eigen::SparseMatrix<EdgeCost> adjacencyMatrix(numberOfNodes, numberOfNodes);
  adjacencyMatrix.setFromTriplets(entries.begin(), entries.end());

  // check the graph is biconnected
  // calculate proper ear decomposition
  // find approximation
}
