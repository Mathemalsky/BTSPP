#include "approximation.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"

// DEBUG
#include <iostream>

namespace approximation {

using Entry = Eigen::Triplet<EdgeWeight>;

class Index {
public:
  Index(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

  size_t numberOfEdges() const { return pNumberOfNodes * (pNumberOfNodes - 1) / 2; }
  size_t edgeIndex(const size_t i, const size_t j) const {
    return i > j ? 0.5 * (i * i - i) + j : 0.5 * (j * j - j) + i;
  }

  Edge edge(const unsigned int k) const {
    const unsigned int i = std::floor(std::sqrt(0.25 + 2 * k) + 0.5);
    return Edge{i, k - (i * i - i) / 2};
  }

private:
  const size_t pNumberOfNodes;
};

Result approximate(const Euclidean& euclidean, const ProblemType problemType) {
  if (problemType == ProblemType::BTSP_approx) {
    const Index index(euclidean.numberOfNodes());
    const size_t numberOfNodes = euclidean.numberOfNodes();
    std::vector<unsigned int> edgeIndeces(index.numberOfEdges());
    std::iota(edgeIndeces.begin(), edgeIndeces.end(), 0);
    std::sort(edgeIndeces.begin(), edgeIndeces.end(), [euclidean, index](const unsigned int a, const unsigned int b) {
      return euclidean.weight(index.edge(a)) < euclidean.weight(index.edge(b));
    });

    // add the first numberOfNodes many edges
    std::vector<Entry> entries;
    entries.reserve(numberOfNodes);
    for (unsigned int i = 0; i < numberOfNodes; ++i) {
      const Edge e = index.edge(edgeIndeces[i]);

      entries.push_back(Entry(e.u, e.v, euclidean.weight(e)));
    }

    // DEBUG
    std::cerr << "#entries in approximate: " << entries.size() << std::endl;

    // DEBUG
    std::cerr << "#nodes in old graph: " << euclidean.numberOfNodes() << std::endl;

    // create an undirected graph from that
    AdjacencyMatrixGraph<Directionality::Undirected> graph(numberOfNodes, entries);

    // continue adding edges until it is biconnected
    unsigned int edgeCounter = numberOfNodes;

    // DEBUG
    std::cerr << "#nodes in new graph: " << graph.numberOfNodes() << std::endl;

    while (!graph.biconnected()) {
      assert(edgeCounter < euclidean.numberOfEdges() && "We cannot add more edges than existing.");

      const Edge e = index.edge(edgeIndeces[edgeCounter]);
      ++edgeCounter;

      graph.addEdge(e, euclidean.weight(e));
    }

    // take the square of the graph
    // graph.square();

    // calculate proper ear decomposition
    const OpenEarDecomposition ears = schmidt(graph);
    // find approximation

    return Result{graph, ears};
  }
  else {
    throw std::runtime_error("Unknown problem type");
  }
}
}  // namespace approximation
