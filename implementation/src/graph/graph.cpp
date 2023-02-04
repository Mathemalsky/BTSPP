#include "graph/graph.hpp"

#include <cassert>
#include <exception>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connected() const {
  std::vector<bool> visited(numberOfNodes(), false);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      for (Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>::InnerIterator it(pAdjacencyMatrix, top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());
        }
      }
      for (unsigned int i = top + 1; i < numberOfNodes(); ++i) {
        if (pAdjacencyMatrix.coeff(i, top) != 0 && !visited[i]) {
          nodeStack.push(i);
        }
      }
      visited[top] = true;
    }
  }

  for (const bool& test : visited) {
    if (!test) {
      return false;
    }
  }
  return true;
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connectedWhithout(const size_t vertex) const {
  std::vector<bool> visited(numberOfNodes(), false);
  std::stack<size_t> nodeStack;
  const size_t rootNode = (vertex != 0 ? 0 : 1);
  nodeStack.push(rootNode);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      for (Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>::InnerIterator it(pAdjacencyMatrix, top); it; ++it) {
        if (!visited[it.index()] && (size_t) it.index() != vertex) {
          nodeStack.push(it.index());
        }
      }
      for (unsigned int i = top + 1; i < numberOfNodes(); ++i) {
        if (pAdjacencyMatrix.coeff(i, top) != 0 && !visited[i] && i != vertex) {
          nodeStack.push(i);
        }
      }
      visited[top] = true;
    }
  }

  for (size_t i = 0; i < visited.size(); ++i) {
    if (!visited[i] && vertex != i) {
      return false;
    }
  }
  return true;
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::biconnected() const {
  for (size_t i = 0; i < numberOfNodes(); ++i) {
    if (!connectedWhithout(i)) {
      return false;
    }
  }
  return true;
}

static AdjacencyListGraph<Directionality::Undirected> findBackedges(
    const AdjacencyMatrixGraph<Directionality::Undirected>& graph, const DfsTree& tree) {
  AdjacencyListGraph<Directionality::Undirected> backedges(graph.numberOfNodes());
  for (Edge e : graph) {
    if (!tree.adjacent(e) && !tree.adjacent(e.reverse())) {
      backedges.addEdge(e);
    }
  }
  return backedges;
}

OpenEarDecomposition schmidt(const AdjacencyMatrixGraph<Directionality::Undirected>& graph) {
  const DfsTree tree                                             = dfs(graph);
  const AdjacencyListGraph<Directionality::Undirected> backedges = findBackedges(graph, tree);
  const size_t numberOfNodes                                     = graph.numberOfNodes();

  std::vector<bool> visited(numberOfNodes, false);
  std::vector<std::vector<size_t>> ears;
  for (size_t v : tree.exploratioOrder()) {     // iterate over all nodes in the order they appeared in dfs
    for (size_t u : backedges.neighbours(v)) {  // for every backedge starting at v
      if (!visited[u]) {
        std::vector<size_t> chain{v, u};
        visited[v] = true;
        while (!visited[u]) {
          visited[u] = true;
          u          = tree.parent(u);
          chain.push_back(u);
        }
        assert((u != v || v == 0) && "Graph is not biconnected!");
        ears.push_back(chain);
      }
    }
  }
  return OpenEarDecomposition{ears};
}

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

AdjacencyMatrixGraph<Directionality::Undirected> biconnectedSpanningGraph(const Euclidean& euclidean) {
  const size_t numberOfNodes = euclidean.numberOfNodes();

  const Index index(numberOfNodes);
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

  // create an undirected graph from that
  AdjacencyMatrixGraph<Directionality::Undirected> graph(numberOfNodes, entries);

  // continue adding edges until it is biconnected
  unsigned int edgeCounter = numberOfNodes;

  while (!graph.biconnected()) {
    assert(edgeCounter < euclidean.numberOfEdges() && "We cannot add more edges than existing.");

    const Edge e = index.edge(edgeIndeces[edgeCounter++]);
    graph.addEdge(e, euclidean.weight(e));
  }

  graph.compressMatrix();  // matrix became uncommpressed when adding edges
  return graph;
}
