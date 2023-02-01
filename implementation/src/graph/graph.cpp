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
      for (unsigned int k = 0; k < top; ++k) {
        if (pAdjacencyMatrix.coeff(top, k) != 0.0) {
          nodeStack.push(k);
        }
      }
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());
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
      for (unsigned int k = 0; k < top; ++k) {
        if (pAdjacencyMatrix.coeff(top, k) != 0.0 && !visited[k] && k != vertex) {
          nodeStack.push(k);
        }
      }
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, top); it; ++it) {
        if (!visited[it.index()] && (size_t) it.index() != vertex) {
          nodeStack.push(it.index());
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

AdjacencyListGraph<Directionality::Undirected> findBackedges(
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
        assert(u != v || v == 0 && "Graph is not biconnected!");
        ears.push_back(chain);
      }
    }
  }
  return OpenEarDecomposition{ears};
}
