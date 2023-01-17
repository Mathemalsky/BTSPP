#include "graph/graph.hpp"

#include <stack>
#include <vector>

#include <Eigen/SparseCore>

bool UndirectedGraph::connected() const {
  std::vector<bool> component(pNumberOfNodes, false);
  component[0] = true;  // without loss of generality we consider the connected component containing node 0
  for (int k = 0; k < pAdjacencyMatrix.outerSize(); ++k) {
    if (!component[k]) {
      return false;  // node k is not connected to any node < k, because adjacency matrix is lower triangular
    }
    for (Eigen::SparseMatrix<EdgeCost>::InnerIterator it(pAdjacencyMatrix, k); it; ++it) {
      component[it.index()] = true;
    }
  }
  return true;
}

bool UndirectedGraph::connectedWhithout(const size_t vertex) const {
  std::vector<bool> component(pNumberOfNodes, false);
  component[0] = true;         // without loss of generality we consider the connected component containing node 0
  component[1] = vertex == 0;  // if we check connectivity whithout 0, we start by one

  for (unsigned int k = 0; k < pAdjacencyMatrix.outerSize(); ++k) {
    if (k != vertex) {
      if (!component[k]) {
        return false;  // node k is not connected to any node < k, because adjacency matrix is lower triangular
      }
      for (Eigen::SparseMatrix<EdgeCost>::InnerIterator it(pAdjacencyMatrix, k); it; ++it) {
        component[it.index()] = true;
      }
    }
  }
  return true;
}

bool UndirectedGraph::biconnected() const {
  for (size_t i = 0; i < pNumberOfNodes; ++i) {
    if (!connectedWhithout(i)) {
      return false;
    }
  }
  return true;
}

OpenEarDecomposition schmidt(const UndirectedGraph& graph) {
  // dfs
  std::vector<bool> visited(graph.numberOfNodes(), false);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);  // w. l. o. g. the root node is 0
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      visited[top] = true;
      for (Eigen::SparseMatrix<EdgeCost>::InnerIterator it(graph.matrix(), top); it; ++it) {
        nodeStack.push(it.index());
      }
    }
  }
}
