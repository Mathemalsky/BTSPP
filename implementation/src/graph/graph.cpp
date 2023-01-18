#include "graph/graph.hpp"

#include <stack>
#include <vector>

#include <Eigen/SparseCore>

bool AdjacencyMatrixGraph::connected() const {
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

bool AdjacencyMatrixGraph::connectedWhithout(const size_t vertex) const {
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

bool AdjacencyMatrixGraph::biconnected() const {
  for (size_t i = 0; i < pNumberOfNodes; ++i) {
    if (!connectedWhithout(i)) {
      return false;
    }
  }
  return true;
}

OpenEarDecomposition schmidt(const UndirectedGraph<AdjacencyMatrixGraph>& graph) {
  // dfs
  const size_t numberOfNodes = graph.numberOfNodes();
  DfsTree tree(numberOfNodes);
  // graph with adjacency list
  size_t indexCounter = 0;
  std::vector<bool> visited(numberOfNodes, false);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);  // w. l. o. g. the root node is 0
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      visited[top] = true;
      for (Eigen::SparseMatrix<EdgeCost>::InnerIterator it(graph.matrix(), top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());     // do not push already visited nodes
          tree.parent(it.index()) = top;  // update the parent
        }
      }
      tree.index(top) = indexCounter;
      ++indexCounter;
    }
  }

  // decompose by iterating over every backedge outgoing from every node
  // for()
  std::vector<std::vector<size_t>> ears;
  for (size_t v = 0; v < numberOfNodes; ++v) {
  }
}
