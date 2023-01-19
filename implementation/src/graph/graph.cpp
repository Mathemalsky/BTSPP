#include "graph/graph.hpp"

#include <stack>
#include <vector>

#include <Eigen/SparseCore>

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connected() const {
  std::vector<bool> component(pNumberOfNodes, false);
  component[0] = true;  // without loss of generality we consider the connected component containing node 0
  for (int k = 0; k < pAdjacencyMatrix.outerSize(); ++k) {
    if (!component[k]) {
      return false;  // node k is not connected to any node < k, because adjacency matrix is lower triangular
    }
    for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, k); it; ++it) {
      component[it.index()] = true;
    }
  }
  return true;
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connectedWhithout(const size_t vertex) const {
  std::vector<bool> component(pNumberOfNodes, false);
  component[0] = true;         // without loss of generality we consider the connected component containing node 0
  component[1] = vertex == 0;  // if we check connectivity whithout 0, we start by one

  for (unsigned int k = 0; k < pAdjacencyMatrix.outerSize(); ++k) {
    if (k != vertex) {
      if (!component[k]) {
        return false;  // node k is not connected to any node < k, because adjacency matrix is lower triangular
      }
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, k); it; ++it) {
        component[it.index()] = true;
      }
    }
  }
  return true;
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::biconnected() const {
  for (size_t i = 0; i < pNumberOfNodes; ++i) {
    if (!connectedWhithout(i)) {
      return false;
    }
  }
  return true;
}

template <>
bool AdjacencyListGraph<Directionality::Directed>::connected() const {
  std::vector<bool> component(pNumberOfNodes, false);
  std::vector<size_t> nodesToCheck;
  nodesToCheck.reserve(pNumberOfNodes);
  nodesToCheck.push_back(0);  // we check if the component containing vertex 0 is the whole graph
  while (!nodesToCheck.empty()) {
    const size_t v = nodesToCheck.back();
    nodesToCheck.pop_back();
    component[v] = true;
    for (const size_t& neighbour : pAdjacencyList[v]) {
      if (!component[neighbour]) {
        nodesToCheck.push_back(neighbour);
      }
    }
  }
  for (const bool& node : component) {
    if (!node) {
      return false;
    }
  }
  return true;
}

OpenEarDecomposition schmidt(const AdjacencyMatrixGraph<Directionality::Undirected>& graph) {
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
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(graph.matrix(), top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());     // do not push already visited nodes
          tree.parent(it.index()) = top;  // update the parent
        }
      }
      tree.index(top) = indexCounter;
      ++indexCounter;
    }
  }

  // compute G \ T
  AdjacencyListGraph<Directionality::Directed> backedges;

  // decompose by iterating over every backedge outgoing from every node
  std::vector<std::vector<size_t>> ears;
  for (size_t v = 0; v < numberOfNodes; ++v) {
    // for every backedge starting at v
  }
}
