#include "graph/graph.hpp"

#include <exception>
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
  DfsTree tree(numberOfNodes);  // graph with adjacency list
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
  for (size_t k = 0; k < numberOfNodes; ++k) {
    for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(graph.matrix(), k); it; ++it) {
      const size_t u = std::max(k, static_cast<size_t>(it.index()));
      const size_t v = std::min(k, static_cast<size_t>(it.index()));
      if (!tree.adjacent(u, v)) {  // if edge (k, it.index()) not in T
        backedges.addEdge(v, u);   // add the reverse directed edge to backedges
      }
    }
  }

  // decompose by iterating over every backedge outgoing from every node
  visited = std::vector<bool>(numberOfNodes, false);
  std::vector<std::vector<size_t>> ears;
  for (size_t v = 0; v < numberOfNodes; ++v) {
    // for every backedge starting at v
    for (const size_t& u : backedges.neighbours(v)) {
      if (!visited[u]) {
        std::vector<size_t> chain{v, u};
        visited[v] = true;
        size_t x   = u;
        while (!visited[x]) {
          chain.push_back(tree.parent(x));
          visited[x] = true;
          x          = tree.parent(x);
        }
        if (x == v && v != 0) {
          throw std::runtime_error("Graph is not biconnected!");
        }
        ears.push_back(chain);
      }
    }
  }
  return OpenEarDecomposition{ears};
}
