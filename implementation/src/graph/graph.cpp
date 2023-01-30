#include "graph/graph.hpp"

#include <exception>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connected() const {
  std::vector<bool> visited(pNumberOfNodes, false);
  std::stack<size_t> nodeStack;
  const size_t rootNode = 0;
  nodeStack.push(rootNode);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      for (unsigned int k = 0; k < top; ++k) {
        if (pAdjacencyMatrix.coeff(k, top) != 0.0) {
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

// DEBUG
#include <iostream>
#include "utility/utils.hpp"

/*
template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connectedWhithout(const size_t vertex) const {
  // DEBUG
  std::cerr << "#nodes in connectedWithout: " << pNumberOfNodes << std::endl;
  std::vector<bool> component(pNumberOfNodes, false);
  component[0] = true;         // without loss of generality we consider the connected component containing node 0
  component[1] = vertex == 0;  // if we check connectivity whithout 0, we start by one

  for (unsigned int k = 0; k < pAdjacencyMatrix.outerSize(); ++k) {
    if (k != vertex) {
      if (!component[k]) {
        for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, k); it; ++it) {
          if (component[it.index()]) {
            component[k] = true;
            break;
          }
        }

        if (!component[k]) {
          // DBEUG
          std::cerr << "#smallest not connected index: " << k << std::endl;

          return false;  // node k is not connected to any node < k, because adjacency matrix is lower triangular
        }
      }
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, k); it; ++it) {
        component[it.index()] = true;
      }
    }
  }
  return true;
}
*/

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connectedWhithout(const size_t vertex) const {
  std::vector<bool> visited(pNumberOfNodes, false);
  std::stack<size_t> nodeStack;
  const size_t rootNode = (vertex != 0 ? 0 : 1);
  nodeStack.push(rootNode);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();

    // DEBUG
    std::cerr << "top node: " << top << std::endl;

    if (!visited[top]) {
      for (unsigned int k = 0; k < top; ++k) {
          // DEBUG
          std::cerr << "(" << k << ", " << top << ") " << pAdjacencyMatrix.coeff(top, k).cost() << std::endl;
        if (pAdjacencyMatrix.coeff(top, k) != 0.0) {
          nodeStack.push(k);

          // DEBUG
          std::cerr << "added node: " << k << " (row iterating)" << std::endl;
        }
      }
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(pAdjacencyMatrix, top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());

          // DEBUG
          std::cerr << "added node: " << it.index() << " (column iterating)" << std::endl;
        }
      }
      visited[top] = true;
    }

    // DEBUG
    std::cerr << "stack size: " << nodeStack.size() << std::endl;
  }

  std::cerr << visited;

  for (const bool& test : visited) {
    if (!test) {
      return false;
    }
  }
  return true;
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::biconnected() const {
  // DEBUG
  std::cerr << "#nodes in biconnected: " << pNumberOfNodes << std::endl;

  for (size_t i = 0; i < pNumberOfNodes; ++i) {
    if (!connectedWhithout(i)) {
      // DEBUG
      std::cerr << "cut vertex: " << i << std::endl;

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

AdjacencyListGraph<Directionality::Directed> findBackedges(
    const AdjacencyMatrixGraph<Directionality::Undirected>& graph, const DfsTree& tree) {
  AdjacencyListGraph<Directionality::Directed> backedges;
  for (size_t k = 0; k < graph.numberOfNodes(); ++k) {
    for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(graph.matrix(), k); it; ++it) {
      const size_t u = std::max(k, static_cast<size_t>(it.index()));
      const size_t v = std::min(k, static_cast<size_t>(it.index()));
      if (!tree.adjacent(u, v)) {  // if edge (k, it.index()) not in T
        backedges.addEdge(v, u);   // add the reverse directed edge to backedges
      }
    }
  }
  return backedges;
}

OpenEarDecomposition schmidt(const AdjacencyMatrixGraph<Directionality::Undirected>& graph) {
  const DfsTree tree                                           = dfs(graph);
  const AdjacencyListGraph<Directionality::Directed> backedges = findBackedges(graph, tree);

  // decompose by iterating over every backedge outgoing from every node
  const size_t numberOfNodes = graph.numberOfNodes();
  std::vector<bool> visited(numberOfNodes, false);
  std::vector<std::vector<size_t>> ears;
  for (size_t v = 0; v < numberOfNodes; ++v) {
    for (const size_t& u : backedges.neighbours(v)) {  // for every backedge starting at v
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
