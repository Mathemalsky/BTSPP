#include "graph/graph.hpp"

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

bool UndirectedGraph::biconnected() const {
  for (size_t i = 0; i < pNumberOfNodes; ++i) {
    Eigen::SparseMatrix<EdgeCost> adjacencyMatrix = pAdjacencyMatrix;
    // remove one column and one row from the matrix
    // check if the graph is still connected
  }
}
