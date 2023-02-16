#include "graph/graph.hpp"

#include <cassert>
#include <exception>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/algorithm.hpp"

/***********************************************************************************************************************
 *                                                adjacency list graph
 **********************************************************************************************************************/

bool AdjacencyListGraph::connected() const {
  std::vector<bool> visited(numberOfNodes(), false);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);
  while (!nodeStack.empty()) {
    const size_t v = nodeStack.top();
    nodeStack.pop();
    if (!visited[v]) {
      visited[v] = true;
      for (const size_t& neighbour : this->neighbours(v)) {
        if (!visited[neighbour]) {
          nodeStack.push(neighbour);
        }
      }
    }
  }
  for (const bool& node : visited) {
    if (!node) {
      return false;
    }
  }
  return true;
}

AdjacencyListGraph AdjacencyListGraph::removeUncriticalEdges() const {
  AdjacencyListGraph saveCopy         = *this;
  AdjacencyListGraph experimentalCopy = *this;

  for (const Edge& e : this->edgesToLowerIndex()) {
    experimentalCopy.removeEdge(e);
    if (experimentalCopy.biconnected()) {
      saveCopy = experimentalCopy;
    }
    else {
      experimentalCopy = saveCopy;
    }
  }
  return saveCopy;
}

/***********************************************************************************************************************
 *                                               adjacency list digraph
 **********************************************************************************************************************/

AdjacencyListGraph AdjacencyListDigraph::undirected() const {
  const size_t numberOfNodes = this->numberOfNodes();
  AdjacencyListGraph undirected(numberOfNodes);
  for (size_t i = 0; i < numberOfNodes; ++i) {
    for (size_t j = 0; j < degree(i); ++j) {
      if (std::find(undirected.neighbours(i).begin(), undirected.neighbours(i).end(), j)
          == undirected.neighbours(i).end()) {
        undirected.addEdge(i, j);
      }
    }
  }
  return undirected;
}

/***********************************************************************************************************************
 *                                               adjacency matrix graph
 **********************************************************************************************************************/

template <>
AdjacencyMatrixGraph<Directionality::Undirected> AdjacencyMatrixGraph<Directionality::Directed>::undirected() const {
  return AdjacencyMatrixGraph<Directionality::Undirected>(
      pAdjacencyMatrix + Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>(pAdjacencyMatrix.transpose()));
}

template <>
AdjacencyMatrixGraph<Directionality::Undirected> AdjacencyMatrixGraph<Directionality::Undirected>::undirected() const {
  assert("You are converting an undirected graph into an undirected graph!");
  return *this;
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::connected() const {
  std::vector<bool> visited(numberOfNodes(), false);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      for (const size_t u : this->neighbours(top)) {
        if (!visited[u]) {
          nodeStack.push(u);
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
bool AdjacencyMatrixGraph<Directionality::Directed>::connected() const {
  return undirected().connected();
}

template <>
bool AdjacencyMatrixGraph<Directionality::Undirected>::biconnected() const {
  return checkBiconnectivity(*this);
}
