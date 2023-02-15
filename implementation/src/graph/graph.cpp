#include "graph/graph.hpp"

#include <cassert>
#include <exception>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/algorithm.hpp"

template <>
AdjacencyListGraph<Directionality::Undirected> AdjacencyListGraph<Directionality::Directed>::undirected() const {
  const size_t numberOfNodes = this->numberOfNodes();
  AdjacencyListGraph<Directionality::Undirected> undirected(numberOfNodes);
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

template <>
AdjacencyListGraph<Directionality::Undirected> AdjacencyListGraph<Directionality::Undirected>::undirected() const {
  assert("You are converting an undirected graph into an undirected graph!");
  return *this;
}

template <>
bool AdjacencyListGraph<Directionality::Undirected>::connected() const {
  std::vector<bool> visited(numberOfNodes(), false);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);
  while (!nodeStack.empty()) {
    const size_t v = nodeStack.top();
    nodeStack.pop();
    if (!visited[v]) {
      visited[v] = true;
      for (const size_t& neighbour : pAdjacencyList[v]) {
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

template <>
bool AdjacencyListGraph<Directionality::Directed>::connected() const {
  return undirected().connected();
}

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
      for (Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>::InnerIterator it(pAdjacencyMatrix, top); it; ++it) {
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
bool AdjacencyMatrixGraph<Directionality::Directed>::connected() const {
  return undirected().connected();
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

template <>
bool AdjacencyListGraph<Directionality::Undirected>::biconnected() const {
  return checkBiconnectivity(*this);
}

// DEBUG
#include <iostream>
#include "graph/ostream.hpp"

/*
template <>
AdjacencyMatrixGraph<Directionality::Undirected>
    AdjacencyMatrixGraph<Directionality::Undirected>::removeUncriticalEdges() const {
  AdjacencyMatrixGraph<Directionality::Undirected> saveCopy        = *this;
  AdjacencyMatrixGraph<Directionality::Undirected> experimetalCopy = *this;
  for (const Edge& e : this->edgesToLowerIndex()) {
    experimetalCopy.removeEdge(e);
    experimetalCopy.prune();
    if (schmidt(experimetalCopy).open()) {
      saveCopy = experimetalCopy;
    }
    else {
      experimetalCopy = saveCopy;
    }
  }
  return saveCopy;
}
*/

template <>
AdjacencyListGraph<Directionality::Undirected> AdjacencyListGraph<Directionality::Undirected>::removeUncriticalEdges()
    const {
  AdjacencyListGraph<Directionality::Undirected> saveCopy         = *this;
  AdjacencyListGraph<Directionality::Undirected> experimentalCopy = *this;
  for (const Edge& e : this->edgesToLowerIndex()) {
    experimentalCopy.removeEdge(e);

    // DEBUG
    std::cerr << "experimental copy\n" << experimentalCopy << std::endl;
    std::cerr << "save         copy\n" << saveCopy << std::endl;

    if (experimentalCopy.biconnected()) {
      saveCopy = experimentalCopy;
    }
    else {
      experimentalCopy = saveCopy;
    }
  }
  return saveCopy;
}
