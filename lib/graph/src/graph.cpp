/*
 * GRAPH is a library to store and manipulate graphs as adjacency list or
 * as sparse eigen matrix. Different specialized types of graphs are
 * supported.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "graph.hpp"

#include <cassert>
#include <exception>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

#include "algorithm.hpp"

namespace graph {

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

/***********************************************************************************************************************
 *                                               adjacency list digraph
 **********************************************************************************************************************/

AdjacencyListGraph AdjacencyListDigraph::undirected() const {
  std::vector<std::vector<size_t>> adjacencyList = pAdjacencyList;
  for (const Edge& edge : edges()) {
    adjacencyList[edge.v].push_back(edge.u);
  }
  return AdjacencyListGraph(adjacencyList);
}

/***********************************************************************************************************************
 *                                               adjacency matrix graph
 **********************************************************************************************************************/

bool AdjacencyMatrixGraph::connected() const {
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

/***********************************************************************************************************************
 *                                              adjacency matrix digraph
 **********************************************************************************************************************/

AdjacencyMatrixGraph AdjacencyMatrixDigraph::undirected() const {
  return AdjacencyMatrixGraph(pAdjacencyMatrix + Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>(pAdjacencyMatrix.transpose()));
}
}  // namespace graph
