#pragma once

#include "graph/graph.hpp"

// DEBUG
#include <iostream>
#include "utility/utils.hpp"
#include "graph/ostream.hpp"

template <typename G>
DfsTree dfs(const G& graph, const size_t rootNode = 0) {
  const size_t numberOfNodes = graph.numberOfNodes();

  DfsTree tree(numberOfNodes);
  std::vector<bool> visited(numberOfNodes, false);
  std::stack<size_t> nodeStack;
  nodeStack.push(rootNode);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      visited[top] = true;
      tree.explorationOrder().push_back(top);  // store order of node exploration
      for (const size_t u : graph.neighbours(top)) {
        if (!visited[u]) {
          nodeStack.push(u);     // do not push already visited nodes
          tree.parent(u) = top;  // update the parent
        }
      }
    }
  }
  return tree;
}

template <typename G>
AdjacencyListGraph<Directionality::Undirected> findBackedges(const G& graph, const DfsTree& tree) {
  AdjacencyListGraph<Directionality::Undirected> backedges(graph.numberOfNodes());
  for (Edge e : graph.edges()) {
    if (!tree.adjacent(e) && !tree.adjacent(e.reverse())) {
      backedges.addEdge(e);
    }
  }
  return backedges;
}

/*!
 * \brief schmidt computes an open ear decomposition
 * \details This function requires the input graph to be biconnected (strongly biconnected in case of directed graph).
 * The algorithm is related to
 * SCHMIDT, Jens M. A simple test on 2-vertex-and 2-edge-connectivity.
 * Information Processing Letters, 2013, 113. Jg., Nr. 7, S. 241-244.
 * \param graph Instance of a graph, which provides the required member functions.
 * \return Open Ear decomposition as std::vector of ears which are std::vector<size_t>.
 */
template <typename G>
EarDecomposition schmidt(const G& graph) {
  assert(checkBiconnectivity(graph) && "The graph is assumed to be biconnected!");

  const DfsTree tree                                             = dfs(graph);
  const AdjacencyListGraph<Directionality::Undirected> backedges = findBackedges(graph, tree);
  const size_t numberOfNodes                                     = graph.numberOfNodes();

  std::vector<bool> visited(numberOfNodes, false);
  std::vector<std::vector<size_t>> ears;
  for (size_t u : tree.explorationOrder()) {    // iterate over all nodes in the order they appeared in dfs
    for (size_t v : backedges.neighbours(u)) {  // for every backedge starting at v
      if (!visited[v]) {
        std::vector<size_t> chain{u, v};
        visited[u] = true;
        while (!visited[v]) {
          visited[v] = true;
          v          = tree.parent(v);
          chain.push_back(v);
        }
        ears.push_back(chain);
      }
    }
  }
  return EarDecomposition{ears};
}

template <typename G>
bool checkBiconnectivity(const G& graph) {
  const DfsTree tree                                             = dfs(graph);
  const AdjacencyListGraph<Directionality::Undirected> backedges = findBackedges(graph, tree);
  const size_t numberOfNodes                                     = graph.numberOfNodes();

  if (backedges.degree(tree.root()) == 0) {
    return false;  // root isn't part of any cycle
  }

  std::vector<bool> visited(numberOfNodes, false);
  for (size_t u : tree.explorationOrder()) {    // iterate over all nodes in the order they appeared in dfs
    for (size_t v : backedges.neighbours(u)) {  // for every backedge starting at v
      if (!visited[v]) {
        visited[u] = true;
        while (!visited[v]) {
          visited[v] = true;
          v          = tree.parent(v);
        }
        if (v == u && u != tree.root()) {
          return false;  // v is an articulation point
        }
      }
    }
  }

  for (const bool b : visited) {
    if (!b) {
      return false;
    }
  }
  return true;
}

AdjacencyMatrixGraph<Directionality::Undirected> biconnectedSpanningGraph(const Euclidean& euclidean);
AdjacencyListGraph<Directionality::Undirected> earDecompToAdjacencyListGraph(
    const EarDecomposition& earDecomposition, const size_t numberOfNodes);
