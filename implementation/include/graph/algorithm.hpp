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

template <typename G>
EarDecomposition schmidt(const G& graph) {
  const DfsTree tree                                             = dfs(graph);
  const AdjacencyListGraph<Directionality::Undirected> backedges = findBackedges(graph, tree);
  const size_t numberOfNodes                                     = graph.numberOfNodes();

  // DEBUG
  std::cerr << "tree called in schmidt\n" << tree.explorationOrder() << std::endl;

  std::vector<bool> visited(numberOfNodes, false);
  std::vector<std::vector<size_t>> ears;
  std::vector<size_t> articulationPoints;
  for (size_t v : tree.explorationOrder()) {    // iterate over all nodes in the order they appeared in dfs
    for (size_t u : backedges.neighbours(v)) {  // for every backedge starting at v
      if (!visited[u]) {
        std::vector<size_t> chain{v, u};
        visited[v] = true;
        while (!visited[u]) {
          visited[u] = true;
          u          = tree.parent(u);
          chain.push_back(u);
        }
        if (u == v && v != tree.root()) {
          articulationPoints.push_back(u);
        }
        ears.push_back(chain);
      }
    }
  }
  return EarDecomposition{ears, articulationPoints};
}

AdjacencyMatrixGraph<Directionality::Undirected> earDecompToGraph(const EarDecomposition& earDecomposition);
AdjacencyMatrixGraph<Directionality::Undirected> biconnectedSpanningGraph(const Euclidean& euclidean);
AdjacencyListGraph<Directionality::Undirected> earDecompToAdjacencyListGraph(
    const EarDecomposition& earDecomposition, const size_t numberOfNodes);
