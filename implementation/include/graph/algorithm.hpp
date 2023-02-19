#pragma once

#include <cassert>
#include <cstddef>
#include <stack>
#include <stdexcept>
#include <vector>

#include "graph/graph.hpp"

/***********************************************************************************************************************
 *                                             types for graph algorithms
 **********************************************************************************************************************/

struct EarDecomposition {
  std::vector<std::vector<size_t>> ears;
};

/***********************************************************************************************************************
 *                                                  graph algorithms
 **********************************************************************************************************************/

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
AdjacencyListGraph findBackedges(const G& graph, const DfsTree& tree) {
  AdjacencyListGraph backedges(graph.numberOfNodes());
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

  const DfsTree tree                 = dfs(graph);
  const AdjacencyListGraph backedges = findBackedges(graph, tree);
  const size_t numberOfNodes         = graph.numberOfNodes();

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

/*!
 * \brief checkBiconnectivity checks a graph is biconnected
 * \details This function works similar to schmidt() but does not track the ears and beforms checks for biconnectivity
 * instead. \param graph The graph object should provide the neighbours() function \return true if the graph is
 * biconnected, false otherwise.
 */
template <typename G>
bool checkBiconnectivity(const G& graph) {
  const DfsTree tree                 = dfs(graph);
  const AdjacencyListGraph backedges = findBackedges(graph, tree);
  const size_t numberOfNodes         = graph.numberOfNodes();

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

/*!
 * \brief biconnectedSpanningGraph computes a bottleneck optimal biconnected spanning subgraph.
 * \details First some edges definitely not increasing the bottleneck are added. Then the other edges are sortet
 * increasing in their length and successively added until the graph is biconnected.
 * \param euclidean complete graph, providing distances between nodes.
 * \return undirected AdjacencyMatrixGraph
 */
AdjacencyMatrixGraph biconnectedSpanningGraph(const Euclidean& euclidean);

/*!
 * \brief earDecompToAdjacencyListGraph puts all edges from ears into an undirected AdjacencyListGraph
 * \param earDecomposition open ear decomposition
 * \param numberOfNodes number of nodes appaering there
 * \return undirected AdjacencyListGraph containing numberOfNodes + number of ears - 1 edges.
 */
AdjacencyListGraph earDecompToAdjacencyListGraph(const EarDecomposition& earDecomposition, const size_t numberOfNodes);

template <typename G>
size_t findNonIsolatedNode(const G& graph) requires(std::is_base_of_v<Graph, G>) {
  for (size_t u = 0; u < graph.numberOfNodes(); ++u) {
    if (graph.degree(u) > 0) {
      return u;
    }
  }
  throw std::runtime_error("All nodes in graph are isolated!");
}

template <typename G>
std::vector<size_t> eulertour(const G& graph) requires(std::is_base_of_v<UndirectedGraph, G>) {
  G workingCopy = graph;
  size_t u      = findNonIsolatedNode(graph);
  std::vector<size_t> eulertour;
  eulertour.reserve(graph.numberOfEdges());
  std::stack<size_t> nodeStack;
  nodeStack.push(u);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    if (workingCopy.degree(top) == 0) {
      // do something
    }
    else {
      // remove an edge adjacent to top
      // put the node at the other end of the vertex on top of the stack
    }
  }
}
