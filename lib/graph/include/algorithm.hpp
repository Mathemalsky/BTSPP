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
#pragma once

#include <cassert>
#include <cstddef>
#include <stack>
#include <stdexcept>
#include <vector>

#include "exceptions.hpp"
#include "graph.hpp"

namespace graph {

/***********************************************************************************************************************
 *                                             types for graph algorithms
 **********************************************************************************************************************/

/*!
 * @brief EarDecomposition is a vector of ears (= vector of nodes)
 */
struct EarDecomposition {
  std::vector<std::vector<size_t>> ears; /**< list of ears */
};

/***********************************************************************************************************************
 *                                                  graph algorithms
 **********************************************************************************************************************/

/*!
 * @brief performs a depth first search on given graph with given node as root node
 * @tparam G type of graph
 * @param graph input graph
 * @param rootNode root node index
 * @return DfsTree providing the order the nodes are explored and the parent of each node
 */
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

/*!
 * \brief Creates a graph consisting of all edges which are not in the dfs tree.
 * \return AdjacencyListGraph with same nodes as dfs tree and complement of edges
 */
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
 * The algorithm is related to SCHMIDT, Jens M. A simple test on 2-vertex-and 2-edge-connectivity.
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
  size_t numberOfEars = 0;
  for (size_t u : tree.explorationOrder()) {    // iterate over all nodes in the order they appeared in dfs
    for (size_t v : backedges.neighbours(u)) {  // for every backedge starting at v
      if (!visited[v]) {                        // skip trivial backedges
        visited[u] = true;
        while (!visited[v]) {
          visited[v] = true;
          v          = tree.parent(v);
        }
        if (v == u && numberOfEars != 0) {
          return false;  // v is an articulation point
        }
        ++numberOfEars;
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
 * @brief finds a node with at least one neighbour
 * @tparam G type of graph
 * @param graph
 * @return index of non isolated node
 */
template <typename G>
size_t findNonIsolatedNode(const G& graph)
  requires(std::is_base_of_v<Graph, G>)
{
  for (size_t u = 0; u < graph.numberOfNodes(); ++u) {
    if (graph.degree(u) > 0) {
      return u;
    }
  }
  throw InfesableRequest("All nodes in graph are isolated!");
}

/*!
 * @brief finds an euler tour in graph using Hierholzers algorithm
 * @tparam G type of graph
 * @param graph an instance of a simple graph, all nodes are required to have even degree, all nodes with nonzero degree
 * must be connected (in fact 2-edge-connected)
 * @attention no check if input graph is indeed eulerian is performed
 * @return vector containing the nodes in the order they appear in the tour. The first node isn't added as last node
 * again.
 */
template <typename G>
std::vector<size_t> hierholzer(const G& graph)
  requires(std::is_base_of_v<Graph, G>)
{
  G workingCopy = graph;
  std::vector<size_t> tour;
  tour.reserve(graph.numberOfEdges() + 1);
  std::stack<size_t> nodeStack;
  nodeStack.push(findNonIsolatedNode(graph));

  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    if (workingCopy.degree(top) == 0) {
      tour.push_back(top);
      nodeStack.pop();
    }
    else {
      const size_t v = workingCopy.neighbourAny(top);
      workingCopy.removeEdge(Edge{top, v});  // remove an edge adjacent to top
      nodeStack.push(v);                     // put the node at the other end of the edge on the stack
    }
  }

  return tour;
}

/*!
 * @brief biconnectedSpanningGraph computes a bottleneck optimal biconnected subgraph.
 * @details First some edges definitely not increasing the bottleneck are added. Then the other edges are sortet
 * increasing in their length and successively added until the graph is biconnected.
 * @param euclidean complete graph, providing distances between nodes
 * @param maxEdgeWeight
 * @return undirected AdjacencyMatrixGraph
 */
AdjacencyMatrixGraph biconnectedSubgraph(const Euclidean& euclidean, double& maxEdgeWeight);

/*!
 * @brief edgeAugmentedBiconnectedSubGraph computes a bottle neck optimal subgraph, that is biconnected when adding the augemtation edge.
 * @details First the augmentation edge is added then some edges definitely not increasing the bottleneck are added. Then the other edges
 * are sortet increasing in their length and successively added until the graph is biconnected. Then the augmentation edge is removed.
 * @param euclidean complete graph, providing distances between nodes
 * @param augmentationEdge
 * @param maxEdgeWeight
 * @return undirected AdjacencyMatrixGraph
 */
AdjacencyMatrixGraph edgeAugmentedBiconnectedSubgraph(const Euclidean& euclidean, const Edge augmentationEdge, double& maxEdgeWeight);

/*!
 * \brief earDecompToAdjacencyListGraph puts all edges from ears into an undirected AdjacencyListGraph
 * \param earDecomposition open ear decomposition
 * \param numberOfNodes number of nodes appaering there
 * \return undirected AdjacencyListGraph containing numberOfNodes + number of ears - 1 edges.
 */
AdjacencyListGraph earDecompToAdjacencyListGraph(const EarDecomposition& earDecomposition, const size_t numberOfNodes);

/*!
 * @brief makes graph minimally biconnected
 * @details removes edges until all edges left are 2-essential
 * @param graph input graph
 * @return minimally biconnected subgraph of input graph
 */
AdjacencyListGraph minimallyBiconnectedSubgraph(const AdjacencyListGraph& graph);

/*!
 * @brief makes graph minimally biconnected but leaves given edge in the graph
 * @param graph input graph
 * @param keepEdge edge to keep
 * @return minimally biconnected graph the contains the edge keepEdge
 */
AdjacencyListGraph edgeKeepingMinimallyBiconectedSubgraph(const AdjacencyListGraph& graph, const Edge& keepEdge);

}  // namespace graph