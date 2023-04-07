#include "solve/approximation.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/algorithm.hpp"
#include "graph/graph.hpp"
#include "graph/utils.hpp"

#include "solve/commonfunctions.hpp"

namespace approximation {

struct GraphPair {
  graph::AdjacencyListDigraph digraph;
  graph::AdjacencyListGraph graph;
};

static std::unordered_set<size_t> reduceToGminus(const graph::AdjacencyListDigraph& digraph, graph::AdjacencyListGraph& graph) {
  graph::AdjacencyListDigraph reverseDirgaph(digraph.numberOfNodes());
  for (const graph::Edge& e : digraph.edges()) {
    reverseDirgaph.addEdge(e.reverse());
  }

  std::unordered_set<size_t> cuttedNodes;
  for (const size_t u : reverseDirgaph.nodes()) {
    if (reverseDirgaph.degree(u) == 2) {
      const size_t v = reverseDirgaph.neighbours(u)[0];
      const size_t w = reverseDirgaph.neighbours(u)[1];
      graph.removeEdge(u, v);
      graph.removeEdge(u, w);
      graph.addEdge(v, w);
      cuttedNodes.insert(u);
    }
  }
  return cuttedNodes;
}

static void insertNodecuts(std::vector<size_t>& eulertourInGMinus, std::unordered_set<size_t>& cuttedNodes,
                           const graph::AdjacencyListDigraph& digraph) {
  eulertourInGMinus.reserve(eulertourInGMinus.size() + cuttedNodes.size());
  for (size_t i = eulertourInGMinus.size() - 1; i > 0; --i) {
    const size_t v = eulertourInGMinus[i];
    const size_t w = eulertourInGMinus[i - 1];
    for (const size_t u : digraph.neighbours(v)) {
      if (cuttedNodes.contains(u) && digraph.adjacent(w, u)) {
        eulertourInGMinus.insert(eulertourInGMinus.begin() + i, u);
        cuttedNodes.erase(u);
        break;
      }
    }
  }
}

static std::vector<size_t> findEulertour(graph::AdjacencyListGraph& graph, const graph::AdjacencyListDigraph& digraph) {
  std::unordered_set<size_t> cuttedNodes = reduceToGminus(digraph, graph);
  std::vector<size_t> eulertourInGMinus  = eulertour(graph);
  insertNodecuts(eulertourInGMinus, cuttedNodes, digraph);
  return eulertourInGMinus;
}

static std::vector<unsigned int> shortcutToHamiltoncycle(const std::vector<unsigned int>& longEulertour,
                                                         graph::AdjacencyListDigraph& digraph) {
  std::vector<unsigned int> hamiltoncycle;
  hamiltoncycle.reserve(digraph.numberOfNodes() + 1);
  hamiltoncycle.push_back(longEulertour[0]);

  for (size_t i = 1; i < longEulertour.size() - 1; ++i) {
    const size_t u = longEulertour[i - 1];
    const size_t w = longEulertour[i];
    const size_t v = longEulertour[i + 1];
    if (digraph.adjacent(w, u) && digraph.adjacent(w, v) && u != v) {
      hamiltoncycle.push_back(v);
      digraph.removeEdge(w, u);  // prevent the node from beeing skipped again between the same 2 nodes
      digraph.removeEdge(w, v);
      ++i;  // do not consider next node for skipping. It's alredy added.
    }
    else {
      hamiltoncycle.push_back(w);
    }
  }
  return hamiltoncycle;
}

static GraphPair constructGraphPair(const graph::EarDecomposition& ears, const size_t numberOfNodes) {
  graph::AdjacencyListGraph graph = earDecompToAdjacencyListGraph(ears, numberOfNodes);
  graph::AdjacencyListDigraph digraph(numberOfNodes);

  // in the first (or last) ear there is never an edge to be deleted
  const std::vector<size_t>& firstEar = ears.ears.back();
  for (size_t i = 0; i < firstEar.size() - 2; ++i) {
    digraph.addEdge(firstEar[i], firstEar[i + 1]);
  }
  digraph.addEdge(firstEar.back(), firstEar[firstEar.size() - 2]);

  // the other ears deletion of at most one edge can occur
  for (long j = ears.ears.size() - 2; j >= 0; --j) {
    const std::vector<size_t>& ear = ears.ears[j];

    digraph.addEdge(ear[0], ear[1]);  // add the first edge directing into the ear
    size_t earPosOfLastDoubledEdge = 0;

    struct EdgeIndex {
      graph::Edge e;
      size_t index;
    };
    std::vector<EdgeIndex> edgesToBeDirected;
    for (size_t i = 1; i < ear.size() - 1; ++i) {
      const size_t u = ear[i];
      const size_t v = ear[i + 1];

      if (graph.degree(u) % 2 == 1) {
        graph.addEdge(u, v);
        digraph.addEdge(u, v);
        digraph.addEdge(v, u);
        earPosOfLastDoubledEdge = i;
      }
      else {
        edgesToBeDirected.push_back(EdgeIndex{
            graph::Edge{u, v},
            i
        });
      }
    }

    // implicitly detect a node y and direct the edges according to paper
    if (earPosOfLastDoubledEdge != 0) {
      const graph::Edge lastDoubledEdge{ear[earPosOfLastDoubledEdge], ear[earPosOfLastDoubledEdge + 1]};
      graph.removeEdge(lastDoubledEdge);
      graph.removeEdge(lastDoubledEdge);

      digraph.removeEdge(lastDoubledEdge);
      digraph.removeEdge(lastDoubledEdge.reverse());

      for (const EdgeIndex& edge : edgesToBeDirected) {
        if (edge.index < earPosOfLastDoubledEdge) {
          digraph.addEdge(edge.e);
        }
        else {
          // the case of equality is implicitly excluded, because then the edge would have been doubled
          digraph.addEdge(edge.e.reverse());
        }
      }
    }
    else {
      for (const EdgeIndex& edge : edgesToBeDirected) {
        if (edge.index != ear.size() - 2) {
          digraph.addEdge(edge.e);
        }
        else {
          digraph.addEdge(edge.e.reverse());
        }
      }
    }
  }

  return GraphPair{digraph, graph};
}

Result approximateBTSP(const graph::Euclidean& euclidean, const bool printInfo) {
  double maxEdgeWeight;
  std::vector<unsigned int> tour, longEulertour;
  const graph::AdjacencyMatrixGraph biconnectedGraph = biconnectedSubgraph(euclidean, maxEdgeWeight);
  const graph::EarDecomposition ears                 = schmidt(biconnectedGraph);
  const graph::AdjacencyListGraph fromEars           = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());
  const graph::AdjacencyListGraph minimal            = minimallyBiconnectedSubgraph(fromEars);
  const graph::EarDecomposition openEars             = schmidt(minimal);  // calculate proper ear decomposition

  if (openEars.ears.size() == 1) {
    tour = std::vector<unsigned int>(openEars.ears[0].begin(), openEars.ears[0].end() - 1);  // do not repeat first node
  }
  else {
    GraphPair graphpair           = constructGraphPair(openEars, euclidean.numberOfNodes());
    const std::vector<size_t> tmp = findEulertour(graphpair.graph, graphpair.digraph);
    tour                          = shortcutToHamiltoncycle(std::vector<unsigned int>(tmp.begin(), tmp.end()), graphpair.digraph);
  }
  const graph::Edge bottleneckEdge = findBottleneck(euclidean, tour, true);
  const double objective           = euclidean.weight(bottleneckEdge);

  if (printInfo) {
    std::cout << "objective           : " << objective << std::endl;
    std::cout << "lower bound on OPT  : " << maxEdgeWeight << std::endl;
    std::cout << "a fortiori guarantee: " << objective / maxEdgeWeight << std::endl;
    assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
  }

  return Result{biconnectedGraph, openEars, tour, objective, bottleneckEdge};
}

static bool isCopyOf(const size_t node, const size_t compare, const size_t numberOfNodes) {
  return (node % numberOfNodes == compare && node < 5 * numberOfNodes);
}

Result approximateBTSPP(const graph::Euclidean& euclidean, const size_t s, const size_t t, const bool printInfo) {
  const graph::Edge st_Edge{s, t};

  double maxEdgeWeight;

  // find graph s.t. G = (V,E) + (s,t) is biconnected
  const graph::AdjacencyMatrixGraph biconnectedGraph = edgeAugmentedBiconnectedSubgraph(euclidean, st_Edge, maxEdgeWeight);
  const graph::EarDecomposition ears                 = schmidt(biconnectedGraph);
  graph::AdjacencyListGraph fromEars                 = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());
  if (!fromEars.adjacent(s, t)) {  // if the s-t edge is one of the removed ones,
    fromEars.addEdge(st_Edge);     // add it again.
  }
  graph::AdjacencyListGraph minimal = edgeKeepingMinimallyBiconectedSubgraph(fromEars, st_Edge);
  minimal.removeEdge(st_Edge);

  // create a graph consisting of 5 times the graph + nodes x and y connected to all copies of s resp. t
  const size_t numberOfNodes = euclidean.numberOfNodes();
  std::vector<std::vector<size_t>> adjacencyList;
  adjacencyList.reserve(numberOfNodes);

  for (size_t i = 0; i < 5; ++i) {
    adjacencyList.insert(adjacencyList.end(), minimal.adjacencyList().begin(), minimal.adjacencyList().end());
  }

  const size_t x = 5 * numberOfNodes;
  const size_t y = 5 * numberOfNodes + 1;
  adjacencyList.emplace_back(std::vector<size_t>{});  // add empty adjacency vector for x
  adjacencyList.emplace_back(std::vector<size_t>{});  // add empty adjacency vector for y
  graph::AdjacencyListGraph fiveFoldGraph(adjacencyList);

  // connect all copies of s to x and all copies of t to y
  for (size_t i = 0; i < 5; i++) {
    fiveFoldGraph.addEdge(i * numberOfNodes + s, x);
    fiveFoldGraph.addEdge(i * numberOfNodes + t, y);
  }

  const graph::EarDecomposition openEars = schmidt(minimal);  // calculate proper ear decomposition
  std::vector<unsigned int> tour, longEulertour;

  if (openEars.ears.size() == 1) {
    tour = std::vector<unsigned int>(openEars.ears[0].begin(), openEars.ears[0].end() - 1);  // do not repeat first node
  }
  else {
    GraphPair graphpair           = constructGraphPair(openEars, euclidean.numberOfNodes());
    const std::vector<size_t> tmp = findEulertour(graphpair.graph, graphpair.digraph);
    tour                          = shortcutToHamiltoncycle(std::vector<unsigned int>(tmp.begin(), tmp.end()), graphpair.digraph);
  }

  // extract s-t-path from solution
  const size_t pos_x = std::distance(tour.begin(), std::find(tour.begin(), tour.end(), x));
  const size_t pos_y = std::distance(tour.begin(), std::find(tour.begin(), tour.end(), y));

  // std::array<bool, 5> graphCopyIsSolution{true, true, true, true, true};
  // graphCopyIsSolution[]

  for (size_t i = 0; i < tour.size(); ++i) {
    if (isCopyOf(tour[i], s, numberOfNodes)) {
      if (graph::previousInCycle(tour, i) != x && graph::successiveInCycle(tour, i) != x) {
        // i is potential start node
      }
    }
  }

  const graph::Edge bottleneckEdge = findBottleneck(euclidean, tour, false);
  const double objective           = euclidean.weight(bottleneckEdge);

  if (printInfo) {
    std::cout << "objective           : " << objective << std::endl;
    std::cout << "lower bound on OPT  : " << maxEdgeWeight << std::endl;
    std::cout << "a fortiori guarantee: " << objective / maxEdgeWeight << std::endl;
    assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
  }

  return Result{biconnectedGraph, openEars, tour, objective, bottleneckEdge};
}
}  // namespace approximation
