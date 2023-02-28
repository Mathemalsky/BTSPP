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

#include "graph/graph.hpp"
#include "graph/algorithm.hpp"

#include "solve/commonfunctions.hpp"

namespace approximation {

static std::unordered_set<size_t> reduceToGminus(const AdjacencyListDigraph& digraph, AdjacencyListGraph& graph) {
  AdjacencyListDigraph reverseDirgaph(digraph.numberOfNodes());
  for (const Edge& e : digraph.edges()) {
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

static void insertNodecuts(
    std::vector<size_t>& eulertourInGMinus, std::unordered_set<size_t>& cuttedNodes,
    const AdjacencyListDigraph& digraph) {
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

static std::vector<size_t> findEulertour(AdjacencyListGraph& graph, const AdjacencyListDigraph& digraph) {
  std::unordered_set<size_t> cuttedNodes = reduceToGminus(digraph, graph);
  std::vector<size_t> eulertourInGMinus  = eulertour(graph);
  insertNodecuts(eulertourInGMinus, cuttedNodes, digraph);
  return eulertourInGMinus;
}

static std::vector<unsigned int> shortcutToHamiltoncycle(
    const std::vector<size_t>& longEulertour, AdjacencyListDigraph& digraph) {
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

static std::vector<unsigned int> hamiltonCycleInSquare(const EarDecomposition& ears, const size_t numberOfNodes) {
  if (ears.ears.size() == 1) {
    // cast down to unsigned int and cut off last node which is same as first
    return std::vector<unsigned int>(ears.ears[0].begin(), ears.ears[0].end() - 1);
  }
  else {
    AdjacencyListGraph graph = earDecompToAdjacencyListGraph(ears, numberOfNodes);
    AdjacencyListDigraph digraph(numberOfNodes);

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
        Edge e;
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
          edgesToBeDirected.push_back(EdgeIndex{Edge{u, v}, i});
        }
      }

      // implicitly detect a node y and direct the edges according to paper
      if (earPosOfLastDoubledEdge != 0) {
        const Edge lastDoubledEdge{ear[earPosOfLastDoubledEdge], ear[earPosOfLastDoubledEdge + 1]};
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

    const std::vector<size_t> longEulertour = findEulertour(graph, digraph);
    return shortcutToHamiltoncycle(longEulertour, digraph);
  }
}

Result approximate(const Euclidean& euclidean, const ProblemType problemType, const bool printInfo) {
  if (problemType == ProblemType::BTSP_approx) {
    double maxEdgeWeight;
    const AdjacencyMatrixGraph biconnectedGraph = biconnectedSpanningGraph(euclidean, maxEdgeWeight);
    const EarDecomposition ears                 = schmidt(biconnectedGraph);
    const AdjacencyListGraph fromEars           = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());
    const AdjacencyListGraph minimal            = fromEars.removeUncriticalEdges();
    const EarDecomposition openEars             = schmidt(minimal);  // calculate proper ear decomposition
    const std::vector<unsigned int> tour        = hamiltonCycleInSquare(openEars, euclidean.numberOfNodes());
    const Edge bottleneckEdge                   = findBottleneck(euclidean, tour, true);
    const double objective                      = euclidean.weight(bottleneckEdge);

    if (printInfo) {
      std::cout << "objective           : " << objective << std::endl;
      std::cout << "lower bound on OPT  : " << maxEdgeWeight << std::endl;
      std::cout << "a fortiori guarantee: " << objective / maxEdgeWeight << std::endl;
      assert(objective / maxEdgeWeight <= 2 && objective / maxEdgeWeight >= 1 && "A fortiori guarantee is nonsense!");
    }

    return Result{biconnectedGraph, openEars, tour, objective, bottleneckEdge};
  }
  else {
    throw std::runtime_error("Unknown problem type");
  }
}
}  // namespace approximation
