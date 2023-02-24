#include "solve/approximation.hpp"

#include <cassert>
#include <cmath>
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

// DEBUG
#include <iostream>
#include "utility/utils.hpp"
#include "graph/ostream.hpp"

namespace approximation {

/*
static std::vector<unsigned int> shortcutToHamiltoncycle(
    const std::vector<size_t>& eulertourInEarDecomposition, const size_t numberOfNodes) {
  std::vector<bool> visited(numberOfNodes, false);
  std::vector<unsigned int> hamiltoncycle;
  hamiltoncycle.reserve(numberOfNodes);

  for (const size_t u : eulertourInEarDecomposition) {
    if (!visited[u]) {
      hamiltoncycle.push_back(static_cast<unsigned int>(u));
      visited[u] = true;
    }
  }
  return hamiltoncycle;
}
*/

/*
static std::vector<size_t> backAndForth(const std::vector<size_t>& vecIn) {
  std::vector<size_t> vecOut;
  vecOut.reserve(vecIn.size());

  for (size_t i = 0; i <= (vecIn.size() - 1) / 2; ++i) {
    vecOut.push_back(vecIn[2 * i]);  // even positions in ascending order
  }
  for (size_t i = vecIn.size() / 2; i > 0; --i) {
    vecOut.push_back(vecIn[2 * i - 1]);  // odd positions in descending order
  }

  return vecOut;
}
*/

/*
static std::vector<unsigned int> assembleTour(
    const std::vector<size_t>& hamiltonSubcycle, const std::vector<std::vector<size_t>>&
earPositionsAdjacentToDoubleEdges, const size_t numberOfNodes) { std::unordered_map<size_t, size_t> doubleEdgeStarts;
  for (size_t i = 0; i < earPositionsAdjacentToDoubleEdges.size(); ++i) {
    doubleEdgeStarts.insert({earPositionsAdjacentToDoubleEdges[i][0], i});
  }

  std::vector<unsigned int> hamiltoncycle;
  hamiltoncycle.reserve(numberOfNodes);

  for (const size_t u : hamiltonSubcycle) {
    hamiltoncycle.push_back(u);
    if (doubleEdgeStarts.contains(u)) {
      const std::vector<unsigned int> tmp = backAndForth(earPositionsAdjacentToDoubleEdges[doubleEdgeStarts[u]]);
      hamiltoncycle.insert(hamiltoncycle.end(), tmp.begin(), tmp.end());
      doubleEdgeStarts.erase(u);
    }
  }
  assert(doubleEdgeStarts.empty() && "Couldn't insert all subtours!");
  return hamiltoncycle;
}
*/

/*
static void resolveDoubleEdges(
    AdjacencyListGraph& graph, std::vector<std::vector<size_t>>& earPositionsAdjacentToDoubleEdges,
    const std::vector<size_t>& ear) {
  bool removedPair = false;
  for (std::vector<size_t>& doubleEdgeSegment : earPositionsAdjacentToDoubleEdges | std::views::reverse) {
    if (!removedPair) {
      // graph.removeEdge(doubleEdgeSegment[doubleEdgeSegment.size() - 2], doubleEdgeSegment.back());
      removedPair = true;
      if (doubleEdgeSegment.size() > 2) {  // if the last edge is not the only one in this segment
        doubleEdgeSegment.pop_back();


        // shortcut edge connecting successor of u in ear an v
        const size_t segmentStartIndex = doubleEdgeSegment[0];
        const size_t u                 = ear[segmentStartIndex];
        const size_t v                 = graph.neighbourAnyExcept(
            u, std::unordered_set<size_t>{ear[segmentStartIndex - 1], ear[segmentStartIndex + 1]});
        graph.removeEdge(u, v);
        graph.addEdge(ear[segmentStartIndex + 1], v);

        // add back and forth edges
        const std::vector<size_t> order = backAndForth(doubleEdgeSegment);
        for (size_t i = 0; i < order.size() - 1; ++i) {
          graph.addEdge(ear[order[i]], ear[order[i + 1]]);
        }
      }
    }
    else if (doubleEdgeSegment.size() == 2) {
      // shortcut consecutive edges in ear
      graph.removeEdge(ear[doubleEdgeSegment.back()], ear[doubleEdgeSegment.back() + 1]);
      graph.addEdge(ear[doubleEdgeSegment.front()], ear[doubleEdgeSegment.back() + 1]);

      // revert former deletion
      graph.addEdge(ear[doubleEdgeSegment.front()], ear[doubleEdgeSegment.back()]);
    }
    else {
      // remove all edges except for the last

      // shortcut edge connecting successor of u in ear an v
      const size_t segmentStartIndex = doubleEdgeSegment[0];
      const size_t u                 = ear[segmentStartIndex];
      const std::unordered_set<size_t> forbiddenNeighbours{ear[segmentStartIndex - 1], ear[segmentStartIndex + 1]};
      const size_t v = graph.neighbourAnyExcept(u, [&](const size_t w) { return forbiddenNeighbours.contains(w); });
      graph.removeEdge(u, v);
      graph.addEdge(ear[segmentStartIndex + 1], v);

      for (size_t i = 0; i < doubleEdgeSegment.size() - 2; ++i) {
        graph.addEdge(ear[doubleEdgeSegment[i]], ear[doubleEdgeSegment[i + 2]]);
      }

      // revert former deletion of last edge
      graph.addEdge(ear[doubleEdgeSegment.back() - 1], ear[doubleEdgeSegment.back()]);
    }
  }
}
*/

/*
static size_t yInEar(const AdjacencyListGraph& graph, const std::vector<size_t>& ear) {
  for (size_t i = 1; i < ear.size() - 1; ++i) {
    if (graph.degree(ear[i]) == 2) {
      return i;
    }
  }
  throw std::runtime_error("There is no degree 2 node in this ear!");
}
*/

/*!
 * @brief
 * @details Differs from eulertour in prefering edges with special properties
 * @param graph
 * @param digraph
 * @return
 */
static std::vector<size_t> findEulertour(const AdjacencyListGraph& graph, const AdjacencyListDigraph& digraph) {
  AdjacencyListGraph workingCopy = graph;
  std::vector<size_t> eulertour;
  eulertour.reserve(graph.numberOfEdges() + 1);
  std::stack<size_t> nodeStack;
  nodeStack.push(0);  // the graph is connected so we start at node 0

  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    if (workingCopy.degree(top) == 0) {
      eulertour.push_back(top);
      nodeStack.pop();
    }
    else {
      const size_t v = workingCopy.neighbourAnyPrefer(top, [&](const size_t v) { return digraph.adjacent(v, top); });
      workingCopy.removeEdge(Edge{top, v});  // remove an edge adjacent to top
      nodeStack.push(v);                     // put the node at the other end of the edge on the stack
    }
  }

  // eulertour.pop_back();  // cut off the last node which is same as first
  return eulertour;
}

static std::vector<unsigned int> shortcutToHamiltoncycle(
    const std::vector<size_t>& longEulertour, AdjacencyListDigraph& digraph, AdjacencyListGraph& eulertourGraph) {
  AdjacencyListGraph graph(digraph.numberOfNodes());

  // DEBUG
  std::cerr << "longEulertour\n" << longEulertour;
  // std::cerr << "digraph\n" << digraph << std::endl;

  for (size_t i = 1; i < longEulertour.size() - 1; ++i) {
    const size_t u = longEulertour[i - 1];
    const size_t w = longEulertour[i];
    const size_t v = longEulertour[i + 1];

    if (digraph.adjacent(w, u) && digraph.adjacent(w, v) && u != v) {
      // DEBUG
      // std::cerr << "u: " << u << "\nw: " << w << "\nv: " << v << "\n";

      // std::cerr << "digraph\n" << digraph << std::endl;

      eulertourGraph.removeEdge(w, u);
      eulertourGraph.removeEdge(w, v);
      eulertourGraph.addEdge(u, v);

      digraph.removeEdge(w, u);
      digraph.removeEdge(w, v);
    }
  }

  std::vector<size_t> hamiltoncycle = eulertour(eulertourGraph);
  return std::vector<unsigned int>(hamiltoncycle.begin(), hamiltoncycle.end());
}

// DEBUG
struct EdgeIndex {
  Edge e;
  size_t index;
};

/*
std::ostream& operator<<(std::ostream& os, const EdgeIndex& ei) {
  return os << ei.e << " index: " << ei.index;
}

template<typename Type>
void printVec(const std::vector<Type>& vec) {
  const size_t size = vec.size();
  for (size_t i = 0; i < size; ++i) {
    std::cerr << vec[i] << (i == size - 1 ? "" : " ");
  }
}
*/

static std::vector<unsigned int> hamiltonCycleInSquare(const EarDecomposition& ears, const size_t numberOfNodes) {
  if (ears.ears.size() == 1) {
    // cast down to unsigned int and cut off last node which is same as first
    return std::vector<unsigned int>(ears.ears[0].begin(), ears.ears[0].end() - 1);
  }
  else {
    AdjacencyListGraph graph = earDecompToAdjacencyListGraph(ears, numberOfNodes);
    AdjacencyListDigraph digraph(numberOfNodes);

    // first ear is trivial
    const std::vector<size_t>& firstEar = ears.ears.back();
    for (size_t i = 0; i < firstEar.size() - 2; ++i) {
      digraph.addEdge(firstEar[i], firstEar[i + 1]);
    }
    digraph.addEdge(firstEar.back(), firstEar[firstEar.size() - 2]);

    // the other ears are more complicated
    for (long j = ears.ears.size() - 2; j >= 0; --j) {
      const std::vector<size_t>& ear = ears.ears[j];
      // const size_t y                 = yInEar(graph, ear);

      digraph.addEdge(ear[0], ear[1]);  // add the first edge directing into the ear
      size_t earPosOfLastDoubledEdge = 0;

      std::vector<EdgeIndex> edgesToBeDirected;
      for (size_t i = 1; i < ear.size() - 1; ++i) {
        const size_t u = ear[i];
        const size_t v = ear[i + 1];

        // DEBUG
        /*
        if (i == ear.size() - 2) {
          std::cerr << "u: " << u << " degree(u): " << graph.degree(u) << std::endl;
          std::cerr << "v: " << v << std::endl;
        }
        */
        if (graph.degree(u) % 2 == 1) {
          graph.addEdge(u, v);
          digraph.addEdge(u, v);
          digraph.addEdge(v, u);
          earPosOfLastDoubledEdge = i;
        }
        else {
          edgesToBeDirected.push_back(EdgeIndex{Edge{u, v}, i});

          // DEBUG
          /*
          std::cerr << "edges to be directed\n";
          for(EdgeIndex& ei : edgesToBeDirected) {
            std::cerr << ei.e << " index: " << ei.index << std::endl;
          }*/
          // printVec(edgesToBeDirected);
          // std::cerr << Edge{u, v} << " added to edgesToBeDirected" << std::endl;
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

            // DEBUG
            // std::cerr << "should add reverse edge " << edge.e.reverse() << std::endl;
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

    const std::vector<size_t> eulertour = findEulertour(graph, digraph);
    return shortcutToHamiltoncycle(eulertour, digraph, graph);

    // DEPRECATED
    /*
    for (long i = ears.ears.size() - 2; i >= 0; --i) {
      const std::vector<size_t>& ear = ears.ears[i];
      bool doublePreviousEdge        = false;
      std::vector<std::vector<size_t>> earPositionsAdjacentToDoubleEdges;

      for (size_t i = 1; i < ear.size() - 1; ++i) {
        const size_t u = ear[i];
        const size_t v = ear[i + 1];
        if (graph.degree(u) % 2 == 1) {
          graph.removeEdge(Edge{u, v});
          if (!doublePreviousEdge) {
            earPositionsAdjacentToDoubleEdges.push_back(std::vector<size_t>{i, i + 1});
          }
          else {
            earPositionsAdjacentToDoubleEdges.back().push_back(i + 1);
          }
          doublePreviousEdge = true;
        }
        else {
          doublePreviousEdge = false;
        }
      }
      resolveDoubleEdges(graph, earPositionsAdjacentToDoubleEdges, ear);

      // DEBUG
      std::cerr << "earPositionsAdjacentToDoubleEdges\n " << earPositionsAdjacentToDoubleEdges;
    }

    // find euler tour
    std::vector<size_t> eulerSubtour = eulertour(graph);

    // DEBUG
    std::cerr << "eulerSubtour\n" << eulerSubtour;

    return shortcutToHamiltoncycle(eulerSubtour, numberOfNodes);

    // DEBUG
    // std::cerr << "hamiltonSubcycle\n" << hamiltonSubcycle;

    // return assembleTour(hamiltonSubcycle, earPositionsAdjacentToDoubleEdges, numberOfNodes);
    */
  }
}

Result approximate(const Euclidean& euclidean, const ProblemType problemType, const bool printInfo) {
  if (problemType == ProblemType::BTSP_approx) {
    double maxEdgeWeight;
    AdjacencyMatrixGraph biconnectedGraph = biconnectedSpanningGraph(euclidean, maxEdgeWeight);

    // DEBUG
    // std::cerr << "undirected\n" << biconnectedGraph;

    const EarDecomposition ears = schmidt(biconnectedGraph);  // calculate proper ear decomposition

    // DEBUG
    // std::cerr << "ears\n" << ears.ears;

    const AdjacencyListGraph fromEars = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());

    // DBEUG
    // std::cerr << "fromEars\n" << fromEars;

    const AdjacencyListGraph minimal = fromEars.removeUncriticalEdges();

    // DBEUG
    // std::cerr << "minimal\n" << minimal;

    const EarDecomposition openEars = schmidt(minimal);

    // DEBUG
    std::cerr << "openEars\n" << openEars.ears;

    const std::vector<unsigned int> tour = hamiltonCycleInSquare(openEars, euclidean.numberOfNodes());

    std::cout << "hamiltoncycle\n" << tour;

    const Edge bottleneckEdge = findBottleneck(euclidean, tour, true);
    const double objective    = euclidean.weight(bottleneckEdge);

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
