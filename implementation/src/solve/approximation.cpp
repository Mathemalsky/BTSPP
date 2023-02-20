#include "solve/approximation.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"
#include "graph/algorithm.hpp"

// DEBUG
#include <iostream>
#include "utility/utils.hpp"
#include "graph/ostream.hpp"

namespace approximation {

static std::vector<size_t> shortcutToHamiltoncycle(
    const std::vector<size_t>& eulertourInEarDecomposition, const size_t numberOfNodes) {
  std::vector<bool> visited(numberOfNodes, false);
  std::vector<size_t> hamiltoncycle;
  hamiltoncycle.reserve(numberOfNodes);

  for (const size_t u : eulertourInEarDecomposition) {
    if (!visited[u]) {
      hamiltoncycle.push_back(static_cast<unsigned int>(u));
    }
  }
  return hamiltoncycle;
}

static std::vector<unsigned int> backAndForth(const std::vector<size_t>& vecIn) {
  std::vector<unsigned int> vecOut;
  vecOut.reserve(vecIn.size() - 1);

  for (size_t i = 1; i <= (vecIn.size() - 1) / 2; ++i) {
    vecOut.push_back(vecIn[2 * i]);  // even positions in ascending order
  }
  for (size_t i = vecIn.size() / 2; i > 0; --i) {
    vecOut.push_back(vecIn[2 * i - 1]);  // odd positions in descending order
  }

  return vecOut;
}

static std::vector<unsigned int> assembleTour(
    const std::vector<size_t>& hamiltonSubcycle, const std::vector<std::vector<size_t>>& nodesPreceedingDoubleEdges,
    const size_t numberOfNodes) {
  std::unordered_set<size_t> doubleEdgeStarts;
  for (const std::vector<size_t>& doubleEdgeSegment : nodesPreceedingDoubleEdges) {
    doubleEdgeStarts.insert(doubleEdgeSegment[0]);
  }

  std::vector<unsigned int> hamiltoncycle;
  hamiltoncycle.reserve(numberOfNodes);

  for (const size_t u : hamiltonSubcycle) {
    hamiltoncycle.push_back(u);
    if (doubleEdgeStarts.contains(u)) {
      doubleEdgeStarts.erase(u);
      const std::vector<unsigned int> tmp = backAndForth(nodesPreceedingDoubleEdges[u]);
      hamiltoncycle.insert(hamiltoncycle.end(), tmp.begin(), tmp.end());
    }
  }
  return hamiltoncycle;
}

static std::vector<unsigned int> hamiltonCycleInSquare(const EarDecomposition& ears, const size_t numberOfNodes) {
  if (ears.ears.size() == 1) {
    return std::vector<unsigned int>(ears.ears[0].begin(), ears.ears[0].end());  // cast down to unsigned int
  }
  else {
    AdjacencyListGraph graph = earDecompToAdjacencyListGraph(ears, numberOfNodes);
    std::vector<std::vector<size_t>> nodesPreceedingDoubleEdges;

    // make degrees even
    for (long i = ears.ears.size() - 2; i >= 0; --i) {
      const std::vector<size_t>& ear = ears.ears[i];
      bool removedPreviousEdge       = false;

      for (size_t i = 1; i < ear.size() - 1; ++i) {
        const size_t u = ear[i];
        if (graph.degree(u) % 2 == 1) {
          graph.removeEdge(Edge{u, ear[i + 1]});
          if (!removedPreviousEdge) {
            nodesPreceedingDoubleEdges.push_back(std::vector<size_t>{u});
          }
          else {
            nodesPreceedingDoubleEdges.back().push_back(u);
          }
          removedPreviousEdge = true;
        }
        else {
          removedPreviousEdge = false;
        }
      }
    }

    // find euler tour
    std::vector<size_t> eulerSubtour     = eulertour(graph);
    std::vector<size_t> hamiltonSubcycle = shortcutToHamiltoncycle(eulerSubtour, numberOfNodes);
    return assembleTour(hamiltonSubcycle, nodesPreceedingDoubleEdges, numberOfNodes);
  }
}

Result approximate(const Euclidean& euclidean, const ProblemType problemType) {
  if (problemType == ProblemType::BTSP_approx) {
    AdjacencyMatrixGraph biconnectedGraph = biconnectedSpanningGraph(euclidean);

    // DEBUG
    std::cerr << "undirected\n" << biconnectedGraph;

    const EarDecomposition ears = schmidt(biconnectedGraph);  // calculate proper ear decomposition

    // DEBUG
    std::cerr << "ears\n" << ears.ears;

    const AdjacencyListGraph fromEars = earDecompToAdjacencyListGraph(ears, biconnectedGraph.numberOfNodes());

    // DBEUG
    std::cerr << "fromEars\n" << fromEars;

    const AdjacencyListGraph minimal = fromEars.removeUncriticalEdges();

    // DBEUG
    std::cerr << "minimal\n" << minimal;

    const EarDecomposition openEars = schmidt(minimal);

    // DEBUG
    std::cerr << "openEars\n" << openEars.ears;

    // find hamiltonian cycle
    const std::vector<unsigned int> tour = hamiltonCycleInSquare(openEars, euclidean.numberOfNodes());

    return Result{biconnectedGraph, openEars, tour};
  }
  else {
    throw std::runtime_error("Unknown problem type");
  }
}
}  // namespace approximation
