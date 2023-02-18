#include "solve/approximation.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"
#include "graph/algorithm.hpp"

// DEBUG
#include <iostream>
#include "utility/utils.hpp"
#include "graph/ostream.hpp"

namespace approximation {

template <typename G>
size_t findDegree2Node(const G& graph, const std::vector<size_t>& ear) {
  for (size_t i = 1; i < ear.size() - 1; ++i) {
    const size_t u = ear[i];
    if (graph.degree(u) == 2) {
      return u;
    }
  }
  assert(false && "There should be a degree 2 node in the ear!");
}

std::vector<unsigned int> findTour(const EarDecomposition& ears, const size_t numberOfNodes) {
  if (ears.ears.size() == 1) {
    return std::vector<unsigned int>(ears.ears[0].begin(), ears.ears[0].end());  // cast down to unsigned int
  }
  else {
    AdjacencyListGraph graph    = earDecompToAdjacencyListGraph(ears, numberOfNodes);
    std::vector<size_t> degrees = graph.degrees();

    // make degrees even
    for (long i = ears.ears.size() - 2; i >= 0; --i) {
      const std::vector<size_t>& ear = ears.ears[i];
      const size_t y                 = findDegree2Node(graph, ear);
      bool lastEdgeRemoved           = false;
      for (size_t i = 1; i < ear.size() - 1; ++i) {
        const size_t u = ear[i];
        if (graph.degree(u) % 2 == 1) {
          graph.removeEdge(Edge{u, ear[i + 1]});
          if (i == ear.size() - 1) {
            lastEdgeRemoved = true;
          }
          else {
            // store the deleted edge elsewhere
          }
        }
      }
    }

    // find euler tour
    // shortcut to hamiltonian cycle

    // DEBUG dummy return value
    std::vector<unsigned int> tour(numberOfNodes);

    return tour;
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
    const std::vector<unsigned int> tour = findTour(openEars, euclidean.numberOfNodes());

    return Result{biconnectedGraph, openEars, tour};
  }
  else {
    throw std::runtime_error("Unknown problem type");
  }
}
}  // namespace approximation
