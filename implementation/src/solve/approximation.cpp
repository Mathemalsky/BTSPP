#include "solve/approximation.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"

// DEBUG
#include <iostream>
#include "utility/utils.hpp"
#include "graph/ostream.hpp"

namespace approximation {

Result approximate(const Euclidean& euclidean, const ProblemType problemType) {
  if (problemType == ProblemType::BTSP_approx) {
    AdjacencyMatrixGraph<Directionality::Undirected> graph = biconnectedSpanningGraph(euclidean);

    const EarDecomposition ears = schmidt(graph);  // calculate proper ear decomposition

    // DEBUG
    std::cerr << "ears\n" << ears.ears;

    // find approximation
    const AdjacencyMatrixGraph<Directionality::Undirected> fromEars = earDecompToGraph(ears);

    // DBEUG
    std::cerr << "fromEars\n" << fromEars;

    /*
    // TEST
    std::cerr << "neighbours of 0\n";
    for (const size_t u : fromEars.neighbours(0)) {
      std::cerr << u << " ";
    }
    std::cerr << std::endl;

    std::cerr << "neighbours of 3\n";
    for (const size_t u : fromEars.neighbours(3)) {
      std::cerr << u << " ";
    }
    std::cerr << std::endl;
    */

    const AdjacencyMatrixGraph<Directionality::Undirected> minimal = fromEars.removeUncriticalEdges();

    // DBEUG
    std::cerr << "minimal\n" << minimal;

    const EarDecomposition openEars = schmidt(minimal);

    // DEBUG
    std::cerr << "openEars\n" << openEars.ears;

    return Result{graph, openEars};
  }
  else {
    throw std::runtime_error("Unknown problem type");
  }
}
}  // namespace approximation
