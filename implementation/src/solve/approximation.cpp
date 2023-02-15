#include "solve/approximation.hpp"

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

Result approximate(const Euclidean& euclidean, const ProblemType problemType) {
  if (problemType == ProblemType::BTSP_approx) {
    AdjacencyMatrixGraph<Directionality::Undirected> graph = biconnectedSpanningGraph(euclidean);

    // DEBUG
    std::cerr << "undirected\n";
    for (const Edge& e : graph.edgesToLowerIndex()) {
      std::cerr << e << std::endl;
    }

    const EarDecomposition ears = schmidt(graph);  // calculate proper ear decomposition

    // DEBUG
    std::cerr << "ears\n" << ears.ears;

    const AdjacencyListGraph<Directionality::Undirected> fromEars =
        earDecompToAdjacencyListGraph(ears, graph.numberOfNodes());

    // DBEUG
    std::cerr << "fromEars\n" << fromEars;

    const AdjacencyListGraph<Directionality::Undirected> minimal = fromEars.removeUncriticalEdges();

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
