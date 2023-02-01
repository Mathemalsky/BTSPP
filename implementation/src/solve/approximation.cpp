#include "solve/approximation.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"

namespace approximation {

Result approximate(const Euclidean& euclidean, const ProblemType problemType) {
  if (problemType == ProblemType::BTSP_approx) {
    AdjacencyMatrixGraph<Directionality::Undirected> graph = biconnectedSpanningGraph(euclidean);

    const OpenEarDecomposition ears = schmidt(graph);  // calculate proper ear decomposition

    // find approximation

    return Result{graph, ears};
  }
  else {
    throw std::runtime_error("Unknown problem type");
  }
}
}  // namespace approximation
