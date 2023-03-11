#include <cstdlib>

#include "draw/visualization.hpp"

#include "graph/graph.hpp"

#include "solve/euclideandistancegraph.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  graph::Euclidean euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]));

#if (VISUALIZATION)
  return visualize(euclidean);
#else
#endif
}
