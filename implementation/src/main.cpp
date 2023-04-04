#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "draw/visualization.hpp"

#include "graph/graph.hpp"

#include "solve/euclideandistancegraph.hpp"

#include "utility/utils.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
#if (VISUALIZATION)
  try {
    if (argc != 2) {
      throw std::invalid_argument("[MAIN] Got " + std::to_string(argc) + " arguments, but expected exactly 2!");
    }
    if (std::atoi(argv[1]) < 3) {
      throw std::invalid_argument("[MAIN] Graph must have at least 3 vertices!");
    }
    graph::Euclidean euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]));
    visualize(euclidean);
  }
  catch (std::exception& error) {
    printLightred("Error");
    std::cerr << ": " << error.what() << std::endl;
    return -1;
  }

#else
  try {
    if (argc < 3) {
      throw std::invalid_argument("To few arguments! Got " + std::to_string(argc) + " but expected at least 3!");
    }
    if (std::atoi(argv[1]) < 3) {
      throw std::invalid_argument("[MAIN] Graph must have at least 3 vertices!");
    }
    graph::Euclidean euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]));
  }
  catch (...) {
  }

#endif

  return 0;
}
