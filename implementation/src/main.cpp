#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "draw/visualization.hpp"

#include "graph/graph.hpp"

#include "solve/definitions.hpp"
#include "solve/euclideandistancegraph.hpp"

#include "utility/utils.hpp"

#include "commandinterpreter.hpp"
#include "exceptions.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
#if (VISUALIZATION)
  try {
    if (argc < 2 || argc != 2 + 1 + SEED_LENGTH) {
      throw InvalidArgument("[MAIN] Got " + std::to_string(argc) + " arguments, but expected exactly 2 or additional seed!");
    }
    if (std::atoi(argv[1]) < 3) {
      throw InvalidArgument("[MAIN] Graph must have at least 3 vertices!");
    }

    graph::Euclidean euclidean;
    if (argc == 2) {
      euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]));
    }
    else {
      std::array<uint_fast32_t, SEED_LENGTH> seed;
      int i = 2;
      if (findSeed(seed, argv, i)) {
        euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]), seed);
      }
      else {
        throw InvalidArgument("Expected a seed!");
      }
    }
    visualize(euclidean);
  }
  catch (Exception& error) {
    printLightred("Error");
    std::cerr << ": " << error.what() << std::endl;
    return -1;
  }

#else
  try {
    if (argc < 3) {
      throw InvalidArgument("To few arguments! Got " + std::to_string(argc) + " but expected at least 3!");
    }
    if (std::atoi(argv[1]) < 3) {
      throw InvalidArgument("[MAIN] Graph must have at least 3 vertices!");
    }
    interpretCommandLine(argc, argv);
  }
  catch (Exception& error) {
    printLightred("Error");
    std::cerr << ": " << error.what() << std::endl;
    return -1;
  }

#endif

  return 0;
}
