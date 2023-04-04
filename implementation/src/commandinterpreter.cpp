#include "commandinterpreter.hpp"

#include <array>
#include <iostream>
#include <string>
#include <unordered_set>

#include "graph/graph.hpp"

#include "solve/approximation.hpp"
#include "solve/definitions.hpp"
#include "solve/euclideandistancegraph.hpp"
#include "solve/exactsolver.hpp"

#include "utility/utils.hpp"

bool findSeed(std::array<uint_fast32_t, SEED_LENGTH>& seed, char* argv[], int& i) {
  if (std::string(argv[i]) == "-seed") {
    for (unsigned int j = 0; j < SEED_LENGTH; ++j) {
      seed[j] = std::atoi(argv[++i]);
    }
    return true;
  }
  return false;
}

void interpretCommandLine(const int argc, char* argv[]) {
  std::unordered_set<std::string> arguments;
  bool seeded = false;
  std::array<uint_fast32_t, SEED_LENGTH> seed;
  for (int i = 2; i < argc; ++i) {
    if (findSeed(seed, argv, i)) {
      seeded = true;
      continue;
    }
    arguments.insert(std::string(argv[i]));
  }

  graph::Euclidean euclidean;
  if (seeded) {
    euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]), seed);
  }
  else {
    euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]));
  }

  if (arguments.contains("-btsp")) {
    const approximation::Result res = approximation::approximate(euclidean, ProblemType::BTSP_approx);
    arguments.erase("-btsp");
  }
  if (arguments.contains("-btsp-e")) {
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::BTSP_exact);
    arguments.erase("-btsp-e");
  }
  if (arguments.contains("-btspp-e")) {
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::BTSPP_exact);
    arguments.erase("-btspp-e");
  }
  if (arguments.contains("-tsp-e")) {
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::TSP_exact);
    arguments.erase("-tsp-e");
  }

  for (const std::string& str : arguments) {
    printYellow("Warning");
    std::cout << ": Unknown argument <" << str << ">!" << std::endl;
  }
}
