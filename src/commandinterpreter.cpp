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
/***********************************************************************************************************************
 *                                                      general
 **********************************************************************************************************************/

static bool findSeed(std::array<uint_fast32_t, SEED_LENGTH>& seed, char* argv[], int& i) {
  if (std::string(argv[i]) == "-seed") {
    for (unsigned int j = 0; j < SEED_LENGTH; ++j) {
      seed[j] = std::atoi(argv[++i]);
    }
    return true;
  }
  return false;
}

static void printSeedAdvice() {
  std::cout << "<-seed> <int1> ... <int" << SEED_LENGTH << "> to pass a seed for generating the graph\n";
}

static void printSyntax() {
  std::cout << "Syntax\n";
  std::cout << "./<NameOfExecutable> <NumberOfNodesInGraph> <arg1> <arg2> ...\n";
}

/***********************************************************************************************************************
 *                                                  visual program
 **********************************************************************************************************************/

#if (VISUALIZATION)

  #include "draw/visualization.hpp"
void interpretCommandLine(const int argc, char* argv[]) {
  if (argc == 2 && std::string(argv[1]) == "help") {
    printSyntax();
    printSeedAdvice();
    return;
  }
  if (argc != 2 && argc != 2 + 1 + SEED_LENGTH) {
    throw InvalidArgument("[COMMAND INTERPRETER] Got " + std::to_string(argc) + " arguments, but expected exactly 2 or additional seed!");
  }
  if (std::atoi(argv[1]) < 3) {
    throw InvalidArgument("[COMMAND INTERPRETER] Invalid argument, graph must have at least 3 vertices!");
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
      throw InvalidArgument("[COMMAND INTERPRETER] Invalid argument, expected a seed!");
    }
  }
  visualize(euclidean);
}
#endif
/***********************************************************************************************************************
 *                                               command line program
 **********************************************************************************************************************/

#if not(VISUALIZATION)
static void readArguments(const int argc, char* argv[]) {
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
    const approximation::Result res = approximation::approximateBTSP(euclidean);
    arguments.erase("-btsp");
  }
  if (arguments.contains("-btspp")) {
    const approximation::Result res = approximation::approximateBTSPP(euclidean);
    arguments.erase("-btspp");
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

static void printArgumentList() {
  std::cout << "Valid command line arguments are: \n";
  std::cout << "<-btsp> to approximate BTSP\n";
  std::cout << "<-btspp> to approximate BTSP\n";
  std::cout << "<-btsp-e> to solve exact BTSP\n";
  std::cout << "<-btspp-e> to solve exact BTSPP\n";
  std::cout << "<-tsp-e> to solve exact TSP\n";
}

void interpretCommandLine(const int argc, char* argv[]) {
  if (argc == 2 && std::string(argv[1]) == "help") {
    printSyntax();
    printArgumentList();
    printSeedAdvice();
    return;
  }
  if (argc < 3) {
    throw InvalidArgument("[COMMAND INTERPRETER] To few arguments! Got " + std::to_string(argc) + " but expected at least 3!");
  }
  if (std::atoi(argv[1]) < 3) {
    throw InvalidArgument("[COMMAND INTERPRETER] Invalid argument, graph must have at least 3 vertices!");
  }
  readArguments(argc, argv);
}
#endif