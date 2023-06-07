/*
 * pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
 * BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "commandinterpreter.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>

// graph library
#include "graph.hpp"

#include "exception/exceptions.hpp"

#include "solve/approximation.hpp"
#include "solve/definitions.hpp"
#include "solve/euclideandistancegraph.hpp"
#include "solve/exactsolver.hpp"

#include "utility/utils.hpp"
/***********************************************************************************************************************
 *                                                      general
 **********************************************************************************************************************/

/*!
 * @brief extracts seed from list of arguments
 * @details if the keyword <-seed> is found, the following SEED_LENGTH elements from the argument list are read as seed,
 * the position index is moved forward by SEED_LENGTH positions
 * @param seed array to write the seed to
 * @param argv argument list as passed to main
 * @param i position index in argument list
 * @return true if a seed was found
 */
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

static void printNoCrossingAdvice() {
  std::cout << "<-no-crossing> if <-btsp-e> is set, to find a solution without crossing\n";
}

[[maybe_unused]] static void printLogfileAdvice() {
  std::cout << "<-logfile:=<filename>> to write infos to <filename>\n";
}

[[maybe_unused]] static void printRepetitionAdvice() {
  std::cout << "<-repetitions:=<numberOfRepetitions>> to compute several instances serial in one execution.\n";
}

static void printSyntax() {
  std::cout << "Syntax\n";
  std::cout << "./<NameOfExecutable> <NumberOfNodesInGraph> <arg1> <arg2> ...\n";
}

[[maybe_unused]] static void printInfo(const approximation::Result& res, const ProblemType problemType, const double runtime) {
  std::cout << "-------------------------------------------------------\n";
  std::cout << "Approximated an instance of " << INSTANCE_TYPES[std::to_underlying(problemType)] << std::endl;
  std::cout << "objective                            : " << res.objective << std::endl;
  std::cout << "lower bound on OPT                   : " << res.lowerBoundOnOPT << std::endl;
  std::cout << "a fortiori guarantee                 : " << res.objective / res.lowerBoundOnOPT << std::endl;
  std::cout << "edges in biconnected graph           : " << res.biconnectedGraph.numberOfEdges() << std::endl;
  std::cout << "edges in minimally biconnected graph : " << res.numberOfEdgesInMinimallyBiconectedGraph << std::endl;
  std::cout << "elapsed time                         : " << runtime << " ms\n";
}

[[maybe_unused]] static void printInfo(const exactsolver::Result& res, const ProblemType problemType, const double runtime) {
  std::cout << "-------------------------------------------------------\n";
  std::cout << "Solved an instance of " << INSTANCE_TYPES[std::to_underlying(problemType)] << std::endl;
  std::cout << "OPT                                  : " << res.opt << std::endl;
  std::cout << "elapsed time                         : " << runtime << " ms\n";
}

/***********************************************************************************************************************
 *                                                  visual program
 **********************************************************************************************************************/

#if (VISUALISATION)

  #include "draw/visualization.hpp"
void interpretCommandLine(const int argc, char* argv[]) {
  if (argc == 2 && std::string(argv[1]) == "help") {
    printSyntax();
    printSeedAdvice();
    printNoCrossingAdvice();
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

#if not(VISUALISATION)
constexpr std::string_view LOG_FILE_IDENTIFIER   = "-logfile:=";
constexpr std::string_view REPETITION_IDENTIFIER = "-repetitions:=";

static void writeStatsToFile(const approximation::Result& res, const ProblemType type, const std::string& filename, const double runtime) {
  std::ofstream outputfile;
  outputfile.open(filename, std::ios::out | std::ios::app);
  if (!outputfile) {
    throw InvalidFileOperation("Failed to open <" + filename + ">!");
  }
  outputfile << std::to_underlying(type) << ",";
  outputfile << res.biconnectedGraph.numberOfNodes() << ",";
  outputfile << res.objective << ",";
  outputfile << res.lowerBoundOnOPT << ",";
  outputfile << res.objective / res.lowerBoundOnOPT << ",";
  outputfile << res.biconnectedGraph.numberOfEdges() << ",";
  outputfile << res.numberOfEdgesInMinimallyBiconectedGraph << ",";
  outputfile << runtime << std::endl;
}

static void handleApproxOutput(const approximation::Result& res,
                               const ProblemType type,
                               const std::string& filename,
                               const double runtime) {
  printInfo(res, type, runtime);
  if (filename.length() > 0) {
    writeStatsToFile(res, type, filename, runtime);
  }
}

static graph::Euclidean adaptSeededGeneration(const size_t numberOfNodes,
                                              const std::array<uint_fast32_t, SEED_LENGTH>& seed,
                                              const bool seeded) {
  if (seeded) {
    return generateEuclideanDistanceGraph(numberOfNodes, seed);
  }
  else {
    return generateEuclideanDistanceGraph(numberOfNodes);
  }
}

static void readArguments(const int argc, char* argv[]) {
  std::string filename = "";
  size_t repetitions   = 1;
  std::unordered_set<std::string> arguments;
  bool seeded = false;
  std::array<uint_fast32_t, SEED_LENGTH> seed;
  for (int i = 2; i < argc; ++i) {
    if (findSeed(seed, argv, i)) {
      seeded = true;
      continue;
    }
    if (std::string(argv[i]).starts_with(LOG_FILE_IDENTIFIER)) {
      filename = std::string(argv[i]).substr(LOG_FILE_IDENTIFIER.length(), std::string(argv[i]).length() - LOG_FILE_IDENTIFIER.length());
      continue;
    }
    if (std::string(argv[i]).starts_with(REPETITION_IDENTIFIER)) {
      repetitions = std::stoul(
          std::string(argv[i]).substr(REPETITION_IDENTIFIER.length(), std::string(argv[i]).length() - REPETITION_IDENTIFIER.length()));
      continue;
    }
    arguments.insert(std::string(argv[i]));
  }

  if (arguments.empty()) {
    printYellow("Warning");
    std::cout << ": No problem type given. Nothing to do." << std::endl;
  }

  Stopwatch stopWatch;  // create stop watch

  if (arguments.contains("-btsp")) {
    for (size_t i = 0; i < repetitions; ++i) {
      graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded);
      stopWatch.reset();
      const approximation::Result res = approximation::approximateBTSP(euclidean);
      const double runtime            = stopWatch.elapsedTimeInMilliseconds();
      handleApproxOutput(res, ProblemType::BTSP_approx, filename, runtime);
    }
    arguments.erase("-btsp");
  }
  if (arguments.contains("-btspp")) {
    for (size_t i = 0; i < repetitions; ++i) {
      graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded);
      stopWatch.reset();
      const approximation::Result res = approximation::approximateBTSPP(euclidean);
      const double runtime            = stopWatch.elapsedTimeInMilliseconds();
      handleApproxOutput(res, ProblemType::BTSPP_approx, filename, runtime);
    }
    arguments.erase("-btspp");
  }
  if (arguments.contains("-btsp-e")) {
    graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded);
    stopWatch.reset();
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::BTSP_exact, arguments.contains("-no-crossing"));
    const double runtime          = stopWatch.elapsedTimeInMilliseconds();
    printInfo(res, ProblemType::BTSP_exact, runtime);
    arguments.erase("-btsp-e");
    arguments.erase("-no-crossing");
  }
  if (arguments.contains("-btspp-e")) {
    graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded);
    stopWatch.reset();
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::BTSPP_exact);
    const double runtime          = stopWatch.elapsedTimeInMilliseconds();
    printInfo(res, ProblemType::BTSPP_exact, runtime);
    arguments.erase("-btspp-e");
  }
  if (arguments.contains("-tsp-e")) {
    graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded);
    stopWatch.reset();
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::TSP_exact);
    const double runtime          = stopWatch.elapsedTimeInMilliseconds();
    printInfo(res, ProblemType::TSP_exact, runtime);
    arguments.erase("-tsp-e");
  }

  if (filename.length() > 0) {
    printLightgreen("Info");
    std::cout << ": Output has been written to <" << filename << ">.\n";
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
    printNoCrossingAdvice();
    printLogfileAdvice();
    printRepetitionAdvice();
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