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

static void printSyntax() {
  std::cout << "Syntax\n";
  std::cout << "./<NameOfExecutable> <NumberOfNodesInGraph> <arg1> <arg2> ...\n";
}

/***********************************************************************************************************************
 *                                                  visual program
 **********************************************************************************************************************/

#if (VISUALISATION)

  #include "draw/visualisation.hpp"
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
  visualise(euclidean);
}
#endif
/***********************************************************************************************************************
 *                                               command line program
 **********************************************************************************************************************/

#if not(VISUALISATION)
constexpr std::string_view LOG_FILE_IDENTIFIER   = "-logfile:=";
constexpr std::string_view REPETITION_IDENTIFIER = "-repetitions:=";
constexpr std::string_view SUPPRESS_INFO_TAG     = "-no-info";
constexpr std::string_view SUPPRESS_SEED_TAG     = "-no-seed";
constexpr std::string_view NO_CROSSING_TAG       = "-no-crossing";
constexpr std::string_view BTSP_APPROX_TAG       = "-btsp";
constexpr std::string_view BTSPP_APPROX_TAG      = "-btspp";
constexpr std::string_view BTSP_EXACT_TAG        = "-btsp-e";
constexpr std::string_view BTSPP_EXACT_TAG       = "-btspp-e";
constexpr std::string_view TSP_EXACT_TAG         = "-tsp-e";

static void printAdvices() {
  std::cout << "<" << NO_CROSSING_TAG << "> if <-btsp-e> is set, to find a solution without crossing.\n";
  std::cout << "<" << LOG_FILE_IDENTIFIER << "<filename>> to write infos to <filename>.\n";
  std::cout << "<" << REPETITION_IDENTIFIER << "<numberOfRepetitions>> to compute several instances serial in one execution.\n";
  std::cout << "<" << SUPPRESS_INFO_TAG << " to suppress detailed terminal output.\n";
  std::cout << "<" << SUPPRESS_SEED_TAG << " to suppress seed output in terminal.\n";
}
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
                               const double runtime,
                               const bool suppressInfo) {
  if (!suppressInfo) {
    printInfo(res, type, runtime);
  }
  if (filename.length() > 0) {
    writeStatsToFile(res, type, filename, runtime);
  }
}

static graph::Euclidean adaptSeededGeneration(const size_t numberOfNodes,
                                              const std::array<uint_fast32_t, SEED_LENGTH>& seed,
                                              const bool seeded,
                                              bool suppressSeed) {
  if (seeded) {
    return generateEuclideanDistanceGraph(numberOfNodes, seed, suppressSeed);
  }
  else {
    return generateEuclideanDistanceGraph(numberOfNodes, suppressSeed);
  }
}

static void readArguments(const int argc, char* argv[]) {
  std::string filename = "";
  size_t repetitions   = 1;
  bool suppressInfo = false, suppressSeed = false;
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
    if (std::string(argv[i]) == SUPPRESS_INFO_TAG) {
      suppressInfo = true;
      continue;
    }
    if (std::string(argv[i]) == SUPPRESS_SEED_TAG) {
      suppressSeed = true;
      continue;
    }
    arguments.insert(std::string(argv[i]));
  }

  if (arguments.empty()) {
    printYellow("Warning");
    std::cout << ": No problem type given. Nothing to do." << std::endl;
  }

  Stopwatch stopWatch;  // create stop watch

  if (arguments.contains(std::string(BTSP_APPROX_TAG))) {
    for (size_t i = 0; i < repetitions; ++i) {
      graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded, suppressSeed);
      stopWatch.reset();
      const approximation::Result res = approximation::approximateBTSP(euclidean);
      const double runtime            = stopWatch.elapsedTimeInMilliseconds();
      handleApproxOutput(res, ProblemType::BTSP_approx, filename, runtime, suppressInfo);
    }
    arguments.erase(std::string(BTSP_APPROX_TAG));
  }
  if (arguments.contains(std::string(BTSPP_APPROX_TAG))) {
    for (size_t i = 0; i < repetitions; ++i) {
      graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded, suppressSeed);
      stopWatch.reset();
      const approximation::Result res = approximation::approximateBTSPP(euclidean);
      const double runtime            = stopWatch.elapsedTimeInMilliseconds();
      handleApproxOutput(res, ProblemType::BTSPP_approx, filename, runtime, suppressInfo);
    }
    arguments.erase(std::string(BTSPP_APPROX_TAG));
  }
  if (arguments.contains(std::string(BTSP_EXACT_TAG))) {
    graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded, suppressSeed);
    stopWatch.reset();
    const exactsolver::Result res =
        exactsolver::solve(euclidean, ProblemType::BTSP_exact, arguments.contains(std::string(NO_CROSSING_TAG)));
    const double runtime = stopWatch.elapsedTimeInMilliseconds();
    printInfo(res, ProblemType::BTSP_exact, runtime);
    arguments.erase(std::string(BTSP_EXACT_TAG));
    arguments.erase(std::string(NO_CROSSING_TAG));
  }
  if (arguments.contains(std::string(BTSPP_EXACT_TAG))) {
    graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded, suppressSeed);
    stopWatch.reset();
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::BTSPP_exact);
    const double runtime          = stopWatch.elapsedTimeInMilliseconds();
    printInfo(res, ProblemType::BTSPP_exact, runtime);
    arguments.erase(std::string(BTSPP_EXACT_TAG));
  }
  if (arguments.contains(std::string(TSP_EXACT_TAG))) {
    graph::Euclidean euclidean = adaptSeededGeneration(std::atoi(argv[1]), seed, seeded, suppressSeed);
    stopWatch.reset();
    const exactsolver::Result res = exactsolver::solve(euclidean, ProblemType::TSP_exact);
    const double runtime          = stopWatch.elapsedTimeInMilliseconds();
    printInfo(res, ProblemType::TSP_exact, runtime);
    arguments.erase(std::string(TSP_EXACT_TAG));
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
  std::cout << "<" << BTSP_APPROX_TAG << "> to approximate BTSP\n";
  std::cout << "<" << BTSPP_APPROX_TAG << "> to approximate BTSP\n";
  std::cout << "<" << BTSP_EXACT_TAG << "> to solve exact BTSP\n";
  std::cout << "<" << BTSPP_EXACT_TAG << "> to solve exact BTSPP\n";
  std::cout << "<" << TSP_EXACT_TAG << "> to solve exact TSP\n";
}

void interpretCommandLine(const int argc, char* argv[]) {
  if (argc == 2 && std::string(argv[1]) == "help") {
    printSyntax();
    printArgumentList();
    printSeedAdvice();
    printAdvices();
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