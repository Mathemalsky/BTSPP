#include "commandinterpreter.hpp"

#include <iostream>
#include <string>
#include <unordered_set>

#include "graph/graph.hpp"

#include "solve/approximation.hpp"
#include "solve/definitions.hpp"
#include "solve/euclideandistancegraph.hpp"
#include "solve/exactsolver.hpp"

#include "utility/utils.hpp"

void interpretCommandLine(const int argc, char* argv[]) {
  std::unordered_set<std::string> arguments;
  for (int i = 2; i < argc; ++i) {
    arguments.insert(std::string(argv[i]));
  }

  graph::Euclidean euclidean = generateEuclideanDistanceGraph(std::atoi(argv[1]));

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
