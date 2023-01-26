#pragma once

#include <vector>

#include <draw/definitions.hpp>

#include "graph/graph.hpp"

struct Result {
  std::vector<unsigned int> tour;
  Edge bottleneck;
};

/*!
 * \brief solveExact solve the tsp to optimality using Miller-Tucker-Zemlin formulation
 * \param euclidean
 * \details x_ij belongs to the column j* numberOfNodes +i
 * hereby we ommit th x_jj entries and store u_j in that place
 */
std::vector<unsigned int> solve(const Euclidean& euclidean, const ProblemType problemType);
