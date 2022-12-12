#pragma once

#include <vector>

#include "graph/graph.hpp"

/*!
 * \brief solveExact solve the tsp to optimality using Miller-Tucker-Zemlin formulation
 * \param euclidean
 * \details x_ij belongs to the column j* numberOfNodes +i
 * hereby we ommit th x_jj entries and store u_j in that place
 */
std::vector<size_t> solveTSP(const Euclidean& euclidean);
