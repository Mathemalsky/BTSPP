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
#pragma once

#include <vector>

#include <solve/definitions.hpp>

// graph library
#include "graph.hpp"

namespace exactsolver {

struct Result {
  std::vector<unsigned int> tour;
  double opt;
  graph::Edge bottleneckEdge;
};

/*!
 * @brief solves an instance of BTSP, BTSPP or TSP
 * @param problemType type of instance
 * @param noCrossing if BTSP the solution can be forced to have no crossings
 */
Result solve(const graph::Euclidean& euclidean, const ProblemType problemType, const bool noCrossing = false);

}  // namespace exactsolver
