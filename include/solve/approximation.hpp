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

// graph library
#include "algorithm.hpp"
#include "graph.hpp"

#include "definitions.hpp"

namespace approximation {
struct Result {
  graph::AdjacencyListGraph biconnectedGraph;
  graph::EarDecomposition openEarDecomposition;
  std::vector<unsigned int> tour;
  double opt;
  graph::Edge bottleneckEdge;
};

/*!
 * \brief approximates the BTSP
 * \param euclidean complete graph, providing distances between nodes
 * \param printInfo controls if objective, lower bound on OPT and a fortiori guarantee are printed to console
 * \return Result
 */
Result approximateBTSP(const graph::Euclidean& euclidean, const bool printInfo = true);

/*!
 * @brief approximates the BTSPP
 * @param euclidean complete graph, providing distances between nodes
 * @param s start node
 * @param t end node
 * @param printInfo controls if objective, lower bound on OPT and a fortiori guarantee are printed to console
 * @return Result
 */
Result approximateBTSPP(const graph::Euclidean& euclidean, const size_t s = 0, const size_t t = 1, const bool printInfo = true);
}  // namespace approximation
