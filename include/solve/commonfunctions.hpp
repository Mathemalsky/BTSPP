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

#include <fstream>
#include <vector>

// graph library
#include "graph.hpp"

#include "solve/definitions.hpp"

template <typename G>
  requires(std::is_base_of_v<graph::CompleteGraph, G> && std::is_base_of_v<graph::WeightedGraph, G>)
graph::Edge findBottleneck(const G& completeGraph, const std::vector<size_t>& tour, const bool isCycle) {
  size_t bottleneckEdgeEnd = 0;
  double bottleneckWeight  = completeGraph.weight(tour[0], tour[1]);
  for (size_t i = 1; i < completeGraph.numberOfNodes() - 1; ++i) {
    if (completeGraph.weight(tour[i], tour[i + 1]) > bottleneckWeight) {
      bottleneckEdgeEnd = i;
      bottleneckWeight  = completeGraph.weight(tour[i], tour[i + 1]);
    }
  }
  if (isCycle && completeGraph.weight(tour.back(), tour[0]) > bottleneckWeight) {
    return graph::Edge{tour.back(), 0};
  }
  else {
    return graph::Edge{tour[bottleneckEdgeEnd], tour[bottleneckEdgeEnd + 1]};
  }
}

inline std::ostream& operator<<(std::ostream& os, const ProblemType type) {
  if (type == ProblemType::BTSP_approx || type == ProblemType::BTSP_exact) {
    os << "BTSP";
  }
  else if (type == ProblemType::BTSPP_approx || type == ProblemType::BTSPP_exact) {
    os << "BTSPP";
  }
  else if (type == ProblemType::TSP_exact) {
    os << "TSP";
  }
  return os;
}
