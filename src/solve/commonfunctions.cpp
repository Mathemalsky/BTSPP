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
#include "solve/commonfunctions.hpp"

#include <vector>

// graph library
#include "graph.hpp"

graph::Edge findBottleneck(const graph::Euclidean& euclidean, const std::vector<unsigned int>& tour, const bool isCycle) {
  unsigned int bottleneckEdgeEnd = 0;
  double bottleneckWeight        = euclidean.weight(tour[0], tour[1]);
  for (unsigned int i = 1; i < euclidean.numberOfNodes() - 1; ++i) {
    if (euclidean.weight(tour[i], tour[i + 1]) > bottleneckWeight) {
      bottleneckEdgeEnd = i;
      bottleneckWeight  = euclidean.weight(tour[i], tour[i + 1]);
    }
  }
  if (isCycle && euclidean.weight(tour.back(), tour[0]) > bottleneckWeight) {
    return graph::Edge{tour.back(), 0};
  }
  else {
    return graph::Edge{tour[bottleneckEdgeEnd], tour[bottleneckEdgeEnd + 1]};
  }
}
