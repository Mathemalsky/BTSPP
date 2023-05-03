/*
 * GRAPH is a library to store and manipulate graphs as adjacency list or
 * as sparse eigen matrix. Different specialized types of graphs are
 * supported.
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

#include <ostream>

#include "graph.hpp"

namespace graph {
inline std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  return os << "(" << edge.u << ", " << edge.v << ") ";
}

template <typename G>
inline std::ostream& operator<<(std::ostream& os, const G& graph)
  requires(std::is_base_of_v<Graph, G>)
{
  os << "Number of nodes in graph: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph.edges()) {
    os << e << std::endl;
  }
  return os;
}
}  // namespace graph
