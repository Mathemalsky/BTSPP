#pragma once

#include <ostream>

#include "graph/graph.hpp"

namespace graph {
inline std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  return os << "(" << edge.u << ", " << edge.v << ") ";
}

template <typename G>
inline std::ostream& operator<<(std::ostream& os, const G& graph) requires(std::is_base_of_v<Graph, G>) {
  os << "Number of nodes in graph: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph.edges()) {
    os << e << std::endl;
  }
  return os;
}
}  // namespace graph
