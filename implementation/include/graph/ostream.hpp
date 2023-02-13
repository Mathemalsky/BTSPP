#pragma once

#include <ostream>

#include "graph/graph.hpp"

inline std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  return os << "(" << edge.u << ", " << edge.v << ") ";
}

inline std::ostream& operator<<(std::ostream& os, const AdjMatGraph& graph) {
  os << "Number of nodes: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph.edges()) {
    os << e << " " << graph.matrix().coeff(e.u, e.v).cost() << std::endl;
  }
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const DfsTree& tree) {
  os << "Number of nodes: " << tree.numberOfNodes() << std::endl;
  for (const Edge& e : tree) {
    os << e << std::endl;
  }
  return os;
}

template <Directionality DIRECT>
inline std::ostream& operator<<(std::ostream& os, const AdjacencyListGraph<DIRECT>& graph) {
  os << "Number of nodes: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph) {
    os << e << std::endl;
  }
  return os;
}
