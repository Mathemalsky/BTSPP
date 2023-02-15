#pragma once

#include <ostream>

#include "graph/graph.hpp"

inline std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  return os << "(" << edge.u << ", " << edge.v << ") ";
}

inline std::ostream& operator<<(std::ostream& os, const AdjMatGraph& graph) {
  os << "Number of nodes matrix graph: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph.edges()) {
    os << e << " " << graph.matrix().coeff(e.u, e.v).cost() << std::endl;
  }
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const DfsTree& tree) {
  os << "Number of nodes in tree: " << tree.numberOfNodes() << std::endl;
  for (const Edge& e : tree.edges()) {
    os << e << std::endl;
  }
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const AdjListGraph& graph) {
  os << "Number of nodes in list graph: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph.edgesToLowerIndex()) {
    os << e << std::endl;
  }
  return os;
}
