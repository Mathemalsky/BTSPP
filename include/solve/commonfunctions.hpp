#pragma once

#include <fstream>
#include <vector>

#include "graph/graph.hpp"

#include "solve/definitions.hpp"

graph::Edge findBottleneck(const graph::Euclidean& euclidean, const std::vector<unsigned int>& tour, const bool cycle);

template <typename Type>
Type previousInCycle(const std::vector<Type>& vec, const size_t position) {
  return (position != 0 ? vec[position - 1] : vec.back());
}

template <typename Type>
Type successiveInCycle(const std::vector<Type>& vec, const size_t position) {
  return (position != vec.size() - 1 ? vec[position + 1] : vec.front());
}

inline size_t previousModulo(const size_t number, const size_t modulus) {
  return (number == 0 ? modulus - 1 : number - 1);
}

inline size_t successiveModulo(const size_t number, const size_t modulus) {
  return (number == modulus - 1 ? 0 : number + 1);
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
