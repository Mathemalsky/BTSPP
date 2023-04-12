#pragma once

#include <vector>

#include "graph/graph.hpp"

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
