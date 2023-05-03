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
