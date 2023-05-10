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

#include <fstream>
#include <vector>

namespace graph {
template <typename Type>
bool removeAnyElementByValue(std::vector<Type>& vec, const Type& val) {
  typename std::vector<Type>::iterator it = std::find(vec.begin(), vec.end(), val);
  if (it == vec.end()) {
    return false;
  }
  std::swap(*it, vec.back());
  vec.pop_back();
  return true;
}

template <typename Type>
size_t findPosition(const std::vector<Type>& vec, const Type& element) {
  return std::distance(vec.begin(), std::find(vec.begin(), vec.end(), element));
}

template <typename Type, size_t N>
size_t findPosition(const std::array<Type, N>& ary, const Type& element) {
  return std::distance(ary.begin(), std::find(ary.begin(), ary.end(), element));
}
}  // namespace graph