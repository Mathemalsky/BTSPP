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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename Type>
size_t bytes_of(const std::vector<Type>& vec) {
  return vec.size() * sizeof(Type);
}

template <typename Type>
std::ostream& operator<<(std::ostream& os, const std::vector<Type>& vec) {
  const size_t size = vec.size();
  for (size_t i = 0; i < size; ++i) {
    os << vec[i] << (i == size - 1 ? "" : " ");
  }
  return os << std::endl;
}

/*!
 * @brief prints the word in red
 * @param word string to be printed
 */
inline void printLightred(std::string word) noexcept {
  std::cout << "\033[1;31m";
  std::cout << word;
  std::cout << "\033[0m";
}

/*!
 * @brief prints the word in yellow
 * @param word string to be printed
 */
inline void printYellow(std::string word) noexcept {
  std::cout << "\033[1;33m";
  std::cout << word;
  std::cout << "\033[0m";
}

class Stopwatch {
public:
  Stopwatch()  = default;
  ~Stopwatch() = default;

  void reset() { pStartTime = std::chrono::high_resolution_clock::now(); }

  double elapsedTimeInMilliseconds() const {
    const std::chrono::time_point<std::chrono::high_resolution_clock> now;
    return std::round(std::chrono::duration_cast<std::chrono::microseconds>(now - pStartTime).count()) / 1000.0;
  }

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> pStartTime;
};