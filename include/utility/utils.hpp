#pragma once

#include <algorithm>
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
