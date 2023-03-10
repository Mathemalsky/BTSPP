#pragma once

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
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
