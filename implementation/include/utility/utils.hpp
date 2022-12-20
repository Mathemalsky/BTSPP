#pragma once

#include <cstddef>
#include <vector>

template <typename Type>
size_t bytes_of(const std::vector<Type>& vec) {
  return vec.size() * sizeof(Type);
}
