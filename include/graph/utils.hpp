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