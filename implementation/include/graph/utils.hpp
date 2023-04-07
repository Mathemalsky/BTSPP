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
Type previousInCycle(const std::vector<Type>& vec, const size_t position) {
  return (position != 0 ? vec[position - 1] : vec.back());
}

template <typename Type>
Type successiveInCycle(const std::vector<Type>& vec, const size_t position) {
  return (position != vec.size() - 1 ? vec[position + 1] : vec.front());
}
}  // namespace graph