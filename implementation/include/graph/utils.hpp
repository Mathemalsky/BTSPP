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
}  // namespace graph