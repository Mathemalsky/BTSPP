#pragma once

#include <cstddef>
#include <vector>

template <typename T>
static inline size_t bisectionSearch(
    typename std::vector<T>::iterator it, typename std::vector<T>::iterator end(), const T& search) {
}

template <typename T>
class Storage {
  Storage(const size_t outer) { pOuterIndeces.reserve(outer); }

  void insert(const size_t outer, const size_t inner, const T& value);

private:
  std::vector<size_t> pOuterIndeces;
  std::vector<size_t> pInnerIndeces;
  std::vector<T> pValues;
};

template <typename T>
void Storage<T>::insert(const size_t outer, const size_t inner, const T& value) {
}
