#pragma once

#include <cassert>
#include <cstddef>
#include <fstream>
#include <vector>

template <typename Data>
class DataIterator {
public:
  using ValueType = typename Data::ValueType;
  using Pointer   = ValueType*;
  using Reference = ValueType&;

  DataIterator(Pointer ptr) : pPtr(ptr) {}

  DataIterator& operator++() {
    ++pPtr;
    return *this;
  }

  DataIterator operator++(int) {
    DataIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  Reference operator*() const { return *pPtr; }

  bool operator==(const DataIterator& other) const { return pPtr == other.pPtr; }

  bool operator!=(const DataIterator& other) const { return !(pPtr == other.pPtr); }

private:
  Pointer pPtr;
};

/*!
 * \brief Data class is designed as simple non-resizable container.
 * \details Data is a pointer of any type and stores the number of objects in there.
 * Additionally it tracks if it should destruct the pointer when leaving scope.
 */
template <typename T>
class Data {
public:
  using ValueType = T;
  using Iterator  = DataIterator<Data<ValueType>>;

  Data() = default;
  Data(ValueType* pointer, const unsigned int size, bool toDelete = false)
    : pData(pointer), pSize(size), pDestructionNeeded(toDelete) {}

  explicit Data(std::vector<ValueType>& vec) : pData(&vec[0]), pSize(vec.size()), pDestructionNeeded(false) {}

  ~Data() {
    if (pDestructionNeeded) {
      delete[] pData;
    }
  }

  Iterator begin() const { return Iterator(pData); }

  Iterator end() const { return Iterator(pData + pSize); }

  ValueType& operator[](unsigned int index) {
    assert(index < pSize && "[Data] trying to access out of range");
    return pData[index];
  }
  const ValueType& operator[](unsigned int index) const {
    assert(index < pSize && "[Data] trying to access out of range");
    return pData[index];
  }

  ValueType* data() { return pData; }
  size_t size() const { return pSize; }
  size_t byteSize() const { return pSize * sizeof(ValueType); }

private:
  ValueType* pData;        /*!< pointer to memory */
  size_t pSize;            /*!< number of objects at this memory address*/
  bool pDestructionNeeded; /*!< stores if this pointer needs to be destructed when leaving scope */
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Data<T>& data) {
  for (const T& dat : data) {
    os << dat << " ";
  }
  return os << std::endl;
}
