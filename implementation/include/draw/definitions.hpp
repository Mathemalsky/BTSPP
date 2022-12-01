#pragma once

#include <cassert>
#include <cstring>
#include <vector>

namespace imguiwindow {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = true;
}

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 720;
constexpr unsigned int INITIAL_WIDTH  = 1280;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow

/*!
 * \brief Data class is designed as simple non-resizable container.
 * \details Data is a pointer of any type and stores the number of objects in there.
 * Additionally it tracks if it should destruct the pointer when leaving scope.
 */
template <typename T>
class Data {
public:
  Data() = default;
  Data(T* pointer, const unsigned int size, bool toDelete = false)
    : pData(pointer), pSize(size), pDestructionNeeded(toDelete) {}

  explicit Data(std::vector<T>& vec) : pData(&vec[0]), pSize(vec.size()), pDestructionNeeded(false) {}

  ~Data() {
    if (pDestructionNeeded) {
      delete[] pData;
    }
  }

  T& operator[](unsigned int index) {
    assert(index < pSize && "[Data] trying to access out of range");
    return pData[index];
  }
  const T& operator[](unsigned int index) const {
    assert(index < pSize && "[Data] trying to access out of range");
    return pData[index];
  }

  T* data() { return pData; }
  unsigned int size() const { return pSize; }
  unsigned int byteSize() const { return pSize * sizeof(T); }

private:
  T* pData;                /*!< pointer to memory */
  unsigned int pSize;      /*!< number of objects at this memory address*/
  bool pDestructionNeeded; /*!< stores if this pointer needs to be destructed when leaving scope */
};
