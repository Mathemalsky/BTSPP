#pragma once

#include <vector>

namespace imguiwindow {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = true;
}

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 720;
constexpr unsigned int INITIAL_WIDTH  = 1280;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow

template <class T>
class Data {
public:
  Data(const T* pointer, const unsigned int size) : pData(pointer), pSize(size), pDestructionNeeded(false) {}
  explicit Data(const std::vector<T>& vec) : pData(&vec[0]), pSize(vec.size()), pDestructionNeeded(false) {}

  ~Data() {
    if (pDestructionNeeded) {
      delete[] pData;
    }
  }

  void operator=(const std::vector<T>& vec) {
    pSize              = vec.size();
    pData              = new T[pSize];
    pDestructionNeeded = true;
    // std::memcpy()
  }

  T& operator[](unsigned int index) { return pData[index]; }
  T operator[](unsigned int index) const { return pData[index]; }

  T* data() { return pData; }
  unsigned int size() const { return pSize; }
  unsigned int byteSize() const { return pSize * sizeof(T); }

private:
  const T* pData;
  unsigned int pSize;
  bool pDestructionNeeded;
};
