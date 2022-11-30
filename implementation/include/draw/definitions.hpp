#pragma once

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

template <class T>
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

  T& operator[](unsigned int index) { return pData[index]; }
  T operator[](unsigned int index) const { return pData[index]; }

  T* data() { return pData; }
  unsigned int size() const { return pSize; }
  unsigned int byteSize() const { return pSize * sizeof(T); }

private:
  T* pData;
  unsigned int pSize;
  bool pDestructionNeeded;
};

template <typename T>
Data<T> toData(const std::vector<T>& vec) {
  const unsigned int size = vec.size();
  T* data                 = new T[size];
  std::memcpy(data, &vec[0], size * sizeof(T));
  return Data(data, size, true);
}
