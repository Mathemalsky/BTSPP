#include <cstdlib>

#include "draw/visualization.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  return visualize(std::atoi(argv[1]));
}
