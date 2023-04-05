#include <iostream>

#include "exception/exceptions.hpp"

#include "utility/utils.hpp"

#include "commandinterpreter.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  try {
    interpretCommandLine(argc, argv);
  }
  catch (const Exception& error) {
    printLightred("Error");
    std::cerr << ": " << error.what() << std::endl;
    return -1;
  }
  return 0;
}
