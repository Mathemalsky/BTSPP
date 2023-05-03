/*
 * pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
 * BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <iostream>

// graph library
#include "exceptions.hpp"

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
  catch (const graph::Exception& error) {
    printLightred("Error");
    std::cerr << ": " << error.what() << std::endl;
    return -1;
  }
  return 0;
}
