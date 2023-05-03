# GRAPH is a library to store and manipulate graphs as adjacency list or
# as sparse eigen matrix. Different specialized types of graphs are
# supported.
# Copyright (C) 2023 Jurek Rostalsky
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
cmake_minimum_required(VERSION 3.20)
set(SOURCES
  "src/algorithm.cpp"
  "src/graph.cpp"
  )

set(HEADERS
  "include/algorithm.hpp"
  "include/exceptions.hpp"
  "include/geometry.hpp"
  "include/graph.hpp"
  "include/ostream.hpp"
  "include/utils.hpp"
  )

# find Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(GRAPH STATIC ${SOURCES} ${HEADERS} )
target_include_directories(GRAPH PUBLIC include)
target_link_libraries(GRAPH PRIVATE Eigen3::Eigen)