# pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
# BTSP and TSP. Drawing is limited to euclidean graphs.
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

project(BTSP_Test_Suite)

# specify where the output should be compiled
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23 -O0 -g -Wall -Wextra -pedantic")

OPTION(Visualization "Visualisation" ON)
message(STATUS "Visualisation=${Visualisation}")

if(${Visualisation} STREQUAL ON)
  # lists all sourcefiles to be compiled with the project
  file(GLOB SOURCES "src/*.cpp" "src/draw/*.cpp" "src/graph/*.cpp" "src/solve/*.cpp")

  # lists all header files to be included in the project
  file(GLOB HEADERS "include/*.hpp" "include/draw/*.hpp" "include/graph/*.hpp" "include/solve/*.hpp" "include/utility/*.hpp")

  # find the libraries GLEW and glfw3
  find_package(GLEW REQUIRED)
  find_package(glfw3 REQUIRED)

  # compile the imgui library
  add_subdirectory(lib/imgui)

  # pass the visalization argument to preprocessor to cut of unnecessary code
  add_compile_definitions(VISUALISATION=1)
else()
  # lists all sourcefiles to be compiled with the project
  file(GLOB SOURCES "src/*.cpp" "src/graph/*.cpp" "src/solve/*.cpp")

  # lists all header files to be included in the project
  file(GLOB HEADERS "include/*.hpp" "include/graph/*.hpp" "include/solve/*.hpp" "include/utility/*.hpp")
endif()

# include graph library
add_subdirectory(lib/graph)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# check include directory for included files
target_include_directories(${PROJECT_NAME} PUBLIC include)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(HIGHS_DIR ${CMAKE_SYSTEM_PREFIX_PATH}/lib/cmake/highs)
find_package(HIGHS REQUIRED)
find_package(Threads)

if(${Visualisation} STREQUAL ON)
  target_link_libraries (${PROJECT_NAME} PRIVATE glfw GLEW::GLEW IMGUI GRAPH Eigen3::Eigen libhighs)
else()
  # glfw, GLEW and IMGUI not needed
  target_link_libraries (${PROJECT_NAME} PRIVATE GRAPH Eigen3::Eigen libhighs)
endif()
