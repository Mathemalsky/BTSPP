cmake_minimum_required(VERSION 3.10)

project(BTSP_Test_Suite)

# specify where the output should be compiled
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0 -g -Wall -Wextra -pedantic -fPIE")

# lists all sourcefiles to be compiled with the project
file(GLOB SOURCES "src/*.cpp" "src/draw/*.cpp" "src/graph/*.cpp")

# lists all header files to be included in the project
file(GLOB HEADERS "include/*.hpp" "include/draw/*.hpp" "include/graph/*.hpp" "include/utility/*.hpp")

add_subdirectory(lib/imgui)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# check include directory for included files
target_include_directories(${PROJECT_NAME} PUBLIC include)

find_package(GLEW REQUIRED)

find_package(glfw3 REQUIRED)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(HIGHS_DIR ${CMAKE_SYSTEM_PREFIX_PATH}/lib/cmake/highs)
find_package(HIGHS REQUIRED)
find_package(Threads)

target_link_libraries (${PROJECT_NAME} PRIVATE glfw GLEW::GLEW IMGUI Eigen3::Eigen libhighs)
