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