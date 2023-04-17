#include "draw/variables.hpp"

#include <array>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

#include "utility/utils.hpp"

#include "solve/exactsolver.hpp"

namespace drawing {
graph::Euclidean EUCLIDEAN;

bool SHOW_DEBUG_WINDOW;
bool SHOW_SETTINGS_WINDOW;

std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;
std::array<RGBA_COLOUR, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> COLOUR;
bool BTSP_DRAW_OPEN_EAR_DECOMPOSITION;
bool BTSP_DRAW_BICONNECTED_GRAPH;
bool BTSP_DRAW_HAMILTON_CYCLE;
bool BTSPP_DRAW_BICONNECTED_GRAPH;
bool BTSPP_DRAW_HAMILTON_PATH;
std::array<std::vector<unsigned int>, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> ORDER;
std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> INITIALIZED;
std::array<float, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> THICKNESS;
RGBA_COLOUR CLEAR_COLOUR;
RGBA_COLOUR VERTEX_COLOUR;

void updateOrder(const std::vector<unsigned int>& order, const ProblemType& type) {
  if (type == ProblemType::BTSP_approx || type == ProblemType::BTSP_exact || type == ProblemType::TSP_exact) {
    ORDER[std::to_underlying(type)].resize(order.size() + PATH_OVERHEAD);
    std::memcpy(ORDER[std::to_underlying(type)].data(), order.data(), bytes_of(order));
    std::memcpy(ORDER[std::to_underlying(type)].data() + order.size(), order.data(), PATH_OVERHEAD * sizeof(unsigned int));
  }
  else if (type == ProblemType::BTSPP_approx || type == ProblemType::BTSPP_exact) {
    ORDER[std::to_underlying(type)].resize(order.size() + PATH_OVERHEAD - 1);  // n-1 path segments to draw
    ORDER[std::to_underlying(type)][0] = order[1];
    std::memcpy(ORDER[std::to_underlying(type)].data() + 1, order.data(), bytes_of(order));
    ORDER[std::to_underlying(type)].back() = order[order.size() - 2];
  }
  INITIALIZED[std::to_underlying(type)] = true;
}
}  // namespace drawing

namespace input {
std::unordered_map<int, bool> STATE;

namespace mouse {
graph::Point2D MOUSE_LEFT_CLICKED;
int NODE_IN_MOTION;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
int HEIGHT;
int WIDTH;
}  // namespace mainwindow

namespace slowEvents {
std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace slowEvents
