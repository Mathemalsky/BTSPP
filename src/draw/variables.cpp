#include "draw/variables.hpp"

#include <array>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

// graph library
#include "geometry.hpp"
#include "graph.hpp"

#include "utility/utils.hpp"

#include "solve/exactsolver.hpp"

namespace drawing {
graph::Euclidean EUCLIDEAN;

bool SHOW_DEBUG_WINDOW;
bool SHOW_SETTINGS_WINDOW;

std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;
bool BTSP_DRAW_OPEN_EAR_DECOMPOSITION;
bool BTSP_DRAW_BICONNECTED_GRAPH;
bool BTSP_DRAW_HAMILTON_CYCLE;
bool BTSPP_DRAW_BICONNECTED_GRAPH;
bool BTSPP_DRAW_HAMILTON_PATH;
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

namespace solve {
bool BTSP_FORBID_CROSSING;
std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace solve
