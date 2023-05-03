#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>
#include <vector>

#include "draw/definitions.hpp"

// graph library
#include "geometry.hpp"
#include "graph.hpp"

#include "solve/approximation.hpp"
#include "solve/exactsolver.hpp"

namespace drawing {
extern graph::Euclidean EUCLIDEAN;

extern bool SHOW_SETTINGS_WINDOW;
extern bool SHOW_DEBUG_WINDOW;

extern std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;

extern bool BTSP_DRAW_OPEN_EAR_DECOMPOSITION;
extern bool BTSP_DRAW_BICONNECTED_GRAPH;
extern bool BTSP_DRAW_HAMILTON_CYCLE;
extern bool BTSPP_DRAW_BICONNECTED_GRAPH;
extern bool BTSPP_DRAW_HAMILTON_PATH;
}  // namespace drawing

namespace input {
extern std::unordered_map<int, bool> STATE;

namespace mouse {
extern graph::Point2D MOUSE_LEFT_CLICKED;
extern int NODE_IN_MOTION;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
extern int HEIGHT;
extern int WIDTH;
}  // namespace mainwindow

namespace solve {
extern bool BTSP_FORBID_CROSSING;
extern std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace solve
