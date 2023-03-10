#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

#include "solve/approximation.hpp"
#include "solve/exactsolver.hpp"

namespace drawing {
extern graph::Euclidean EUCLIDEAN;
extern std::vector<float> POINTS_F;

extern bool SHOW_SETTINGS_WINDOW;
extern bool SHOW_DEBUG_WINDOW;

extern std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;
extern std::array<RGBA_COLOUR, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> COLOUR;
extern bool DRAW_OPEN_EAR_DECOMPOSITION;
extern bool DRAW_BICONNECTED_GRAPH;
extern bool DRAW_HAMILTON_CYCLE;
extern std::array<std::vector<unsigned int>, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> ORDER;
extern std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIALIZED;
extern std::array<float, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> THICKNESS;
extern RGBA_COLOUR CLEAR_COLOUR;
extern RGBA_COLOUR VERTEX_COLOUR;

extern approximation::Result BTSP_APPROX_RESULT;
extern exactsolver::Result BTSP_EXACT_RESULT;
extern exactsolver::Result BTSPP_EXACT_RESULT;

void updatePointsfFromEuclidean();

void updateOrder(const std::vector<unsigned int>& order, const ProblemType& type);
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

namespace slowEvents {
extern std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace slowEvents
