#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

namespace graph {
extern Euclidean EUCLIDEAN;
extern std::vector<float> POINTS_F;
extern std::array<std::vector<unsigned int>, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> ORDER;

void updatePointsfFromEuclidean();

void updateOrder(const std::vector<unsigned int>& order, const ProblemType& type);
}  // namespace graph

namespace imguiwindow {
extern bool SHOW_SETTINGS_WINDOW;
extern std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;
extern std::array<RGBA_COLOUR, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)>
    COLOUR;  // 4th flt crrntly unused
extern std::array<float, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> THICKNESS;
extern RGBA_COLOUR VERTEX_COLOUR;  // 4th float currently unused
}  // namespace imguiwindow

namespace input {
extern std::unordered_map<int, bool> STATE;

namespace mouse {
extern double x;
extern double y;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
extern int HEIGHT;
extern int WIDTH;
}  // namespace mainwindow

namespace slowEvents {
extern std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace slowEvents
