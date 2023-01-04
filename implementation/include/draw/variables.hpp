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
extern std::unordered_map<ProblemType, std::vector<unsigned int>> ORDER;

void updatePointsfFromEuclidean();

void updateOrder(const std::vector<unsigned int>& order, const ProblemType& type);
}  // namespace graph

namespace imguiwindow {
extern bool SHOW_SETTINGS_WINDOW;
extern std::unordered_map<ProblemType, bool> ACTIVE;
extern std::unordered_map<ProblemType, std::array<float, 4>> COLOR;  // 4th float currently unused
extern std::unordered_map<ProblemType, float> THICKNESS;
extern std::array<float, 4> VERTEX_COLOR;  // 4th float currently unused

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
extern std::unordered_map<ProblemType, bool> SOLVE;
}  // namespace slowEvents
