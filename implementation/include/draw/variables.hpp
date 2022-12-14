#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

namespace graph {
extern Data<Point2D> POINTS;
extern Data<float> POINTS_F;
extern std::vector<size_t> TOUR;
extern std::vector<uint32_t> TOUR_32;

void initPointsfFromPoints();
void updatePointsfFromPoints();

void initTour32FromTour();
void updateTour32FromTour();
}  // namespace graph

namespace imguiwindow {
extern bool SHOW_SETTINGS_WINDOW;
}

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
