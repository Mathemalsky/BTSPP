#pragma once

#include <unordered_map>

#include "draw/definitions.hpp"

#include "graph/geometry.hpp"

namespace graph {
extern Data<Point2D> POINTS;
}

namespace imguiwindow {
extern bool SHOW_SETTINGS_WINDOW;
}

namespace input {
extern std::unordered_map<int, bool> STATE;
}  // namespace input

namespace mainwindow {
extern int HEIGHT;
extern int WIDTH;
}  // namespace mainwindow

namespace mouse {
extern double x;
extern double y;
}  // namespace mouse
