#pragma once

#include <unordered_map>

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

namespace graph {
extern Data<Point2D> POINTS;
}

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
