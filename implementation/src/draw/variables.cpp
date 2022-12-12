#include "draw/variables.hpp"

#include <cstddef>
#include <vector>
#include <unordered_map>

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

namespace graph {
Data<Point2D> POINTS;
std::vector<size_t> TOUR;
}  // namespace graph

namespace imguiwindow {
bool SHOW_SETTINGS_WINDOW;
}

namespace input {
std::unordered_map<int, bool> STATE;

namespace mouse {
double x;
double y;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
int HEIGHT;
int WIDTH;
}  // namespace mainwindow
