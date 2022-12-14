#include "draw/variables.hpp"

#include <cstddef>
#include <vector>
#include <unordered_map>

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

namespace graph {
Data<Point2D> POINTS;
Data<float> POINTS_F;
std::vector<size_t> TOUR;

void initPointsfFromPoints() {
  float* pointsF = new float[2 * POINTS.size()];
  for (size_t i = 0; i < POINTS.size(); ++i) {
    pointsF[2 * i]     = (float) POINTS[i].x;
    pointsF[2 * i + 1] = (float) POINTS[i].y;
  }
  POINTS_F.set(pointsF, 2 * POINTS.size(), true);
}

void updatePointsfFromPoints() {
  for (size_t i = 0; i < POINTS.size(); ++i) {
    POINTS_F[2 * i]     = (float) POINTS[i].x;
    POINTS_F[2 * i + 1] = (float) POINTS[i].y;
  }
}
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
