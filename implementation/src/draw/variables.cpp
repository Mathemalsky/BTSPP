#include "draw/variables.hpp"

#include <array>
#include <cstddef>
#include <vector>
#include <unordered_map>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

#include "utility/utils.hpp"

namespace graph {
Euclidean EUCLIDEAN;
std::vector<float> POINTS_F;
std::unordered_map<ProblemType, std::vector<unsigned int>> ORDER;

void updatePointsfFromEuclidean() {
  POINTS_F.resize(2 * EUCLIDEAN.numberOfNodes());
  const std::vector<Point2D> points = EUCLIDEAN.verteces();
  for (size_t i = 0; i < points.size(); ++i) {
    POINTS_F[2 * i]     = (float) points[i].x;
    POINTS_F[2 * i + 1] = (float) points[i].y;
  }
}

void updateOrder(const std::vector<unsigned int>& order, const ProblemType& type) {
  if (type == ProblemType::BTSP_exact || type == ProblemType::TSP_exact) {
    ORDER[type].resize(order.size() + 3);
    std::memcpy(ORDER[type].data(), order.data(), bytes_of(order));
    std::memcpy(ORDER[type].data() + order.size(), order.data(), 3 * sizeof(unsigned int));
  }
}
}  // namespace graph

namespace imguiwindow {
bool SHOW_SETTINGS_WINDOW;
std::unordered_map<ProblemType, bool> ACTIVE;
std::unordered_map<ProblemType, std::array<float, 4>> COLOR;  // 4th float currently unused
std::unordered_map<ProblemType, float> THICKNESS;
std::array<float, 4> VERTEX_COLOR;  // 4th float currently unused
}  // namespace imguiwindow

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

namespace slowEvents {
std::unordered_map<ProblemType, bool> SOLVE;
}  // namespace slowEvents
