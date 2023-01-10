#include "draw/variables.hpp"

#include <array>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <utility>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

#include "utility/utils.hpp"

namespace graph {
Euclidean EUCLIDEAN;
std::vector<float> POINTS_F;
std::array<std::vector<unsigned int>, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> ORDER;

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
    ORDER[(unsigned int) type].resize(order.size() + 3);
    std::memcpy(ORDER[(unsigned int) type].data(), order.data(), bytes_of(order));
    std::memcpy(ORDER[(unsigned int) type].data() + order.size(), order.data(), 3 * sizeof(unsigned int));
  }
}
}  // namespace graph

namespace imguiwindow {
bool SHOW_SETTINGS_WINDOW;
std::array<bool, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> ACTIVE;
std::array<std::array<float, 4>, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> COLOR;  // 4th float currently unused
std::array<float, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> THICKNESS;
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
std::array<bool, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> SOLVE;
}  // namespace slowEvents
