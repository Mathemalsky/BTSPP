#include "draw/variables.hpp"

#include <array>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <utility>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

#include "utility/utils.hpp"

#include "solve/exactsolver.hpp"

namespace drawing {
Euclidean EUCLIDEAN;
std::vector<float> POINTS_F;

bool SHOW_DEBUG_WINDOW;
bool SHOW_SETTINGS_WINDOW;

std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;
std::array<RGBA_COLOUR, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> COLOUR;
bool DRAW_OPEN_EAR_DECOMPOSITION;
bool DRAW_BICONNECTED_GRAPH;
std::array<std::vector<unsigned int>, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> ORDER;
std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIALIZED;
std::array<float, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> THICKNESS;
RGBA_COLOUR VERTEX_COLOUR;

approximation::Result BTSP_APPROX_RESULT;
exactsolver::Result BTSP_EXACT_RESULT;
exactsolver::Result BTSPP_EXACT_RESULT;

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
  if (type == ProblemType::BTSPP_exact) {
    ORDER[(unsigned int) type].resize(order.size() + 2);
    ORDER[(unsigned int) type][0] = order.back();
    std::memcpy(ORDER[(unsigned int) type].data() + 1, order.data(), bytes_of(order));
    ORDER[(unsigned int) type].back() = order[0];
  }
  INITIALIZED[static_cast<unsigned int>(type)] = true;
}
}  // namespace drawing

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
std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace slowEvents
