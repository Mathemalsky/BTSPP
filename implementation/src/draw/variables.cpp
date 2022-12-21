#include "draw/variables.hpp"

#include <cstddef>
#include <vector>
#include <unordered_map>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

#include "utility/utils.hpp"

namespace graph {
Euclidean EUCLIDEAN;
std::vector<float> POINTS_F;
std::vector<unsigned int> TOUR;
std::vector<unsigned int> TOUR_DRAW_CYCLE;

void updatePointsfFromEuclidean() {
  const std::vector<Point2D> points = EUCLIDEAN.verteces();
  for (size_t i = 0; i < points.size(); ++i) {
    POINTS_F[2 * i]     = (float) points[i].x;
    POINTS_F[2 * i + 1] = (float) points[i].y;
  }
}

void initPointsfFromEuclidean() {
  POINTS_F.resize(2 * EUCLIDEAN.numberOfNodes());
  updatePointsfFromEuclidean();
}

void updateTourDrawCycleFromTour() {
  std::memcpy(TOUR_DRAW_CYCLE.data(), graph::TOUR.data(), bytes_of(graph::TOUR));
  std::memcpy(TOUR_DRAW_CYCLE.data() + graph::TOUR.size(), graph::TOUR.data(), 3 * sizeof(unsigned int));
}

void initTourDrawCycleFromTour() {
  TOUR_DRAW_CYCLE.resize(TOUR.size() + 3);
  updateTourDrawCycleFromTour();
}

}  // namespace graph

namespace imguiwindow {
bool SHOW_SETTINGS_WINDOW;
std::unordered_map<ProblemType, bool> ACTIVE;
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
