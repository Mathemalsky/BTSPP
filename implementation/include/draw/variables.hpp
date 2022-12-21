#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "draw/definitions.hpp"

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

namespace graph {
extern Euclidean EUCLIDEAN;
;
extern std::vector<float> POINTS_F;
extern std::vector<unsigned int> TOUR;
extern std::vector<unsigned int> TOUR_DRAW_CYCLE;

void updatePointsfFromEuclidean();
void initPointsfFromEuclidean();

void updateTourDrawCycleFromTour();
void initTourDrawCycleFromTour();
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

namespace slowEvents {
extern std::unordered_map<ProblemType, bool> SOLVE;
extern std::unordered_map<ProblemType, bool> ACTIVE;
}  // namespace slowEvents
