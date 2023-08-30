/*
 * BTSPP is a tool to solve, approximate and draw instances of BTSVPP,
 * BTSPP, BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "draw/variables.hpp"

#include <array>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

// graph library
#include "geometry.hpp"
#include "graph.hpp"

#include "utility/utils.hpp"

#include "solve/exactsolver.hpp"

namespace drawing {
graph::Euclidean EUCLIDEAN;

bool SHOW_DEBUG_WINDOW;
bool SHOW_SETTINGS_WINDOW;
bool COLLAPSE_SETTINGS_WINDOW;

std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;
bool BTSP_DRAW_OPEN_EAR_DECOMPOSITION;
bool BTSP_DRAW_BICONNECTED_GRAPH;
bool BTSP_DRAW_HAMILTON_CYCLE;
bool BTSPP_DRAW_BICONNECTED_GRAPH;
bool BTSPP_DRAW_HAMILTON_PATH;
bool BTSVPP_DRAW_BICONNECTED_GRAPH;
bool BTSVPP_DRAW_HAMILTON_PATH;
}  // namespace drawing

namespace input {
std::unordered_map<int, bool> STATE;

namespace mouse {
graph::Point2D MOUSE_LEFT_CLICKED;
int NODE_IN_MOTION;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
int HEIGHT;
int WIDTH;
}  // namespace mainwindow

namespace solve {
bool BTSP_FORBID_CROSSING;
std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace solve
