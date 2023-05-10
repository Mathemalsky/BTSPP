/*
 * pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
 * BTSP and TSP. Drawing is limited to euclidean graphs.
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
#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>
#include <vector>

// graph library
#include "geometry.hpp"
#include "graph.hpp"

#include "draw/definitions.hpp"

#include "solve/approximation.hpp"
#include "solve/exactsolver.hpp"

namespace drawing {
extern graph::Euclidean EUCLIDEAN;

extern bool SHOW_SETTINGS_WINDOW;
extern bool SHOW_DEBUG_WINDOW;

extern std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> ACTIVE;

extern bool BTSP_DRAW_OPEN_EAR_DECOMPOSITION;
extern bool BTSP_DRAW_BICONNECTED_GRAPH;
extern bool BTSP_DRAW_HAMILTON_CYCLE;
extern bool BTSPP_DRAW_BICONNECTED_GRAPH;
extern bool BTSPP_DRAW_HAMILTON_PATH;
}  // namespace drawing

namespace input {
extern std::unordered_map<int, bool> STATE;

namespace mouse {
extern graph::Point2D MOUSE_LEFT_CLICKED;
extern int NODE_IN_MOTION;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
extern int HEIGHT;
extern int WIDTH;
}  // namespace mainwindow

namespace solve {
extern bool BTSP_FORBID_CROSSING;
extern std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> SOLVE;
}  // namespace solve
