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
#pragma once

/*!
 * \file
 * definitions of constants
 * This file is for defning compile time constants like initial values for global variables and typedefs. It shouldn't
 * depend on any other file from the draw directories.
 */

#include <array>

#include "solve/definitions.hpp"

using RGBA_COLOUR = std::array<float, 4>;  // 4th float currently unused

namespace namedInts {
constexpr int INVALID = -1;
}  // namespace namedInts

namespace problemType {
constexpr std::array<ProblemType, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> PROBLEM_TYPES = {ProblemType::BTSP_approx,
                                                                                                              ProblemType::BTSPP_approx,
                                                                                                              ProblemType::BTSVPP_approx,
                                                                                                              ProblemType::BTSP_exact,
                                                                                                              ProblemType::BTSPP_exact,
                                                                                                              ProblemType::TSP_exact};
}  // namespace problemType

namespace drawing {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW  = true;
constexpr bool INITIAL_SHOW_DEBUG_WINDOW     = false;
constexpr float VETREX_RADIUS                = 0.01f;
constexpr int CIRCLE_STEPS                   = 8;
constexpr float BOTLLENECK_EDGE_WIDTH_FACTOR = 2.0f;
constexpr unsigned int PATH_OVERHEAD         = 3;

constexpr std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_ACTIVENESS{
    false,  // BTSP_approx
    false,  // BTSPP_approx
    false,  // BTSVPP_approx
    false,  // BTSP_exact
    false,  // BTSPP_axact
    false   // TSP_exact
};

constexpr std::array<RGBA_COLOUR, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_COLOUR{
    RGBA_COLOUR{0.0f, 1.0f, 0.0f, 1.0f}, // BTSP_approx
    RGBA_COLOUR{1.0f, 0.0f, 1.0f, 1.0f}, // BTSPP_approx
    RGBA_COLOUR{1.0f, 0.5f, 0.0f, 1.0f}, // BTSVPP_approx
    RGBA_COLOUR{1.0f, 1.0f, 0.0f, 1.0f}, // BTSP_exact
    RGBA_COLOUR{0.3f, 0.7f, 0.2f, 1.0f}, // BTSPP_exact
    RGBA_COLOUR{1.0f, 0.0f, 0.0f, 1.0f}  // TSP_exact
};

constexpr bool INITIAL_BTSP_DRAW_BICONNECTED_GRAPH      = false;
constexpr bool INITIAL_BTSP_DRAW_OPEN_EAR_DECOMPOSITION = false;
constexpr bool INITIAL_BTSP_DRAW_HAMILTON_CYCLE         = false;
constexpr bool INITIAL_BTSPP_DRAW_BICONNECTED_GRAPH     = false;
constexpr bool INITIAL_BTSPP_DRAW_HAMILTON_PATH         = false;
constexpr bool INITIAL_BTSVPP_DRAW_BICONNECTED_GRAPH    = false;
constexpr bool INITIAL_BTSVPP_DRAW_HAMILTON_PATH        = false;

constexpr std::array<float, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_THICKNESS{
    4.0f,  // BTSP_approx
    4.0f,  // BTSPP_approx
    2.0f,  // BTSVPP_approx
    2.0f,  // BTSP_exact
    3.0f,  // BTSPP_exact
    3.0f   // TSP_exact
};

constexpr RGBA_COLOUR INITIAL_CLEAR_COLOUR  = {0.19, 0.19, 0.25, 1.0};
constexpr RGBA_COLOUR INITIAL_VERTEX_COLOUR = {1.0f, 1.0f, 0.0f, 1.0f};
}  // namespace drawing

namespace input {
namespace mouse {
constexpr int INITIAL_NODE_IN_MOTION = namedInts::INVALID;
}  // namespace mouse
}  // namespace input

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 1000;
constexpr unsigned int INITIAL_WIDTH  = 1000;
constexpr const char* NAME            = "BTSPP";
}  // namespace mainwindow
