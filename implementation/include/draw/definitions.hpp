#pragma once

#include <array>

#include "solve/definitions.hpp"

using RGBA_COLOUR = std::array<float, 4>;  // 4th float currently unused

namespace namedInts {
constexpr int INVALID = -1;
}

namespace problemType {
constexpr std::array<ProblemType, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> PROBLEM_TYPES = {
    ProblemType::BTSP_approx, ProblemType::BTSP_exact, ProblemType::BTSPP_exact, ProblemType::TSP_exact};
}

namespace drawing {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = true;
constexpr bool INITIAL_SHOW_DEBUG_WINDOW    = false;
constexpr float VETREX_RADIUS               = 0.01f;
constexpr unsigned int PATH_OVERHEAD        = 3;

constexpr std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_ACTIVENESS{
    false,  // BTSP_approx
    false,  // BTSP_exact
    false,  // BTSPP_axact
    false   // TSP_exact
};

constexpr std::array<RGBA_COLOUR, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_COLOUR{
    RGBA_COLOUR{0.0f, 1.0f, 0.0f, 1.0f},  // BTSP_approx
    RGBA_COLOUR{1.0f, 1.0f, 0.0f, 1.0f},  // BTSP_exact
    RGBA_COLOUR{0.3f, 0.7f, 0.2f, 1.0f},  // BTSPP_exact
    RGBA_COLOUR{1.0f, 0.0f, 0.0f, 1.0f}   // TSP_exact
};

constexpr bool INITIAL_DRAW_BICONNECTED_GRAPH      = false;
constexpr bool INITIAL_DRAW_OPEN_EAR_DECOMPOSITION = false;

constexpr std::array<float, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_THICKNESS{
    6.0f,  // BTSP_approx
    5.0f,  // BTSP_exact
    5.0f,  // BTSPP_exact
    3.0f   // TSP_exact
};

constexpr RGBA_COLOUR INITIAL_CLEAR_COLOUR  = {0.19, 0.19, 0.25, 1.0};
constexpr RGBA_COLOUR INITIAL_VERTEX_COLOUR = {1.0f, 0.0f, 0.0f, 1.0f};
}  // namespace drawing

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 960;
constexpr unsigned int INITIAL_WIDTH  = 960;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow
