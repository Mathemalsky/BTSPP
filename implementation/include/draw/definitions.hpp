#pragma once

#include <array>

using RGBA_COLOUR = std::array<float, 4>;  // 4th float currently unused

enum class ProblemType : unsigned int { BTSP_approx, BTSP_exact, TSP_exact, NUMBER_OF_OPTIONS };
namespace problemType {
constexpr std::array<ProblemType, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> PROBLEM_TYPES = {
    ProblemType::BTSP_approx, ProblemType::BTSP_exact, ProblemType::TSP_exact};
}

namespace drawing {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = true;
constexpr float VETREX_RADIUS               = 0.01f;

constexpr std::array<bool, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_ACTIVENESS{
    false,  // BTSP_approx
    false,  // BTSP_exact
    false   // TSP_exact
};

constexpr std::array<RGBA_COLOUR, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_COLOUR{
    RGBA_COLOUR{0.0f, 1.0f, 0.0f, 1.0f},  // BTSP_approx
    RGBA_COLOUR{1.0f, 1.0f, 0.0f, 1.0f},  // BTSP_exact
    RGBA_COLOUR{1.0f, 0.0f, 0.0f, 1.0f}   // TSP_exact
};

constexpr std::array<float, static_cast<unsigned int>(ProblemType::NUMBER_OF_OPTIONS)> INITIAL_THICKNESS{
    6.0f,  // BTSP_approx
    5.0f,  // BTSP_exact
    3.0f   // TSP_exact
};

constexpr RGBA_COLOUR INITIAL_VERTEX_COLOUR = {1.0f, 0.0f, 0.0f, 1.0f};
}  // namespace drawing

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 960;
constexpr unsigned int INITIAL_WIDTH  = 960;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow
