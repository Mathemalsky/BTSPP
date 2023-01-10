#pragma once

#include <array>

enum class ProblemType : unsigned int { BTSP_approx, BTSP_exact, TSP_exact, NUMBER_OF_OPTIONS };
namespace problemType {
constexpr std::array<ProblemType, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> PROBLEM_TYPES = {
  ProblemType::BTSP_approx, ProblemType::BTSP_exact, ProblemType::TSP_exact};
}

namespace drawing {
constexpr float VETREX_RADIUS = 0.01f;
}

namespace imguiwindow {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = true;

constexpr std::array<bool, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> INITIAL_ACTIVENESS{
  false,  // BTSP_approx
  true,   // BTSP_exact
  false   // TSP_exact
};

constexpr std::array<std::array<float, 4>, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> INITIAL_COLOR{
  std::array<float, 4>{0.0f, 1.0f, 0.0f, 1.0f},  // BTSP_approx
  std::array<float, 4>{1.0f, 1.0f, 0.0f, 1.0f},  // BTSP_exact
  std::array<float, 4>{1.0f, 0.0f, 0.0f, 1.0f}   // TSP_exact
};

constexpr std::array<float, (unsigned int) ProblemType::NUMBER_OF_OPTIONS> INITIAL_THICKNESS{
  6.0f,  // BTSP_approx
  5.0f,  // BTSP_exact
  3.0f   // TSP_exact
};

constexpr std::array<float, 4> INITIAL_VERTEX_COLOR = {1.0f, 0.0f, 0.0f, 1.0f};
}  // namespace imguiwindow

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 960;
constexpr unsigned int INITIAL_WIDTH  = 960;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow
