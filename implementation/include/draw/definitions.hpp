#pragma once

#include <array>
#include <unordered_map>

enum class ProblemType { BTSP_exact, TSP_exact };
namespace problemType {
constexpr std::array<ProblemType, 2> PROBLEM_TYPES = {ProblemType::BTSP_exact, ProblemType::TSP_exact};
}

namespace drawing {
constexpr float VETREX_RADIUS = 0.01f;
}

namespace imguiwindow {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = true;
const std::unordered_map<ProblemType, bool> INITIAL_ACTIVENESS{
  {ProblemType::BTSP_exact, true},
  {ProblemType::TSP_exact, false}};
}  // namespace imguiwindow

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 960;
constexpr unsigned int INITIAL_WIDTH  = 960;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow
