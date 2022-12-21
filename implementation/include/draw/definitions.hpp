#pragma once

enum class Problem { BTSP, TSP };

namespace drawing {
constexpr float VETREX_RADIUS = 0.01f;
}

namespace imguiwindow {
constexpr bool INITIAL_SHOW_SETTINGS_WINDOW = false;
}

namespace mainwindow {
constexpr unsigned int INITIAL_HEIGHT = 960;
constexpr unsigned int INITIAL_WIDTH  = 960;
constexpr const char* NAME            = "BTSP";
}  // namespace mainwindow
