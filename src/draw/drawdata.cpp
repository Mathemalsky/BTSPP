#include "draw/drawdata.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

/***********************************************************************************************************************
 *                                           callback helper functions
 **********************************************************************************************************************/
// copy here is quickfix, should be moved into drawdata.hpp
static graph::Point2D transformCoordinates(const double x, const double y) {
  return graph::Point2D{2.0 * x / mainwindow::WIDTH - 1.0, -2.0 * y / mainwindow::HEIGHT + 1.0};
}

static void toggle(bool& value) {
  value = !value;
}

static void selectNodeToMove(GLFWwindow* window) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  input::mouse::MOUSE_LEFT_CLICKED = transformCoordinates(x, y);
  for (unsigned int i = 0; i < drawing::EUCLIDEAN.vertices().size(); ++i) {
    if (norm2(drawing::EUCLIDEAN.vertices()[i] - input::mouse::MOUSE_LEFT_CLICKED) < drawing::VETREX_RADIUS) {
      input::mouse::NODE_IN_MOTION = i;
      break;  // move only one point at one at a time
    }
  }
}

static void cycleBTSPApproxDisplay() {
  if (drawing::BTSP_DRAW_BICONNECTED_GRAPH) {
    drawing::BTSP_DRAW_BICONNECTED_GRAPH      = false;
    drawing::BTSP_DRAW_OPEN_EAR_DECOMPOSITION = true;
  }
  else if (drawing::BTSP_DRAW_OPEN_EAR_DECOMPOSITION) {
    drawing::BTSP_DRAW_OPEN_EAR_DECOMPOSITION = false;
    drawing::BTSP_DRAW_HAMILTON_CYCLE         = true;
  }
  else if (drawing::BTSP_DRAW_HAMILTON_CYCLE) {
    drawing::BTSP_DRAW_HAMILTON_CYCLE    = false;
    drawing::BTSP_DRAW_BICONNECTED_GRAPH = true;
  }
  else {
    drawing::BTSP_DRAW_BICONNECTED_GRAPH = true;
  }
}

static void cycleBTSPPApproxDisplay() {
  if (drawing::BTSPP_DRAW_BICONNECTED_GRAPH) {
    drawing::BTSPP_DRAW_BICONNECTED_GRAPH = false;
    drawing::BTSPP_DRAW_HAMILTON_PATH     = true;
  }
  else if (drawing::BTSPP_DRAW_HAMILTON_PATH) {
    drawing::BTSPP_DRAW_HAMILTON_PATH     = false;
    drawing::BTSPP_DRAW_BICONNECTED_GRAPH = true;
  }
  else {
    drawing::BTSPP_DRAW_BICONNECTED_GRAPH = true;
  }
}

/***********************************************************************************************************************
 *                                                      callbacks
 **********************************************************************************************************************/

void DrawData::keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action,
                           [[maybe_unused]] int mods) {
  if (key == GLFW_KEY_F1 && action == GLFW_PRESS) {
    toggle(drawing::SHOW_DEBUG_WINDOW);
  }
  if (key == GLFW_KEY_F3 && action == GLFW_PRESS) {
    toggle(drawing::SHOW_SETTINGS_WINDOW);
  }
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }
  if (key == GLFW_KEY_R && action == GLFW_PRESS) {
    for (const ProblemType& type : problemType::PROBLEM_TYPES) {
      if (drawing::ACTIVE[std::to_underlying(type)]) {
        slowEvents::SOLVE[std::to_underlying(type)] = true;
      }
    }
  }
  if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_approx)]);
  }
  if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::BTSPP_approx)]);
  }
  if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_exact)]);
  }
  if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::BTSPP_exact)]);
  }
  if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::TSP_exact)]);
  }
  if (key == GLFW_KEY_T && action == GLFW_PRESS) {
    if (drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_approx)]) {
      cycleBTSPApproxDisplay();
    }
    if (drawing::ACTIVE[std::to_underlying(ProblemType::BTSPP_approx)]) {
      cycleBTSPPApproxDisplay();
    }
  }
}

void DrawData::mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods) {
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    input::STATE[GLFW_MOUSE_BUTTON_LEFT] = true;
    selectNodeToMove(window);
  }
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    input::STATE[GLFW_MOUSE_BUTTON_LEFT] = false;
    input::mouse::NODE_IN_MOTION         = namedInts::INVALID;
  }
}