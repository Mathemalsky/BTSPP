#include "draw/events.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/variables.hpp"

#include "solve/exactsolver.hpp"

/***********************************************************************************************************************
 *                                               fast events
 **********************************************************************************************************************/

static Point2D transformCoordinates(const double x, const double y) {
  return Point2D{2.0 * x / mainwindow::WIDTH - 1.0, -2.0 * y / mainwindow::HEIGHT + 1.0};
}

static void moveNode(GLFWwindow* window, const Buffers& buffers) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  drawing::EUCLIDEAN.verteces()[input::mouse::NODE_IN_MOTION] = transformCoordinates(x, y);
  drawing::updatePointsfFromEuclidean();
  buffers.coordinates.bufferSubData(drawing::POINTS_F);
  buffers.tourCoordinates.bufferSubData(drawing::POINTS_F);
}

static void handleFastEvents(GLFWwindow* window, const Buffers& buffers) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT] && input::mouse::NODE_IN_MOTION != namedInts::INVALID) {
    moveNode(window, buffers);
  }
}

/***********************************************************************************************************************
 *                                               slow events
 **********************************************************************************************************************/

static void handleSlowEvents([[maybe_unused]] const Buffers& buffers) {
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_approx)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_approx)] = false;
    drawing::BTSP_APPROX_RESULT = approximation::approximate(drawing::EUCLIDEAN, ProblemType::BTSP_approx);
    drawing::INITIALIZED[static_cast<unsigned int>(ProblemType::BTSP_approx)] = true;
  }
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_exact)] = false;
    exactsolver::Result res = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSP_exact);
    drawing::updateOrder(res.tour, ProblemType::BTSP_exact);
    drawing::BTSP_EXACT_RESULT = res;
  }
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSPP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSPP_exact)] = false;
    exactsolver::Result res = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSPP_exact);
    drawing::updateOrder(res.tour, ProblemType::BTSPP_exact);
    drawing::BTSPP_EXACT_RESULT = res;
  }
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::TSP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::TSP_exact)] = false;
    exactsolver::Result res = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::TSP_exact);
    drawing::updateOrder(res.tour, ProblemType::TSP_exact);
  }
}

/***********************************************************************************************************************
 *                                           callback helper functions
 **********************************************************************************************************************/

static void toggleSettingsWindow() {
  if (drawing::SHOW_SETTINGS_WINDOW) {
    drawing::SHOW_SETTINGS_WINDOW = false;
  }
  else {
    drawing::SHOW_SETTINGS_WINDOW = true;
  }
}

static void toggleDebugWindow() {
  if (drawing::SHOW_DEBUG_WINDOW) {
    drawing::SHOW_DEBUG_WINDOW = false;
  }
  else {
    drawing::SHOW_DEBUG_WINDOW = true;
  }
}

static void selectNodeToMove(GLFWwindow* window) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  input::mouse::MOUSE_LEFT_CLICKED = transformCoordinates(x, y);
  for (unsigned int i = 0; i < drawing::EUCLIDEAN.verteces().size(); ++i) {
    if (norm2(drawing::EUCLIDEAN.verteces()[i] - input::mouse::MOUSE_LEFT_CLICKED) < drawing::VETREX_RADIUS) {
      input::mouse::NODE_IN_MOTION = i;
      break;  // move only one point at one at a time
    }
  }
}

/***********************************************************************************************************************
 *                                                 callbacks
 **********************************************************************************************************************/

void keyCallback(
    [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action,
    [[maybe_unused]] int mods) {
  if (key == GLFW_KEY_F1 && action == GLFW_PRESS) {
    toggleDebugWindow();
  }
  if (key == GLFW_KEY_F3 && action == GLFW_PRESS) {
    toggleSettingsWindow();
  }
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }
  if (key == GLFW_KEY_R && action == GLFW_PRESS) {
    for (const ProblemType& type : problemType::PROBLEM_TYPES) {
      if (drawing::ACTIVE[(unsigned int) type]) {
        slowEvents::SOLVE[(unsigned int) type] = true;
      }
    }
  }
}

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods) {
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    input::STATE[GLFW_MOUSE_BUTTON_LEFT] = true;
    selectNodeToMove(window);
  }
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    input::STATE[GLFW_MOUSE_BUTTON_LEFT] = false;
    input::mouse::NODE_IN_MOTION         = namedInts::INVALID;
  }
}

/***********************************************************************************************************************
 *                                           event handling
 **********************************************************************************************************************/

void handleEvents(GLFWwindow* window, const Buffers& buffers) {
  handleFastEvents(window, buffers);
  handleSlowEvents(buffers);
}
