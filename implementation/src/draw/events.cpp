#include "draw/events.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/variables.hpp"

#include "exactsolver.hpp"

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
  }
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    input::STATE[GLFW_MOUSE_BUTTON_LEFT] = false;
  }
}

/***********************************************************************************************************************
 *                                               fast events
 **********************************************************************************************************************/

static void moveNode(GLFWwindow* window, const Buffers& buffers) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  const Point2D oldMousePos  = transformCoordinates(input::mouse::x, input::mouse::y);
  const Point2D diffMousePos = transformCoordinates(x, y) - oldMousePos;
  for (Point2D& point : drawing::EUCLIDEAN.verteces()) {
    if (norm2(point - oldMousePos) < drawing::VETREX_RADIUS) {
      point += diffMousePos;
      break;  // move only one point at one at a time
    }
  }
  drawing::updatePointsfFromEuclidean();
  buffers.coordinates.bufferSubData(drawing::POINTS_F);
  buffers.tourCoordinates.bufferSubData(drawing::POINTS_F);
}

static void handleFastEvents(GLFWwindow* window, const Buffers& buffers) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT]) {
    moveNode(window, buffers);
  }
  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}

/***********************************************************************************************************************
 *                                               slow events
 **********************************************************************************************************************/

static void handleSlowEvents([[maybe_unused]] const Buffers& buffers) {
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_exact)] = false;
    exactsolver::Result res = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSP_exact);
    drawing::updateOrder(res.tour, ProblemType::BTSP_exact);
  }
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::TSP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::TSP_exact)] = false;
    exactsolver::Result res = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::TSP_exact);
    drawing::updateOrder(res.tour, ProblemType::TSP_exact);
  }
}

void handleEvents(GLFWwindow* window, const Buffers& buffers) {
  handleFastEvents(window, buffers);
  handleSlowEvents(buffers);
}
