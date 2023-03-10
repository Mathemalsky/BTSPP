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
    drawing::updateOrder(drawing::BTSP_APPROX_RESULT.tour, ProblemType::BTSP_approx);
  }
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSP_exact)] = false;
    drawing::BTSP_EXACT_RESULT = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSP_exact);
    drawing::updateOrder(drawing::BTSP_EXACT_RESULT.tour, ProblemType::BTSP_exact);
  }
  if (slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSPP_exact)]) {
    slowEvents::SOLVE[static_cast<unsigned int>(ProblemType::BTSPP_exact)] = false;
    drawing::BTSPP_EXACT_RESULT = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSPP_exact);
    drawing::updateOrder(drawing::BTSPP_EXACT_RESULT.tour, ProblemType::BTSPP_exact);
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

static void toggle(bool& value) {
  value = !value;
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

static void cycleBTSPApproxDisplay() {
  if (drawing::DRAW_BICONNECTED_GRAPH) {
    drawing::DRAW_BICONNECTED_GRAPH      = false;
    drawing::DRAW_OPEN_EAR_DECOMPOSITION = true;
  }
  else if (drawing::DRAW_OPEN_EAR_DECOMPOSITION) {
    drawing::DRAW_OPEN_EAR_DECOMPOSITION = false;
    drawing::DRAW_HAMILTON_CYCLE         = true;
  }
  else if (drawing::DRAW_HAMILTON_CYCLE) {
    drawing::DRAW_HAMILTON_CYCLE    = false;
    drawing::DRAW_BICONNECTED_GRAPH = true;
  }
  else {
    drawing::DRAW_BICONNECTED_GRAPH = true;
  }
}

/***********************************************************************************************************************
 *                                                      callbacks
 **********************************************************************************************************************/

void keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action,
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
  if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_exact)]);
  }
  if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::BTSPP_exact)]);
  }
  if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
    toggle(drawing::ACTIVE[std::to_underlying(ProblemType::TSP_exact)]);
  }
  if (key == GLFW_KEY_T && action == GLFW_PRESS && drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_approx)]) {
    cycleBTSPApproxDisplay();
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
