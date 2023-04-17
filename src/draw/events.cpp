#include "draw/events.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/drawdata.hpp"
#include "draw/variables.hpp"

#include "solve/exactsolver.hpp"

using namespace drawing;

/***********************************************************************************************************************
 *                                               fast events
 **********************************************************************************************************************/

static graph::Point2D transformCoordinates(const double x, const double y) {
  return graph::Point2D{2.0 * x / mainwindow::WIDTH - 1.0, -2.0 * y / mainwindow::HEIGHT + 1.0};
}

static void moveNode(GLFWwindow* window, DrawData& drawData) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  drawing::EUCLIDEAN.vertices()[input::mouse::NODE_IN_MOTION] = transformCoordinates(x, y);
  drawData.floatVertices.updatePointsfFromEuclidean(drawing::EUCLIDEAN);
  drawData.buffers.coordinates.bufferSubData(drawData.floatVertices.read());
  drawData.buffers.tourCoordinates.bufferSubData(drawData.floatVertices.read());
}

static void handleFastEvents(GLFWwindow* window, DrawData& drawData) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT] && input::mouse::NODE_IN_MOTION != namedInts::INVALID) {
    moveNode(window, drawData);
  }
}

/***********************************************************************************************************************
 *                                               slow events
 **********************************************************************************************************************/

static void handleSlowEvents(DrawData& drawData) {
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSP_approx)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSP_approx)] = false;
    drawData.results.BTSP_APPROX_RESULT                             = approximation::approximateBTSP(drawing::EUCLIDEAN);
    drawData.vertexOrder.updateOrder(drawData.results.BTSP_APPROX_RESULT.tour, ProblemType::BTSP_approx);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_approx)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_approx)] = false;
    drawData.results.BTSPP_APPROX_RESULT                             = approximation::approximateBTSPP(drawing::EUCLIDEAN);
    drawData.vertexOrder.updateOrder(drawData.results.BTSPP_APPROX_RESULT.tour, ProblemType::BTSPP_approx);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSP_exact)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSP_exact)] = false;
    drawData.results.BTSP_EXACT_RESULT                             = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSP_exact);
    drawData.vertexOrder.updateOrder(drawData.results.BTSP_EXACT_RESULT.tour, ProblemType::BTSP_exact);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_exact)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_exact)] = false;
    drawData.results.BTSPP_EXACT_RESULT                             = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSPP_exact);
    drawData.vertexOrder.updateOrder(drawData.results.BTSPP_EXACT_RESULT.tour, ProblemType::BTSPP_exact);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::TSP_exact)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::TSP_exact)] = false;
    exactsolver::Result res                                       = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::TSP_exact);
    drawData.vertexOrder.updateOrder(res.tour, ProblemType::TSP_exact);
  }
}

/***********************************************************************************************************************
 *                                           event handling
 **********************************************************************************************************************/

void handleEvents(GLFWwindow* window, DrawData& drawData) {
  handleFastEvents(window, drawData);
  handleSlowEvents(drawData);
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

void keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods) {
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