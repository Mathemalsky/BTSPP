#include "draw/events.hpp"

#include <memory>

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

static void moveNode(GLFWwindow* window, std::shared_ptr<drawing::DrawData> drawData) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  drawing::EUCLIDEAN.vertices()[input::mouse::NODE_IN_MOTION] = transformCoordinates(x, y);
  drawData->floatVertices.updatePointsfFromEuclidean(drawing::EUCLIDEAN);
  drawData->buffers.coordinates->bufferSubData(drawData->floatVertices.read());
  drawData->buffers.tourCoordinates->bufferSubData(drawData->floatVertices.read());
}

static void handleFastEvents(GLFWwindow* window, std::shared_ptr<drawing::DrawData> drawData) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT] && input::mouse::NODE_IN_MOTION != namedInts::INVALID) {
    moveNode(window, drawData);
  }
}

/***********************************************************************************************************************
 *                                               slow events
 **********************************************************************************************************************/

static void handleSlowEvents(std::shared_ptr<DrawData> drawData) {
  if (drawData->settings.solve[std::to_underlying(ProblemType::BTSP_approx)]) {
    drawData->settings.solve[std::to_underlying(ProblemType::BTSP_approx)] = false;
    drawData->results.BTSP_APPROX_RESULT                                   = approximation::approximateBTSP(drawing::EUCLIDEAN);
    drawData->vertexOrder.updateOrder(drawData->results.BTSP_APPROX_RESULT.tour, ProblemType::BTSP_approx);
  }
  if (drawData->settings.solve[std::to_underlying(ProblemType::BTSPP_approx)]) {
    drawData->settings.solve[std::to_underlying(ProblemType::BTSPP_approx)] = false;
    drawData->results.BTSPP_APPROX_RESULT                                   = approximation::approximateBTSPP(drawing::EUCLIDEAN);
    drawData->vertexOrder.updateOrder(drawData->results.BTSPP_APPROX_RESULT.tour, ProblemType::BTSPP_approx);
  }
  if (drawData->settings.solve[std::to_underlying(ProblemType::BTSP_exact)]) {
    drawData->settings.solve[std::to_underlying(ProblemType::BTSP_exact)] = false;
    drawData->results.BTSP_EXACT_RESULT = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSP_exact, solve::BTSP_FORBID_CROSSING);
    drawData->vertexOrder.updateOrder(drawData->results.BTSP_EXACT_RESULT.tour, ProblemType::BTSP_exact);
  }
  if (drawData->settings.solve[std::to_underlying(ProblemType::BTSPP_exact)]) {
    drawData->settings.solve[std::to_underlying(ProblemType::BTSPP_exact)] = false;
    drawData->results.BTSPP_EXACT_RESULT = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSPP_exact, solve::BTSP_FORBID_CROSSING);
    drawData->vertexOrder.updateOrder(drawData->results.BTSPP_EXACT_RESULT.tour, ProblemType::BTSPP_exact);
  }
  if (drawData->settings.solve[std::to_underlying(ProblemType::TSP_exact)]) {
    drawData->settings.solve[std::to_underlying(ProblemType::TSP_exact)] = false;
    exactsolver::Result res                                              = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::TSP_exact);
    drawData->vertexOrder.updateOrder(res.tour, ProblemType::TSP_exact);
  }
}

/***********************************************************************************************************************
 *                                           event handling
 **********************************************************************************************************************/

void handleEvents(GLFWwindow* window, std::shared_ptr<drawing::DrawData> drawData) {
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

static void cycleBTSPApproxDisplay(Settings& settings) {
  if (settings.BTSP_drawBiconnectedGraph) {
    settings.BTSP_drawBiconnectedGraph     = false;
    settings.BTSP_drawOpenEarDecomposition = true;
  }
  else if (settings.BTSP_drawOpenEarDecomposition) {
    settings.BTSP_drawOpenEarDecomposition = false;
    settings.BTSP_drawHamiltonCycle        = true;
  }
  else if (drawing::BTSP_DRAW_HAMILTON_CYCLE) {
    settings.BTSP_drawHamiltonCycle    = false;
    settings.BTSP_drawBiconnectedGraph = true;
  }
  else {
    settings.BTSP_drawBiconnectedGraph = true;
  }
}

static void cycleBTSPPApproxDisplay(Settings& settings) {
  if (settings.BTSPP_drawBiconnectedGraph) {
    settings.BTSPP_drawBiconnectedGraph = false;
    settings.BTSPP_drawHamiltonPath     = true;
  }
  else if (settings.BTSPP_drawHamiltonPath) {
    settings.BTSPP_drawHamiltonPath     = false;
    settings.BTSPP_drawBiconnectedGraph = true;
  }
  else {
    settings.BTSPP_drawBiconnectedGraph = true;
  }
}

/***********************************************************************************************************************
 *                                                      callbacks
 **********************************************************************************************************************/

void keyCallback([[maybe_unused]] GLFWwindow* window,
                 int key,
                 [[maybe_unused]] int scancode,
                 int action,
                 [[maybe_unused]] int mods,
                 std::shared_ptr<DrawData> drawData) {
  if (key == GLFW_KEY_F1 && action == GLFW_PRESS) {
    toggle(drawData->settings.showDebugWindow);
  }
  if (key == GLFW_KEY_F3 && action == GLFW_PRESS) {
    toggle(drawData->settings.showSettingsWindow);
  }
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }
  if (key == GLFW_KEY_R && action == GLFW_PRESS) {
    for (const ProblemType& type : problemType::PROBLEM_TYPES) {
      if (drawData->settings.activeness[std::to_underlying(type)]) {
        drawData->settings.solve[std::to_underlying(type)] = true;
      }
    }
  }
  if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
    toggle(drawData->settings.activeness[std::to_underlying(ProblemType::BTSP_approx)]);
  }
  if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
    toggle(drawData->settings.activeness[std::to_underlying(ProblemType::BTSPP_approx)]);
  }
  if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
    toggle(drawData->settings.activeness[std::to_underlying(ProblemType::BTSP_exact)]);
  }
  if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
    toggle(drawData->settings.activeness[std::to_underlying(ProblemType::BTSPP_exact)]);
  }
  if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
    toggle(drawData->settings.activeness[std::to_underlying(ProblemType::TSP_exact)]);
  }
  if (key == GLFW_KEY_T && action == GLFW_PRESS) {
    if (drawData->settings.activeness[std::to_underlying(ProblemType::BTSP_approx)]) {
      cycleBTSPApproxDisplay(drawData->settings);
    }
    if (drawData->settings.activeness[std::to_underlying(ProblemType::BTSPP_approx)]) {
      cycleBTSPPApproxDisplay(drawData->settings);
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