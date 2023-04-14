#include "draw/events.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/drawdata.hpp"
#include "draw/variables.hpp"

#include "solve/exactsolver.hpp"

/***********************************************************************************************************************
 *                                               fast events
 **********************************************************************************************************************/

static graph::Point2D transformCoordinates(const double x, const double y) {
  return graph::Point2D{2.0 * x / mainwindow::WIDTH - 1.0, -2.0 * y / mainwindow::HEIGHT + 1.0};
}

static void moveNode(GLFWwindow* window, DrawData& drawData) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  drawData.drawgraph.euclidean().vertices()[input::mouse::NODE_IN_MOTION] = transformCoordinates(x, y);
  drawing::updatePointsfFromEuclidean();
  drawData.buffers.coordinates.bufferSubData(drawing::POINTS_F);
  drawData.buffers.tourCoordinates.bufferSubData(drawing::POINTS_F);
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
    drawing::BTSP_APPROX_RESULT                                     = approximation::approximateBTSP(drawing::EUCLIDEAN);
    drawing::updateOrder(drawing::BTSP_APPROX_RESULT.tour, ProblemType::BTSP_approx);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_approx)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_approx)] = false;
    drawing::BTSPP_APPROX_RESULT                                     = approximation::approximateBTSPP(drawing::EUCLIDEAN);
    drawing::updateOrder(drawing::BTSPP_APPROX_RESULT.tour, ProblemType::BTSPP_approx);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSP_exact)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSP_exact)] = false;
    drawing::BTSP_EXACT_RESULT                                     = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSP_exact);
    drawing::updateOrder(drawing::BTSP_EXACT_RESULT.tour, ProblemType::BTSP_exact);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_exact)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::BTSPP_exact)] = false;
    drawing::BTSPP_EXACT_RESULT                                     = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::BTSPP_exact);
    drawing::updateOrder(drawing::BTSPP_EXACT_RESULT.tour, ProblemType::BTSPP_exact);
  }
  if (slowEvents::SOLVE[std::to_underlying(ProblemType::TSP_exact)]) {
    slowEvents::SOLVE[std::to_underlying(ProblemType::TSP_exact)] = false;
    exactsolver::Result res                                       = exactsolver::solve(drawing::EUCLIDEAN, ProblemType::TSP_exact);
    drawing::updateOrder(res.tour, ProblemType::TSP_exact);
  }
}

/***********************************************************************************************************************
 *                                           event handling
 **********************************************************************************************************************/

void handleEvents(GLFWwindow* window, DrawData& drawData) {
  handleFastEvents(window, drawData);
  handleSlowEvents(drawData);
}
