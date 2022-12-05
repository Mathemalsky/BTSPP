#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglhandler.hpp"
#include "draw/variables.hpp"

#include "graph/geometry.hpp"

void keyCallback(
  [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);

inline Point2D transformCoordinates(const double x, const double y) {
  return Point2D{2.0 * x / mainwindow::WIDTH - 1.0, -2.0 * y / mainwindow::HEIGHT + 1.0};
}

/***********************************************************************************************************************
 *                          implementation of templates needs to be included in header
 **********************************************************************************************************************/

template <typename Type>
void moveNode(GLFWwindow* window, VertexAttributes<Type>& vertexAttributes) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  const Point2D oldMousePos  = transformCoordinates(input::mouse::x, input::mouse::y);
  const Point2D diffMousePos = transformCoordinates(x, y) - oldMousePos;
  for (Point2D& point : graph::POINTS) {
    if (norm2(point - oldMousePos) < 0.01) {  // MAGIC NUMBER
      point += diffMousePos;
      break;  // move only one point at one at a time
    }
  }
  vertexAttributes.updatePointsInVertexBufferData(graph::POINTS);
}

template <typename Type>
void handleFastEvents(GLFWwindow* window, VertexAttributes<Type>& vertexAttributes) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT]) {
    moveNode(window, vertexAttributes);
  }
  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}
