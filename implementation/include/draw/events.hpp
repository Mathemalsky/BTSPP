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

// DEBUG
#include <iostream>
template <typename Type>
void moveNode(GLFWwindow* window, OpenGLHandler<Type>& openGLHandler) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  const Point2D oldMousePos  = transformCoordinates(input::mouse::x, input::mouse::y);
  const Point2D diffMousePos = transformCoordinates(x, y) - oldMousePos;
  for (Point2D& point : graph::POINTS) {
    if (norm2(point - oldMousePos) < 0.01) {
      std::cerr << point << "-> ";  // DEBUG
      point += diffMousePos;
      std::cerr << point << std::endl;  // DBEUG
    }
  }
  openGLHandler.updatePointsInVertexBufferData(graph::POINTS);
  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}

template <typename Type>
void handleFastEvents(GLFWwindow* window, OpenGLHandler<Type>& openGLHandler) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT]) {
    moveNode(window, openGLHandler);
  }

  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}
