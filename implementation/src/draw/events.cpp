#include "draw/events.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/definitions.hpp"
#include "draw/variables.hpp"

#include "utility/datacontainer.hpp"

void toggleSettingsWindow() {
  if (imguiwindow::SHOW_SETTINGS_WINDOW == true) {
    imguiwindow::SHOW_SETTINGS_WINDOW = false;
  }
  else {
    imguiwindow::SHOW_SETTINGS_WINDOW = true;
  }
}

void keyCallback(
  [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods) {
  if (key == GLFW_KEY_F3 && action == GLFW_PRESS) {
    toggleSettingsWindow();
  }
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
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

void moveNode(GLFWwindow* window) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);
  x -= input::mouse::x;
  y -= input::mouse::y;
  for (Point2D& point : graph::POINTS) {
    if (norm2(point.x - input::mouse::x, point.y - input::mouse::y) < 0.01) {
      point.x += x;
      point.y += y;
    }
  }
  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}

void handleFastEvents(GLFWwindow* window) {
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT]) {
    moveNode(window);
  }
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);
  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}
