#include "draw/events.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/definitions.hpp"
#include "draw/variables.hpp"
#include "draw/vertexattributes.hpp"

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

void moveNode(GLFWwindow* window, const VertexBuffer& coordinates) {
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
  graph::updatePointsfFromPoints();
  coordinates.bufferSubData(graph::POINTS_F);
}

void handleFastEvents(GLFWwindow* window, const VertexBuffer& coordinates) {
  glfwGetFramebufferSize(window, &mainwindow::WIDTH, &mainwindow::HEIGHT);  // update window size
  if (input::STATE[GLFW_MOUSE_BUTTON_LEFT]) {
    moveNode(window, coordinates);
  }
  glfwGetCursorPos(window, &input::mouse::x, &input::mouse::y);
}
