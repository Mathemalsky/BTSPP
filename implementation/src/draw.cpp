#include "draw.hpp"

#include <GLFW/glfw3.h>

#include "variables.hpp"

void draw(GLFWwindow* window) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
}
