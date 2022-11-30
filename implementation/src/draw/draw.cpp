#include "draw/draw.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

static void clearWindow(GLFWwindow* window) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(0.6, 0.6, 0.8, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

static void drawCircles() {
  glDrawArrays(GL_POINTS, 0, 20);  // MAGIC NUMBER
}

void draw(GLFWwindow* window) {
  clearWindow(window);
  drawCircles();
}
