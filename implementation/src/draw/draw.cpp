#include "draw/draw.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/shader.hpp"
#include "draw/variables.hpp"

static void clearWindow(GLFWwindow* window) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(0.6, 0.6, 0.8, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

static void drawVerteces(const ShaderProgram& drawCircles) {
  drawCircles.use();  // need to call glUseProgram before setting uniforms
  drawCircles.setUniform("u_steps", 8);
  drawCircles.setUniform("u_radius", 0.1f);
  drawCircles.setUniform("u_color", 1.0f, 0.0f, 0.0f, 1.0f);

  glDrawArrays(GL_POINTS, 0, graph::POINTS.size());  // start at index 0
}

static void drawEdges(const ShaderProgram& drawLineSegments) {
  drawLineSegments.use();
  drawLineSegments.setUniform("u_thickness", 5.0f);  // thickness is a float
  drawLineSegments.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawLineSegments.setUniform("u_color", 1.0f, 0.0f, 0.0f, 1.0f);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6 * (graph::POINTS.size()));
}

void draw(GLFWwindow* window, const ShaderProgramCollection& programs) {
  clearWindow(window);
  drawVerteces(programs.drawCircles);
  drawEdges(programs.drawLineSegments);
}
