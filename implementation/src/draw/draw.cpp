#include "draw/draw.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
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
  drawCircles.setUniform("u_radius", drawing::VETREX_RADIUS);
  drawCircles.setUniform("u_color", imguiwindow::VERTEX_COLOR);

  glDrawArrays(GL_POINTS, 0, graph::EUCLIDEAN.numberOfNodes());  // start at index 0
}

static void drawEdges(const ShaderProgram& drawLineSegments, const ProblemType type) {
  drawLineSegments.use();
  drawLineSegments.setUniform("u_thickness", imguiwindow::THICKNESS[(unsigned int) type]);
  drawLineSegments.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawLineSegments.setUniform("u_color", imguiwindow::COLOR[(unsigned int) type]);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6 * graph::EUCLIDEAN.numberOfNodes());
}

void draw(GLFWwindow* window, const ShaderProgramCollection& programs, Buffers& buffers) {
  clearWindow(window);
  drawVerteces(programs.drawCircles);

  for (const ProblemType& type : problemType::PROBLEM_TYPES) {
    if (imguiwindow::ACTIVE[(unsigned int) type]) {
      buffers.tour.bufferSubData(graph::ORDER[(unsigned int) type]);
      drawEdges(programs.drawLineSegments, type);
    }
  }
}
