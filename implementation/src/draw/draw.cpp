#include "draw/draw.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/graph.hpp"

static void clearWindow(GLFWwindow* window) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(0.6, 0.6, 0.8, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

// the vertex buffer object needs to be bound and the attribute vertex_position needs to be enabled
static void drawVerteces(const ShaderProgram& drawCircles) {
  drawCircles.use();  // need to call glUseProgram before setting uniforms
  drawCircles.setUniform("u_steps", 8);
  drawCircles.setUniform("u_radius", drawing::VETREX_RADIUS);
  drawCircles.setUniform("u_colour", drawing::VERTEX_COLOUR);

  glDrawArrays(GL_POINTS, 0, drawing::EUCLIDEAN.numberOfNodes());  // start at index 0
}

static void drawCycle(
    const ShaderProgram& drawCycleSegments, const ShaderBuffer& shaderBuffer, const std::vector<unsigned int>& order,
    const float thickness, const RGBA_COLOUR& colour) {
  shaderBuffer.bufferSubData(order);
  drawCycleSegments.use();
  drawCycleSegments.setUniform("u_thickness", thickness);
  drawCycleSegments.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawCycleSegments.setUniform("u_colour", colour);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6 * (order.size() - 3));
}

static void drawEdge(const ShaderProgram& drawLine, const Edge& e, const float thickness, const RGBA_COLOUR& colour) {
  const std::array<float, 4> coords{
      drawing::POINTS_F[2 * e.u], drawing::POINTS_F[2 * e.u + 1], drawing::POINTS_F[2 * e.v],
      drawing::POINTS_F[2 * e.v + 1]};
  drawLine.use();
  drawLine.setUniform("u_ends", coords);
  drawLine.setUniform("u_thickness", thickness);
  drawLine.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawLine.setUniform("u_colour", colour);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6);
}

void drawGraph(const ShaderProgram& drawLine, const AdjacencyMatrixGraph<Directionality::Undirected>& graph) {
  for (unsigned int k = 0; k < graph.numberOfNodes(); ++k) {
    for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(graph.matrix(), k); it; ++it) {
      drawEdge(drawLine, Edge{(size_t) it.row(), (size_t) it.col()}, 5.0f, RGBA_COLOUR{1.0, 0.0, 1.0, 1.0});
    }
  }
}

void draw(GLFWwindow* window, const ShaderProgramCollection& programs, const Buffers& buffers) {
  clearWindow(window);
  drawVerteces(programs.drawCircles);

  unsigned int typeInt = static_cast<unsigned int>(ProblemType::BTSP_approx);
  if (drawing::ACTIVE[typeInt] && drawing::INITIALIZED[typeInt]) {
    drawGraph(programs.drawLine, drawing::BTSP_APPROX_RESULT.biconnectedGraph);
  }
  typeInt = static_cast<unsigned int>(ProblemType::BTSP_exact);
  if (drawing::ACTIVE[typeInt] && drawing::INITIALIZED[typeInt]) {
    drawCycle(
        programs.drawCycleSegments, buffers.tour, drawing::ORDER[typeInt], drawing::THICKNESS[typeInt],
        drawing::COLOUR[typeInt]);
    drawEdge(
        programs.drawLine, drawing::BTSP_EXACT_RESULT.bottleneckEdge, drawing::THICKNESS[typeInt] * 1.5f,
        drawing::COLOUR[typeInt]);
  }
  typeInt = static_cast<unsigned int>(ProblemType::TSP_exact);
  if (drawing::ACTIVE[typeInt] && drawing::INITIALIZED[typeInt]) {
    drawCycle(
        programs.drawCycleSegments, buffers.tour, drawing::ORDER[typeInt], drawing::THICKNESS[typeInt],
        drawing::COLOUR[typeInt]);
  }
}
