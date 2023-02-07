#include "draw/draw.hpp"

#include <array>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/graph.hpp"

using namespace drawing;

static RGBA_COLOUR operator*(const RGBA_COLOUR& colour, const float fade) {
  return RGBA_COLOUR{colour[0] * fade, colour[1] * fade, colour[2] * fade, colour[3]};
}

static void clearWindow(GLFWwindow* window) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(CLEAR_COLOUR[0], CLEAR_COLOUR[1], CLEAR_COLOUR[2], CLEAR_COLOUR[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

// the vertex buffer object needs to be bound and the attribute vertex_position needs to be enabled
static void drawVerteces(const ShaderProgram& drawCircles) {
  drawCircles.use();  // need to call glUseProgram before setting uniforms
  drawCircles.setUniform("u_steps", 8);
  drawCircles.setUniform("u_radius", VETREX_RADIUS);
  drawCircles.setUniform("u_colour", VERTEX_COLOUR);

  glDrawArrays(GL_POINTS, 0, EUCLIDEAN.numberOfNodes());  // start at index 0
}

static void drawPath(
    const ShaderProgram& drawPathSegments, const ShaderBuffer& shaderBuffer, const std::vector<unsigned int>& order,
    const float thickness, const RGBA_COLOUR& colour) {
  shaderBuffer.bufferSubData(order);
  drawPathSegments.use();
  drawPathSegments.setUniform("u_thickness", thickness);
  drawPathSegments.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawPathSegments.setUniform("u_colour", colour);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6 * (order.size() - PATH_OVERHEAD));
}

static void drawEdge(const ShaderProgram& drawLine, const Edge& e, const float thickness, const RGBA_COLOUR& colour) {
  drawLine.use();
  drawLine.setUniform("u_ends", POINTS_F[2 * e.u], POINTS_F[2 * e.u + 1], POINTS_F[2 * e.v], POINTS_F[2 * e.v + 1]);
  drawLine.setUniform("u_thickness", thickness);
  drawLine.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawLine.setUniform("u_colour", colour);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6);
}

static void drawGraph(
    const ShaderProgram& drawLine, const AdjacencyMatrixGraph<Directionality::Undirected>& graph,
    const RGBA_COLOUR& colour) {
  for (unsigned int k = 0; k < graph.numberOfNodes(); ++k) {
    for (Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>::InnerIterator it(graph.matrix(), k); it; ++it) {
      drawEdge(drawLine, Edge{(size_t) it.row(), (size_t) it.col()}, 5.0f, colour);
    }
  }
}

static void drawOpenEarDecomposition(const ShaderProgram& drawLine, const OpenEarDecomposition& openEarDecomp) {
  for (unsigned int i = 0; i < openEarDecomp.ears.size(); ++i) {
    const std::vector<size_t>& chain = openEarDecomp.ears[i];
    RGBA_COLOUR colour =
        COLOUR[static_cast<unsigned int>(ProblemType::BTSP_approx)] * ((float) i / (openEarDecomp.ears.size() - 1));
    for (unsigned int j = chain.size() - 1; j > 0; --j) {
      drawEdge(drawLine, Edge{chain[j], chain[j - 1]}, 5.0f, colour);
    }
  }
}

void draw(GLFWwindow* window, const ShaderProgramCollection& programs, const Buffers& buffers) {
  clearWindow(window);
  drawVerteces(programs.drawCircles);

  unsigned int typeInt = static_cast<unsigned int>(ProblemType::BTSP_approx);
  if (DRAW_BICONNECTED_GRAPH && ACTIVE[typeInt] && INITIALIZED[typeInt]) {
    drawGraph(programs.drawLine, BTSP_APPROX_RESULT.biconnectedGraph, COLOUR[typeInt]);
  }
  if (DRAW_OPEN_EAR_DECOMPOSITION && ACTIVE[typeInt] && INITIALIZED[typeInt]) {
    drawOpenEarDecomposition(programs.drawLine, BTSP_APPROX_RESULT.openEarDecomposition);
  }
  typeInt = static_cast<unsigned int>(ProblemType::BTSP_exact);
  if (ACTIVE[typeInt] && INITIALIZED[typeInt]) {
    drawPath(programs.drawPathSegments, buffers.tour, ORDER[typeInt], THICKNESS[typeInt], COLOUR[typeInt]);
    drawEdge(programs.drawLine, BTSP_EXACT_RESULT.bottleneckEdge, THICKNESS[typeInt] * 1.75f, COLOUR[typeInt]);
  }
  typeInt = static_cast<unsigned int>(ProblemType::BTSPP_exact);
  if (ACTIVE[typeInt] && INITIALIZED[typeInt]) {
    drawPath(programs.drawPathSegments, buffers.tour, ORDER[typeInt], THICKNESS[typeInt], COLOUR[typeInt]);
    drawEdge(programs.drawLine, BTSPP_EXACT_RESULT.bottleneckEdge, THICKNESS[typeInt] * 1.75f, COLOUR[typeInt]);
  }
  typeInt = static_cast<unsigned int>(ProblemType::TSP_exact);
  if (ACTIVE[typeInt] && INITIALIZED[typeInt]) {
    drawPath(programs.drawPathSegments, buffers.tour, ORDER[typeInt], THICKNESS[typeInt], COLOUR[typeInt]);
  }
}
