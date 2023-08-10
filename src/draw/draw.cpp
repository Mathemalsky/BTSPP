/*
 * BTSPP is a tool to solve, approximate and draw instances of BTSVPP,
 * BTSPP, BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "draw/draw.hpp"

#include <array>
#include <memory>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/drawdata.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

// graph library
#include "graph.hpp"

namespace drawing {
static RGBA_COLOUR operator*(const RGBA_COLOUR& colour, const float fade) {
  return RGBA_COLOUR{colour[0] * fade, colour[1] * fade, colour[2] * fade, colour[3]};
}

static void clearWindow(GLFWwindow* window, const RGBA_COLOUR& clearColour) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(clearColour[0], clearColour[1], clearColour[2], clearColour[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

// the vertex buffer object needs to be bound and the attribute vertex_position needs to be enabled
static void drawVertices(const ShaderProgram& drawCircles, const size_t numberOfVertices, const RGBA_COLOUR& vertexColour) {
  drawCircles.use();  // need to call glUseProgram before setting uniforms
  drawCircles.setUniform("u_steps", 8);
  drawCircles.setUniform("u_radius", VETREX_RADIUS);
  drawCircles.setUniform("u_colour", vertexColour);

  glDrawArrays(GL_POINTS, 0, numberOfVertices);  // start at index 0
}

static void drawPath(const ShaderProgram& drawPathSegments,
                     const std::shared_ptr<ShaderBuffer> shaderBuffer,
                     const std::vector<uint32_t>& order,
                     const float thickness,
                     const RGBA_COLOUR& colour) {
  shaderBuffer->bufferSubData(order);
  drawPathSegments.use();
  drawPathSegments.setUniform("u_thickness", thickness);
  drawPathSegments.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawPathSegments.setUniform("u_colour", colour);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6 * (order.size() - PATH_OVERHEAD));
}

static void drawEdge(const ShaderProgram& drawLine,
                     const FloatVertices& floatVertices,
                     const graph::Edge& e,
                     const float thickness,
                     const RGBA_COLOUR& colour) {
  drawLine.use();
  drawLine.setUniform("u_ends", floatVertices.xCoord(e.u), floatVertices.yCoord(e.u), floatVertices.xCoord(e.v), floatVertices.yCoord(e.v));
  drawLine.setUniform("u_thickness", thickness);
  drawLine.setUniform("u_resolution", mainwindow::WIDTH, mainwindow::HEIGHT);
  drawLine.setUniform("u_colour", colour);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDrawArrays(GL_TRIANGLES, 0, 6);
}

template <typename G>
static void drawGraph(const ShaderProgram& drawLine,
                      const FloatVertices& floatVertices,
                      const G& graph,
                      const float thickness,
                      const RGBA_COLOUR& colour) {
  for (const graph::Edge& e : graph.edges()) {
    drawEdge(drawLine, floatVertices, e, thickness, colour);
  }
}

static void drawOpenEarDecomposition(const ShaderProgram& drawLine,
                                     const std::shared_ptr<DrawData> drawData,
                                     const graph::EarDecomposition& openEarDecomp) {
  for (size_t i = 0; i < openEarDecomp.ears.size(); ++i) {
    const std::vector<size_t>& chain = openEarDecomp.ears[i];
    RGBA_COLOUR colour =
        drawData->appearance.colour[std::to_underlying(ProblemType::BTSP_approx)] * ((float) i / (openEarDecomp.ears.size() - 1));
    for (size_t j = chain.size() - 1; j > 0; --j) {
      drawEdge(drawLine,
               drawData->floatVertices,
               graph::Edge{chain[j], chain[j - 1]},
               drawData->appearance.thickness[std::to_underlying(ProblemType::BTSP_approx)],
               colour);
    }
  }
}

void draw(GLFWwindow* window, const ShaderProgramCollection& programs, const std::shared_ptr<DrawData> drawData) {
  clearWindow(window, drawData->appearance.clearColour);
  drawVertices(programs.drawCircles, EUCLIDEAN.numberOfNodes(), drawData->appearance.vertexColour);

  // BTSP approx
  unsigned int typeInt = std::to_underlying(ProblemType::BTSP_approx);
  if (BTSP_DRAW_BICONNECTED_GRAPH && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSP_approx)) {
    drawGraph(programs.drawLine,
              drawData->floatVertices,
              drawData->results.BTSP_APPROX_RESULT.biconnectedGraph,
              drawData->appearance.thickness[typeInt],
              drawData->appearance.colour[typeInt]);
  }
  if (BTSP_DRAW_OPEN_EAR_DECOMPOSITION && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSP_approx)) {
    drawOpenEarDecomposition(programs.drawLine, drawData, drawData->results.BTSP_APPROX_RESULT.openEarDecomposition);
  }
  if (BTSP_DRAW_HAMILTON_CYCLE && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSP_approx)) {
    drawPath(programs.drawPathSegments,
             drawData->buffers.tour,
             drawData->vertexOrder[ProblemType::BTSP_approx],
             drawData->appearance.thickness[typeInt],
             drawData->appearance.colour[typeInt]);
    drawEdge(programs.drawLine,
             drawData->floatVertices,
             drawData->results.BTSP_APPROX_RESULT.bottleneckEdge,
             drawData->appearance.thickness[typeInt] * 1.75f,
             drawData->appearance.colour[typeInt]);
  }

  // BTSPP approx
  typeInt = std::to_underlying(ProblemType::BTSPP_approx);
  if (BTSPP_DRAW_BICONNECTED_GRAPH && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSPP_approx)) {
    drawGraph(programs.drawLine,
              drawData->floatVertices,
              drawData->results.BTSPP_APPROX_RESULT.biconnectedGraph,
              drawData->appearance.thickness[typeInt],
              drawData->appearance.colour[typeInt]);
  }
  if (BTSPP_DRAW_HAMILTON_PATH && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSPP_approx)) {
    drawPath(programs.drawPathSegments,
             drawData->buffers.tour,
             drawData->vertexOrder[ProblemType::BTSPP_approx],
             drawData->appearance.thickness[typeInt],
             drawData->appearance.colour[typeInt]);
    drawEdge(programs.drawLine,
             drawData->floatVertices,
             drawData->results.BTSPP_APPROX_RESULT.bottleneckEdge,
             drawData->appearance.thickness[typeInt] * 1.75f,
             drawData->appearance.colour[typeInt]);
  }

  // BTSVPP approx
  typeInt = std::to_underlying(ProblemType::BTSVPP_approx);
  if (BTSVPP_DRAW_BICONNECTED_GRAPH && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSVPP_approx)) {
    drawGraph(programs.drawLine,
              drawData->floatVertices,
              drawData->results.BTSVPP_APPROX_RESULT.biconnectedGraph,
              drawData->appearance.thickness[typeInt],
              drawData->appearance.colour[typeInt]);
  }
  if (BTSVPP_DRAW_HAMILTON_PATH && ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSVPP_approx)) {
    drawPath(programs.drawPathSegments,
             drawData->buffers.tour,
             drawData->vertexOrder[ProblemType::BTSVPP_approx],
             drawData->appearance.thickness[typeInt],
             drawData->appearance.colour[typeInt]);
    drawEdge(programs.drawLine,
             drawData->floatVertices,
             drawData->results.BTSVPP_APPROX_RESULT.bottleneckEdge,
             drawData->appearance.thickness[typeInt] * 1.75f,
             drawData->appearance.colour[typeInt]);
  }

  // BTSP exact
  typeInt = std::to_underlying(ProblemType::BTSP_exact);
  if (ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSP_exact)) {
    drawPath(programs.drawPathSegments,
             drawData->buffers.tour,
             drawData->vertexOrder[ProblemType::BTSP_exact],
             drawData->appearance.thickness[typeInt],
             drawData->appearance.colour[typeInt]);
    drawEdge(programs.drawLine,
             drawData->floatVertices,
             drawData->results.BTSP_EXACT_RESULT.bottleneckEdge,
             drawData->appearance.thickness[typeInt] * 1.75f,
             drawData->appearance.colour[typeInt]);
  }

  // BTSPP exact
  typeInt = std::to_underlying(ProblemType::BTSPP_exact);
  if (ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::BTSPP_exact)) {
    drawPath(programs.drawPathSegments,
             drawData->buffers.tour,
             drawData->vertexOrder[ProblemType::BTSPP_exact],
             drawData->appearance.thickness[typeInt],
             drawData->appearance.colour[typeInt]);
    drawEdge(programs.drawLine,
             drawData->floatVertices,
             drawData->results.BTSPP_EXACT_RESULT.bottleneckEdge,
             drawData->appearance.thickness[typeInt] * 1.75f,
             drawData->appearance.colour[typeInt]);
  }

  // TSP exact
  typeInt = std::to_underlying(ProblemType::TSP_exact);
  if (ACTIVE[typeInt] && drawData->vertexOrder.initialized(ProblemType::TSP_exact)) {
    drawPath(programs.drawPathSegments,
             drawData->buffers.tour,
             drawData->vertexOrder[ProblemType::TSP_exact],
             drawData->appearance.thickness[typeInt],
             drawData->appearance.colour[typeInt]);
  }
}
}  // namespace drawing
