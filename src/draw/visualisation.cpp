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
#include "draw/visualisation.hpp"

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/draw.hpp"
#include "draw/drawdata.hpp"
#include "draw/events.hpp"
#include "draw/gui.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "solve/exactsolver.hpp"

using namespace drawing;

// error callback function which prints glfw errors in case they arise
static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void initDrawingVariables() {
  SHOW_DEBUG_WINDOW                = INITIAL_SHOW_DEBUG_WINDOW;
  SHOW_SETTINGS_WINDOW             = INITIAL_SHOW_SETTINGS_WINDOW;
  ACTIVE                           = INITIAL_ACTIVENESS;
  BTSP_DRAW_BICONNECTED_GRAPH      = INITIAL_BTSP_DRAW_BICONNECTED_GRAPH;
  BTSP_DRAW_OPEN_EAR_DECOMPOSITION = INITIAL_BTSP_DRAW_OPEN_EAR_DECOMPOSITION;
  BTSP_DRAW_HAMILTON_CYCLE         = INITIAL_BTSP_DRAW_HAMILTON_CYCLE;
  BTSPP_DRAW_BICONNECTED_GRAPH     = INITIAL_BTSPP_DRAW_BICONNECTED_GRAPH;
  BTSPP_DRAW_HAMILTON_PATH         = INITIAL_BTSPP_DRAW_HAMILTON_PATH;
  BTSVPP_DRAW_BICONNECTED_GRAPH    = INITIAL_BTSVPP_DRAW_BICONNECTED_GRAPH;
  BTSVPP_DRAW_HAMILTON_PATH        = INITIAL_BTSVPP_DRAW_HAMILTON_PATH;
}

static void initInputVariables() {
  input::mouse::NODE_IN_MOTION = input::mouse::INITIAL_NODE_IN_MOTION;
}

static DrawData setUpBufferMemory(const graph::Euclidean& euclidean) {
  drawing::EUCLIDEAN = euclidean;
  FloatVertices floatVertices;
  floatVertices.updatePointsfFromEuclidean(drawing::EUCLIDEAN);

  std::shared_ptr<VertexBuffer> coordinates = std::make_shared<VertexBuffer>(floatVertices.read(), 2);  // components per vertex
  std::shared_ptr<ShaderBuffer> tourCoordinates =
      std::make_shared<ShaderBuffer>(floatVertices.read());                                  // copy vertex coords to shader buffer
  std::shared_ptr<ShaderBuffer> tour =
      std::make_shared<ShaderBuffer>(std::vector<uint32_t>(euclidean.numberOfNodes() + 3));  // just allocate memory

  return DrawData(Buffers{coordinates, tour, tourCoordinates}, floatVertices);
}

static std::unique_ptr<VertexArray> bindBufferMemory(const Buffers& buffers, const ShaderProgramCollection& programs) {
  std::unique_ptr<VertexArray> vao = std::make_unique<VertexArray>();
  vao->bind();
  vao->mapBufferToAttribute(buffers.coordinates, programs.drawCircles.id(), "vertexPosition");
  vao->enable(programs.drawCircles.id(), "vertexPosition");
  vao->bindBufferBase(buffers.tourCoordinates, 0);
  vao->bindBufferBase(buffers.tour, 1);
  return vao;
}

void visualise(const graph::Euclidean& euclidean) {
  // set error colback function
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit()) {
    throw std::runtime_error("[GLFW] Failed to initialize glfw!");
  }

  // create window in specified size
  GLFWwindow* window = glfwCreateWindow(mainwindow::INITIAL_WIDTH, mainwindow::INITIAL_HEIGHT, mainwindow::NAME, nullptr, nullptr);
  if (window == nullptr) {
    throw std::runtime_error("[GLFW] Failed to create window!");
  }
  glfwMakeContextCurrent(window);

  // initialize glew
  glewExperimental = GL_TRUE;
  if (glewInit()) {
    throw std::runtime_error(" [GLEW] Failed to initialize glew!");
  }

  const char* glsl_version = imguiVersionHints();

  const ShaderCollection collection;
  const ShaderProgram drawCircles      = collection.linkCircleDrawProgram();
  const ShaderProgram drawPathSegments = collection.linkPathSegementDrawProgram();
  const ShaderProgram drawLineProgram  = collection.linkLineDrawProgram();
  const ShaderProgramCollection programs(drawCircles, drawPathSegments, drawLineProgram);

  std::shared_ptr<DrawData> drawData = std::make_shared<DrawData>(setUpBufferMemory(euclidean));
  std::unique_ptr<VertexArray> vao   = bindBufferMemory(drawData->buffers, programs);

  // enable vsync
  glfwSwapInterval(1);

  // set callbacks for keyboard and mouse, must be called before Imgui
  glfwSetKeyCallback(window, keyCallback);
  glfwSetMouseButtonCallback(window, mouseButtonCallback);

  // setup Dear ImGui
  setUpImgui(window, glsl_version);

  // set initial state of variables for drawing and input
  initDrawingVariables();
  initInputVariables();

  // main loop
  while (!glfwWindowShouldClose(window)) {
    // runs only through the loop if something changed
    glfwPollEvents();

    // handle Events triggert by user input, like keyboard etc.
    handleEvents(window, drawData);

    // draw the content
    draw(window, programs, drawData);

    // draw the gui
    drawImgui(drawData->appearance);

    // swap the drawings to the displayed frame
    glfwSwapBuffers(window);
  }

  // clean up Dear ImGui
  cleanUpImgui();

  // clean up glfw window
  glfwDestroyWindow(window);
  glfwTerminate();
}
