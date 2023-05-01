#include "draw/visualization.hpp"

#include <cstdio>
#include <functional>
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
using namespace std::placeholders;

// error callback function which prints glfw errors in case they arise
static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
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
      std::make_shared<ShaderBuffer>(floatVertices.read());  // copy vertex coords to shader buffer
  std::shared_ptr<ShaderBuffer> tour =
      std::make_shared<ShaderBuffer>(std::vector<unsigned int>(euclidean.numberOfNodes() + 3));  // just allocate memory

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

void visualize(const graph::Euclidean& euclidean) {
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
  auto keyCallbackAdapter = std::bind(keyCallback, _1, _2, _3, _4, _5, drawData);
  glfwSetKeyCallback(window, keyCallbackAdapter);
  glfwSetMouseButtonCallback(window, mouseButtonCallback);

  // setup Dear ImGui
  setUpImgui(window, glsl_version);

  // set initial state of variables for drawing and input
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
    drawImgui(drawData->appearance, drawData->settings);

    // swap the drawings to the displayed frame
    glfwSwapBuffers(window);
  }

  // clean up Dear ImGui
  cleanUpImgui();

  // clean up glfw window
  glfwDestroyWindow(window);
  glfwTerminate();
}
