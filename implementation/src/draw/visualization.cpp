#include "draw/visualization.hpp"

#include <cstdio>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/draw.hpp"
#include "draw/events.hpp"
#include "draw/gui.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "euclideandistancegraph.hpp"
#include "exactsolver.hpp"

// error callback function which prints glfw errors in case they arise
static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void initDrawingVariables() {
  drawing::SHOW_SETTINGS_WINDOW = drawing::INITIAL_SHOW_SETTINGS_WINDOW;
  drawing::ACTIVE               = drawing::INITIAL_ACTIVENESS;
  drawing::COLOUR               = drawing::INITIAL_COLOUR;
  drawing::ORDER_INITIALIZED    = drawing::INITIAL_ACTIVENESS;
  drawing::THICKNESS            = drawing::INITIAL_THICKNESS;
  drawing::VERTEX_COLOUR        = drawing::INITIAL_VERTEX_COLOUR;
}

static const Buffers& setUpBufferMemory(const unsigned int numberOfNodes) {
  // generate random vertices in euclidean plane
  drawing::EUCLIDEAN = std::move(generateEuclideanDistanceGraph(numberOfNodes));
  drawing::updatePointsfFromEuclidean();  // convert to 32 bit floats because opengl isn't capable to deal with 64 bit

  const VertexBuffer& coordinates     = *new VertexBuffer(drawing::POINTS_F, 2);  // components per vertex
  const ShaderBuffer& tourCoordinates = *new ShaderBuffer(drawing::POINTS_F);     // copy vertex coords to shader buffer
  const ShaderBuffer& tour = *new ShaderBuffer(std::vector<unsigned int>(numberOfNodes + 3));  // just allocate memory

  return *new Buffers{coordinates, tour, tourCoordinates};
}

static const VertexArray& bindBufferMemory(const Buffers& buffers, const ShaderProgramCollection& programs) {
  const VertexArray& vao = *new VertexArray;
  vao.bind();
  vao.mapBufferToAttribute(buffers.coordinates, programs.drawCircles.id(), "vertexPosition");
  vao.enable(programs.drawCircles.id(), "vertexPosition");
  vao.bindBufferBase(buffers.tourCoordinates, 0);
  vao.bindBufferBase(buffers.tour, 1);
  return vao;
}

int visualize(const unsigned int numberOfNodes) {
  // set error colback function
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit()) {
    return -1;
  }

  // create window in specified size
  GLFWwindow* window =
      glfwCreateWindow(mainwindow::INITIAL_WIDTH, mainwindow::INITIAL_HEIGHT, mainwindow::NAME, nullptr, nullptr);
  if (window == nullptr)
    return -1;
  glfwMakeContextCurrent(window);

  // initialize glew
  glewExperimental = GL_TRUE;
  if (glewInit()) {
    std::cerr << "failed to initialize glew\n";
    return -1;
  }

  const char* glsl_version = imguiVersionHints();

  const ShaderCollection collection;
  const ShaderProgram drawCircles       = collection.linkCircleDrawProgram();
  const ShaderProgram drawCycleSegments = collection.linkCycleSegementDrawProgram();
  const ShaderProgramCollection programs(drawCircles, drawCycleSegments);

  const Buffers& buffers = setUpBufferMemory(numberOfNodes);
  const VertexArray& vao = bindBufferMemory(buffers, programs);

  // enable vsync
  glfwSwapInterval(1);

  // set callbacks for keyboard and mouse, must be called before Imgui
  glfwSetKeyCallback(window, keyCallback);
  glfwSetMouseButtonCallback(window, mouseButtonCallback);

  // setup Dear ImGui
  setUpImgui(window, glsl_version);

  // set initial state of the settings window
  initDrawingVariables();

  // main loop
  while (!glfwWindowShouldClose(window)) {
    // runs only through the loop if something changed
    glfwPollEvents();

    // handle Events triggert by user input, like keyboard etc.
    handleEvents(window, buffers);

    // draw the content
    draw(window, programs, buffers);

    // draw the gui
    drawImgui();

    // swap the drawings to the displayed frame
    glfwSwapBuffers(window);
  }

  // clean up memory
  delete &vao;
  delete &buffers;

  // clean up Dear ImGui
  cleanUpImgui();

  // clean up glfw window
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
