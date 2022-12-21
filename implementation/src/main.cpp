#include <cassert>
#include <cstdio>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/definitions.hpp"
#include "draw/draw.hpp"
#include "draw/events.hpp"
#include "draw/gui.hpp"
#include "draw/openglerrors.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/graph.hpp"

#include "utility/utils.hpp"

#include "euclideandistancegraph.hpp"
#include "exactsolver.hpp"

// error callback function which prints glfw errors in case they arise
static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  // set error colback function
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit()) {
    return -1;
  }

  /* Case distictions needed for imgui.
   * See: https://github.com/ocornut/imgui/blob/master/examples/example_glfw_opengl3/main.cpp */
  // Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
  // GL ES 2.0 + GLSL 100
  const char* glsl_version = "#version 100";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
  // GL 3.2 + GLSL 150
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
  // GL 4.4 + GLSL 440
  const char* glsl_version = "#version 440";
  glfwWindowHint(GLFW_SAMPLES, 4);  // 4x antialiasing
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
#endif

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

  // generate random vertecis in euclidean plane
  graph::EUCLIDEAN = std::move(generateEuclideanDistanceGraph(15));
  graph::TOUR      = solveTSP(graph::EUCLIDEAN);
  graph::initPointsfFromEuclidean();   // convert to 32 bit floats because opengl isn't capable to deal with 64 bit
  graph::initTourDrawCycleFromTour();  // append tour by copy of first 3 entries to make line_segment well defined

  const ShaderCollection collection;
  const ShaderProgram drawCircles      = collection.linkCircleDrawProgram();
  const ShaderProgram drawLineSegments = collection.linkLineSegementDrawProgram();
  const ShaderProgramCollection programs(drawCircles, drawLineSegments);

  Buffers buffers;
  buffers.coordinates.bind();
  buffers.coordinates.bufferData(graph::POINTS_F, 2);   // components per vertex
  buffers.tourCoordinates.bufferData(graph::POINTS_F);  // copy vertex coordinates also into shader buffer
  buffers.tour.bufferData(graph::TOUR_DRAW_CYCLE);      // copy indices of verteces in tour to shader buffer

  VertexArray vao;
  vao.bind();
  vao.mapBufferToAttribute(buffers.coordinates, programs.drawCircles.id(), "vertexPosition");
  vao.enable(programs.drawCircles.id(), "vertexPosition");
  vao.bindBufferBase(buffers.tourCoordinates, 0);
  vao.bindBufferBase(buffers.tour, 1);

  // enable vsync
  glfwSwapInterval(1);

  // setup Dear ImGui
  setUpImgui(window, glsl_version);

  // set initial state of the settings window
  initImGuiWindows();

  // set callbacks for keyboard and mouse
  glfwSetKeyCallback(window, keyCallback);
  glfwSetMouseButtonCallback(window, mouseButtonCallback);

  // set input mode
  glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_TRUE);

  // main loop
  while (!glfwWindowShouldClose(window)) {
    // runs only through the loop if something changed
    glfwPollEvents();

    // handle Events triggert by user input, like keyboard etc.
    handleEvents(window, buffers);

    // draw the content
    draw(window, programs);

    // draw the gui
    drawImgui();

    // swap the drawings to the displayed frame
    glfwSwapBuffers(window);
  }

  // clean up Dear ImGui
  cleanUpImgui();

  // clean up glfw window
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
