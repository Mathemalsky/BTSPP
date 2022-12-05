#include <cassert>
#include <cstdio>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/definitions.hpp"
#include "draw/draw.hpp"
#include "draw/events.hpp"
#include "draw/gui.hpp"
#include "draw/openglerrors.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"
#include "draw/vertexattributes.hpp"

#include "graph/graph.hpp"

#include "utility/datacontainer.hpp"

#include "euclideandistancegraph.hpp"

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

  // generate 20 random vertecis in euclidean plane
  Euclidean euclidean = generateEuclideanDistanceGraph(20);

  graph::POINTS = Data(euclidean.pointer(), euclidean.numberOfNodes());

  vertexArray();
  vertexBuffer();

  ShaderCollection collection;
  ShaderProgram drawCircles = collection.linkCircleDrawProgram();

  VertexAttributes<float> vertexAttributes;
  vertexAttributes.emplaceBack("vertexPosition", 2);
  vertexAttributes.enableAllToShaderProgram(drawCircles.id());

  vertexAttributes.pointsToVertexBufferData(graph::POINTS);
  vertexAttributes.vertexBufferDataToGL();

  drawCircles.use();

  drawCircles.setUniform("u_steps", 6);
  drawCircles.setUniform("u_radius", 0.1f);

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
    handleFastEvents(window, vertexAttributes);

    // draw the content
    draw(window);
    // glUseProgram(programID);

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
