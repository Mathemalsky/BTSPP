#include <cstdio>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/definitions.hpp"
#include "draw/draw.hpp"
#include "draw/events.hpp"
#include "draw/gui.hpp"
#include "draw/shader.hpp"

#include "graph/graph.hpp"

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

  /* begin test example */

  const GLuint shaderProgram = linkShaders();

  // declare vertex array
  GLuint VertexArrayID;
  glGenVertexArrays(1, &VertexArrayID);
  glBindVertexArray(VertexArrayID);

  // An array of 3 vectors which represents 3 vertices
  static const GLfloat g_vertex_buffer_data[] = {
    -0.45f, 0.45f, 0.3f, 20.0f, 0.45f, 0.45f, 0.3f, 20.0f, 0.45f, -0.45f, 0.3f, 20.0f, -0.45f, -0.45f, 0.3f, 20.0f,
  };

  // This will identify our vertex buffer
  GLuint vertexbuffer;
  // Generate 1 buffer, put the resulting identifier in vertexbuffer
  glGenBuffers(1, &vertexbuffer);
  // The following commands will talk about our 'vertexbuffer' buffer
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  // Give our vertices to OpenGL.
  glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

  const GLint vertexPositionAttrib = glGetAttribLocation(shaderProgram, "vertexPosition");
  glEnableVertexAttribArray(vertexPositionAttrib);
  glVertexAttribPointer(
    vertexPositionAttrib,  // attribute 0. No particular reason for 0, but must match the layout in the shader.
    2,                     // size
    GL_FLOAT,              // type
    GL_FALSE,              // normalized?
    4 * sizeof(GLfloat),   // stride
    0                      // array buffer offset
  );

  const GLint sizeAttrib = glGetAttribLocation(shaderProgram, "size");
  glEnableVertexAttribArray(sizeAttrib);
  glVertexAttribPointer(sizeAttrib, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (void*) (2 * sizeof(GLfloat)));

  const GLint stepsAttrib = glGetAttribLocation(shaderProgram, "steps");
  glEnableVertexAttribArray(stepsAttrib);
  glVertexAttribPointer(stepsAttrib, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (void*) (3 * sizeof(GLfloat)));

  /* end test example */

  // enable vsync
  glfwSwapInterval(1);

  // setup Dear ImGui
  setUpImgui(window, glsl_version);

  // set initial state of the settings window
  initImGuiWindows();

  // set callbacks for keyboard and scrolling
  glfwSetKeyCallback(window, keyCallback);

  // test drawing
  // Euclidean euclidean = generateEuclideanDistanceGraph(20);
  // DrawComponents components{&euclidean};

  // main loop
  while (!glfwWindowShouldClose(window)) {
    // runs only through the loop if something changed
    glfwPollEvents();

    // handle Events triggert by user input, like keyboard etc.
    handleFastEvents(window);

    // draw the content
    // draw(window, components);

    // DEBUG
    testdraw(window);
    // glUseProgram(programID);

    // draw the gui
    drawImgui();

    // swap the drawings to the displayed frame
    glfwSwapBuffers(window);
  }

  // delete Opengl data
  glDeleteBuffers(1, &vertexbuffer);
  glDeleteVertexArrays(1, &VertexArrayID);

  // clean up Dear ImGui
  cleanUpImgui();

  // clean up glfw window
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
