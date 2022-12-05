#include <cassert>
#include <cstdio>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/definitions.hpp"
#include "draw/draw.hpp"
#include "draw/events.hpp"
#include "draw/gui.hpp"
#include "draw/openglhandler.hpp"
#include "draw/openglerrors.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/graph.hpp"

#include "utility/datacontainer.hpp"

#include "euclideandistancegraph.hpp"

// DEBUG

static constexpr const char vertexShaderSource[] = R"glsl(
  #version 440 core
  in vec2 vertexPosition;

  void main() {
    gl_Position = vec4(vertexPosition, 0.0, 1.0);
  }
)glsl";

static constexpr const char geometryShaderSource[] = R"glsl(
  #version 440 core
  layout(points) in;
  layout(line_strip, max_vertices = 21) out;

  uniform float u_radius;
  uniform int u_steps;

  const float PI = 3.1415926;

  void main() {
    for(int i = 0; i < u_steps +1; ++i) {
      float ang = 2.0 * PI * i / u_steps;

      vec4 offset = vec4(cos(ang) * u_radius/16, sin(ang) * u_radius/9, 0.0, 0.0);
      gl_Position = gl_in[0].gl_Position + offset;
      EmitVertex();
    }

    EndPrimitive();
}
)glsl";

static constexpr const char fragmentShaderSource[] = R"glsl(
  #version 440 core
  out vec4 color;

  void main() {
    color = vec4(1.0, 0.0, 0.0, 1.0);
  }
)glsl";

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

  // malke sure OpenGL errors will be raised in the correct scope
  /*
  glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);
  glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
  int flags;
  glGetIntegerv(GL_CONTEXT_FLAGS, &flags);

  if (flags & GL_CONTEXT_FLAG_DEBUG_BIT)
  {
      glEnable(GL_DEBUG_OUTPUT);
      glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
      glDebugMessageCallback(glDebugOutput, nullptr);
      glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);
  }
  */

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

  /* begin setting up opengl */

  // generate 20 random vertecis in euclidean plane
  Euclidean euclidean = generateEuclideanDistanceGraph(20);

  graph::POINTS = Data(euclidean.pointer(), euclidean.numberOfNodes());

  /*
  OpenGLHandler<float> openGLHandler;
  //openGLHandler.linkShaderProgram();
  // openGLHandler.linkPrograms();

  //std::cerr << openGLHandler.shaderProgramID() << std::endl;
  // openGLHandler.test();
  openGLHandler.addVertexArray();
  openGLHandler.addVertexBuffer();


  // openGLHandler.emplaceBackAttribute("size", 1);
  // openGLHandler.emplaceBackAttribute("steps", 1);

  openGLHandler.linkPrograms();

  openGLHandler.enableAllVertexAttribArrays();
  openGLHandler.emplaceBackAttribute("vertexPosition", 2);

  openGLHandler.pointsToVertexBufferData(graph::POINTS);
  openGLHandler.vertexBufferDataToGL();



  GLint vertexStepsLocation = glGetUniformLocation(openGLHandler.shaderProgramID(), "u_steps");
  assert(vertexStepsLocation != -1 && "could not find uniform");
  std::cerr <<"shaderProgramID " << openGLHandler.shaderProgramID() << std::endl;
  std::cerr <<"vertexStepslocation " << vertexStepsLocation << std::endl;
  glUniform1i(vertexStepsLocation, 4);

  GLint vertexRadiusLocation = glGetUniformLocation(openGLHandler.shaderProgramID(), "u_radius");
  assert(vertexRadiusLocation != -1 && "could not find uniform");
  glUniform1f(vertexRadiusLocation, 0.1f);

  // openGLHandler.test();

  */
  /* end test setting up opengl */

  /* second test */


  vertexArray();
  vertexBuffer();

  ShaderCollection collection;
  ShaderProgram drawCircles = collection.linkCircleDrawProgram();
  std::cerr << "new ShaderProgram: " << drawCircles.id() << std::endl;
  // drawCircles.link();


  VertexAttributes<float> vertexAttributes;
  vertexAttributes.emplaceBack("vertexPosition", 2);
  vertexAttributes.enableAllToShaderProgram(drawCircles.id());

  /*
  // enable vertex attributes
  GL_CALL(const GLint vertexAttrib = glGetAttribLocation(drawCircles.id(), "vertexPosition");)
  GL_CALL(glEnableVertexAttribArray(vertexAttrib);)
  GL_CALL(glVertexAttribPointer(vertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, (void*) (0));)
  */

  vertexAttributes.pointsToVertexBufferData(graph::POINTS);
  vertexAttributes.vertexBufferDataToGL();

  /*
  // copy data to vertexbuffer
  const size_t size       = graph::POINTS.size();
  float* vertexBufferData = new float[size * 2];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * 2]     = graph::POINTS[i].x;  // cannot use memcpy here because we
    vertexBufferData[i * 2 + 1] = graph::POINTS[i].y;  // need to cast double to float
  }


  GL_CALL(glBufferData(GL_ARRAY_BUFFER, 160, vertexBufferData, GL_DYNAMIC_DRAW);)
  */

  drawCircles.use();

  drawCircles.setUniform("u_steps", 6);
  drawCircles.setUniform("u_radius", 0.1f);


  // set uniforms
  /*
  GL_CALL(GLint vertexStepsLocation2 = glGetUniformLocation(drawCircles.id(), "u_steps");)
  assert(vertexStepsLocation2 != -1 && "could not find uniform");
  std::cerr << "vertexStepslocation2 " << vertexStepsLocation2 << std::endl;
  GL_CALL(glUniform1i(vertexStepsLocation2, 6);)

  GL_CALL(GLint vertexRadiusLocation2 = glGetUniformLocation(drawCircles.id(), "u_radius");)
  assert(vertexRadiusLocation2 != -1 && "could not find uniform");
  GL_CALL(glUniform1f(vertexRadiusLocation2, 0.1f);)
  */

  /* second test end*/

  /* basic test */

  /*
  // copy data to vertexbuffer
  const size_t size        = graph::POINTS.size();
  float* vertexBufferData = new float[size * 2];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * 2] = graph::POINTS[i].x;  // cannot use memcpy here because we
    vertexBufferData[i * 2 + 1] = graph::POINTS[i].y;  // need to cast double to float
  }

  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, 160, vertexBufferData, GL_STATIC_DRAW);

  GLuint vao;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  GLuint shaderProgram = glCreateProgram();

  std::cerr << "shaderprogram : " << shaderProgram << std::endl;

  GLuint vertexshader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
  GLuint geometryshader = compileShader(GL_GEOMETRY_SHADER, geometryShaderSource);
  GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);

  glAttachShader(shaderProgram, vertexshader);
  glAttachShader(shaderProgram, geometryshader);
  glAttachShader(shaderProgram, fragmentShader);

  glLinkProgram(shaderProgram);

  const GLint vertexAttrib = glGetAttribLocation(shaderProgram, "vertexPosition");
  std::cerr <<"vertexAttrib " << vertexAttrib << std::endl;
  glEnableVertexAttribArray(vertexAttrib);
  glVertexAttribPointer(vertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

  int success;
  char infoLog[512];

  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
    std::cerr << "ERROR Linking of shaders failed\n" << infoLog << std::endl;
  }

  glDeleteShader(vertexshader);
  glDeleteShader(geometryshader);
  glDeleteShader(fragmentShader);

  glUseProgram(shaderProgram);

  // set uniforms
  GLint vertexStepsLocation2 = glGetUniformLocation(shaderProgram, "u_steps");
  assert(vertexStepsLocation2 != -1 && "could not find uniform");
  std::cerr <<"vertexStepslocation2 " << vertexStepsLocation2 << std::endl;
  glUniform1i(vertexStepsLocation2, 6);

  GLint vertexRadiusLocation2 = glGetUniformLocation(shaderProgram, "u_radius");
  assert(vertexRadiusLocation2 != -1 && "could not find uniform");
  glUniform1f(vertexRadiusLocation2, 0.1f);

  /* basic test end */

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

  // DEBUG
  //delete[] vertexBufferData;

  // clean up Dear ImGui
  cleanUpImgui();

  // clean up glfw window
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
