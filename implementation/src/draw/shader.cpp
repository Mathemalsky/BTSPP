#include "draw/shader.hpp"

#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

constexpr const char* vertexShaderSource =
  "#version 440 core\n"
  "layout(location = 0) in vec2 vertexPosition;\n"
  "void main()\n"
  "{\n"
  "   gl_Position = vec4(vertexPosition, 0.0, 1.0);\n"
  "}\0";

constexpr const char* fragmentShaderSource =
  "#version 440 core\n"
  "out vec3 color;\n"
  "void main(){\n"
  "  color = vec3(1,0,0);\n"
  "}\0";

GLuint compileShader(const GLenum shaderType, const char* shaderSource) {
  const GLuint shader = glCreateShader(shaderType);
  glShaderSource(shader, 1, &shaderSource, nullptr);
  glCompileShader(shader);
  int success;
  char infoLog[512];
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

  if (!success) {
    glGetShaderInfoLog(shader, 512, nullptr, infoLog);
    std::cout << "ERROR Compilation of shader failed\n" << infoLog << std::endl;
  }
  return shader;
}

GLuint linkShaders() {
  const GLuint shaderProgram  = glCreateProgram();
  const GLuint vertexShader   = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
  const GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);

  int success;
  char infoLog[512];

  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
    std::cout << "ERROR Linking of shaders failed\n" << infoLog << std::endl;
  }

  glUseProgram(shaderProgram);

  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  return shaderProgram;
}
