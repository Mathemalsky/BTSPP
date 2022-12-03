#include "draw/shader.hpp"

#include <cassert>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

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

GLuint compileShader(const GLenum shaderType, const GLchar* shaderSource) {
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

void ShaderProgram::link() {
  glLinkProgram(pProgramID);

  int success;
  char infoLog[512];

  glGetProgramiv(pProgramID, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(pProgramID, 512, nullptr, infoLog);
    std::cerr << "ERROR Linking of shaders failed\n" << infoLog << std::endl;
  }
}

void ShaderProgram::setUniform(const char* name, const float value) {
  const GLint location = glGetUniformLocation(pProgramID, name);
  assert(location != -1 && "could not find uniform");
  glUniform1f(location, value);
}

void ShaderProgram::setUniform(const char* name, const int value) {
  const GLint location = glGetUniformLocation(pProgramID, name);
  assert(location != -1 && "could not find uniform");
  glUniform1i(location, value);
}

ShaderCollection::ShaderCollection()
  : pVertexShader(compileShader(GL_VERTEX_SHADER, vertexShaderSource))
  , pGeometryShader(compileShader(GL_GEOMETRY_SHADER, geometryShaderSource))
  , pFragmentShader(compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource)) {
}

ShaderCollection::~ShaderCollection() {
  glDeleteShader(pVertexShader);
  glDeleteShader(pGeometryShader);
  glDeleteShader(pFragmentShader);
}

ShaderProgram ShaderCollection::linkCircleDrawProgram() const {
  ShaderProgram circleProgram;
  circleProgram.attachShader(pVertexShader);
  circleProgram.attachShader(pGeometryShader);
  circleProgram.attachShader(pFragmentShader);
  circleProgram.link();
  return circleProgram;
}

GLuint linkShaders() {
  const GLuint shaderProgram  = glCreateProgram();
  const GLuint vertexShader   = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
  const GLuint geometryShader = compileShader(GL_GEOMETRY_SHADER, geometryShaderSource);
  const GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, geometryShader);
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
  glDeleteShader(geometryShader);

  return shaderProgram;
}
