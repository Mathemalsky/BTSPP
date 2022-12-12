#include "draw/shader.hpp"

#include <cassert>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglerrors.hpp"

static constexpr const char vertexShaderSource[] = R"glsl(
  #version 440 core
  in vec2 vertexPosition;

  void main() {
    gl_Position = vec4(vertexPosition, 0.0, 1.0);
  }
)glsl";

static constexpr const char lineVertexShader[] = R"glsl(
#version 440 core

layout(std430, binding = 0) buffer LineVertex
{
   vec4 vertex[];
};

uniform float u_thickness;

)glsl";

static constexpr const char circleShaderSource[] = R"glsl(
  #version 440 core
  layout(points) in;
  layout(line_strip, max_vertices = 21) out;

  uniform float u_radius;
  uniform int u_steps;

  const float PI = 3.1415926;

  void main() {
    for(int i = 0; i < u_steps + 1; ++i) {
      float ang = 2.0 * PI * i / u_steps;

      vec4 offset = vec4(cos(ang) * u_radius / 16.0, sin(ang) * u_radius / 9.0, 0.0, 0.0);
      gl_Position = gl_in[0].gl_Position + offset;
      EmitVertex();
    }

    EndPrimitive();
}
)glsl";

static constexpr const char fragmentShaderSource[] = R"glsl(
  #version 440 core
  layout (location = 0) out vec4 fragColor;

  uniform vec4 u_color;

  void main() {
    fragColor = u_color;
  }
)glsl";

void ShaderProgram::link() const {
  GL_CALL(glLinkProgram(pProgramID);)

  int success;
  char infoLog[512];

  GL_CALL(glGetProgramiv(pProgramID, GL_LINK_STATUS, &success);)
  if (!success) {
    GL_CALL(glGetProgramInfoLog(pProgramID, 512, nullptr, infoLog);)
    std::cerr << "ERROR Linking of shaders failed\n" << infoLog << std::endl;
  }
}

void ShaderProgram::setUniform(const char* name, const float value) const {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform1f(location, value);)
}

void ShaderProgram::setUniform(const char* name, const int value) const {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform1i(location, value);)
}

void ShaderProgram::setUniform(
  const char* name, const float val1, const float val2, const float val3, const float val4) const {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform4f(location, val1, val2, val3, val4);)
}

static GLuint compileShader(const GLenum shaderType, const GLchar* shaderSource) {
  GL_CALL(const GLuint shader = glCreateShader(shaderType);)
  GL_CALL(glShaderSource(shader, 1, &shaderSource, nullptr);)
  GL_CALL(glCompileShader(shader);)
  int success;
  char infoLog[512];
  GL_CALL(glGetShaderiv(shader, GL_COMPILE_STATUS, &success);)

  if (!success) {
    GL_CALL(glGetShaderInfoLog(shader, 512, nullptr, infoLog);)
    std::cout << "ERROR Compilation of shader failed\n" << infoLog << std::endl;
  }
  return shader;
}

ShaderCollection::ShaderCollection()
  : pVertexShader(compileShader(GL_VERTEX_SHADER, vertexShaderSource))
  , pCircleShader(compileShader(GL_GEOMETRY_SHADER, circleShaderSource))
  , pFragmentShader(compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource)) {
}

ShaderCollection::~ShaderCollection() {
  GL_CALL(glDeleteShader(pVertexShader);)
  GL_CALL(glDeleteShader(pCircleShader);)
  GL_CALL(glDeleteShader(pFragmentShader);)
}

ShaderProgram ShaderCollection::linkCircleDrawProgram() const {
  const ShaderProgram circleProgram;
  circleProgram.attachShader(pVertexShader);
  circleProgram.attachShader(pCircleShader);
  circleProgram.attachShader(pFragmentShader);
  circleProgram.link();
  return circleProgram;
}
