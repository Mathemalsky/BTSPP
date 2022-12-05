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

void ShaderProgram::link() const {
  glLinkProgram(pProgramID);

  int success;
  char infoLog[512];

  GL_CALL(glGetProgramiv(pProgramID, GL_LINK_STATUS, &success);)
  if (!success) {
    GL_CALL(glGetProgramInfoLog(pProgramID, 512, nullptr, infoLog);)
    std::cerr << "ERROR Linking of shaders failed\n" << infoLog << std::endl;
  }
}

void ShaderProgram::setUniform(const char* name, const float value) {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform1f(location, value);)
}

void ShaderProgram::setUniform(const char* name, const int value) {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform1i(location, value);)
}

ShaderCollection::ShaderCollection()
  : pVertexShader(compileShader(GL_VERTEX_SHADER, vertexShaderSource))
  , pGeometryShader(compileShader(GL_GEOMETRY_SHADER, geometryShaderSource))
  , pFragmentShader(compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource)) {
}

ShaderCollection::~ShaderCollection() {
  GL_CALL(glDeleteShader(pVertexShader);)
  GL_CALL(glDeleteShader(pGeometryShader);)
  GL_CALL(glDeleteShader(pFragmentShader);)
}

ShaderProgram ShaderCollection::linkCircleDrawProgram() const {
  ShaderProgram circleProgram;
  circleProgram.attachShader(pVertexShader);
  circleProgram.attachShader(pGeometryShader);
  circleProgram.attachShader(pFragmentShader);
  circleProgram.link();
  return circleProgram;
}
