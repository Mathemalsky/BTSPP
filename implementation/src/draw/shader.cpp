#include "draw/shader.hpp"

#include <array>
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

static constexpr const char pathVertexShaderSource[] = R"glsl(
  #version 440 core
  layout(std430, binding = 0) buffer lineVertex
  {
     vec2 vertex[];
  };

  layout(std430, binding = 1) buffer lineIndex
  {
    uint index[];
  };

  uniform float u_thickness;
  uniform vec2 u_resolution;

  void main() {
    int line_segment     = gl_VertexID / 6;
    int triangle_vertex  = gl_VertexID % 6;

    vec2 line_direction = normalize(vertex[index[line_segment + 2]] - vertex[index[line_segment + 1]]);
    vec2 line_perpendicular = vec2(-line_direction.y, line_direction.x);

    vec2 pos;
    if(triangle_vertex == 0 || triangle_vertex == 2 || triangle_vertex == 5) {
      vec2 prev_direction = normalize(vertex[index[line_segment + 1]] - vertex[index[line_segment]]);
      vec2 prev_perpendicular = vec2(-prev_direction.y, prev_direction.x);

      vec2 corner_direction = prev_perpendicular + line_perpendicular;
      vec2 offset = u_thickness / dot(corner_direction, line_perpendicular) * corner_direction / u_resolution;

      pos = vertex[index[line_segment + 1]];
      if(triangle_vertex == 0) {
        pos -= offset;
      }
      else {
        pos += offset;
      }
    }
    else {
      vec2 succ_direction = normalize(vertex[index[line_segment + 3]] - vertex[index[line_segment + 2]]);
      vec2 succ_perpendicular = vec2(-succ_direction.y, succ_direction.x);

      vec2 corner_direction = succ_perpendicular + line_perpendicular;
      vec2 offset = u_thickness / dot(corner_direction, line_perpendicular) * corner_direction / u_resolution;

      pos = vertex[index[line_segment + 2]];
      if(triangle_vertex == 4) {
        pos += offset;
      }
      else {
        pos -= offset;
      }
    }
    gl_Position = vec4(pos, 0.0, 1.0);
  }
)glsl";

static constexpr const char lineVertexShaderSource[] = R"glsl(
  #version 440 core

  uniform vec4 u_ends;
  uniform float u_thickness;
  uniform vec2 u_resolution;

  void main() {
    int triangle_vertex  = gl_VertexID % 6;
    vec2 begin = u_ends.xy;
    vec2 end   = u_ends.zw;
    vec2 direction = normalize(end - begin);
    vec2 perpendicular = vec2(-direction.y, direction.x);
    vec2 offset = u_thickness * perpendicular / u_resolution;

    vec2 pos;
    if(triangle_vertex == 0 || triangle_vertex == 2 || triangle_vertex == 5) {
      pos = begin;
      if(triangle_vertex == 0) {
        pos -= offset;
      }
      else {
        pos += offset;
      }
    }
    else {
      pos = end;
      if(triangle_vertex == 4) {
        pos += offset;
      }
      else {
        pos -= offset;
      }
    }
    gl_Position = vec4(pos, 0.0, 1.0);
  }
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

      vec4 offset = vec4(cos(ang) * u_radius, sin(ang) * u_radius, 0.0, 0.0);
      gl_Position = gl_in[0].gl_Position + offset;
      EmitVertex();
    }
    EndPrimitive();
  }
)glsl";

static constexpr const char fragmentShaderSource[] = R"glsl(
  #version 440 core
  layout (location = 0) out vec4 fragColor;

  uniform vec4 u_colour;

  void main() {
    fragColor = u_colour;
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

void ShaderProgram::setUniform(const char* name, const float val1, const float val2) const {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform2f(location, val1, val2);)
}

void ShaderProgram::setUniform(
    const char* name, const float val1, const float val2, const float val3, const float val4) const {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform4f(location, val1, val2, val3, val4);)
}

void ShaderProgram::setUniform(const char* name, const std::array<float, 4>& value) const {
  GL_CALL(const GLint location = glGetUniformLocation(pProgramID, name);)
  assert(location != -1 && "could not find uniform");
  GL_CALL(glUniform4f(location, value[0], value[1], value[2], value[3]);)
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
    std::cout << "ERROR Compilation of " << shaderSource << " failed\n" << infoLog << std::endl;
  }
  return shader;
}

ShaderCollection::ShaderCollection() :
  pVertexShader(compileShader(GL_VERTEX_SHADER, vertexShaderSource)),
  pCircleShader(compileShader(GL_GEOMETRY_SHADER, circleShaderSource)),
  pPathVertexShader(compileShader(GL_VERTEX_SHADER, pathVertexShaderSource)),
  pFragmentShader(compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource)),
  pLineVertexShader(compileShader(GL_VERTEX_SHADER, lineVertexShaderSource)) {
}

ShaderCollection::~ShaderCollection() {
  GL_CALL(glDeleteShader(pVertexShader);)
  GL_CALL(glDeleteShader(pCircleShader);)
  GL_CALL(glDeleteShader(pPathVertexShader);)
  GL_CALL(glDeleteShader(pFragmentShader);)
  GL_CALL(glDeleteShader(pLineVertexShader);)
}

ShaderProgram ShaderCollection::linkCircleDrawProgram() const {
  const ShaderProgram circleProgram;
  circleProgram.attachShader(pVertexShader);
  circleProgram.attachShader(pCircleShader);
  circleProgram.attachShader(pFragmentShader);
  circleProgram.link();
  return circleProgram;
}

ShaderProgram ShaderCollection::linkPathSegementDrawProgram() const {
  const ShaderProgram pathSegmentProgram;
  pathSegmentProgram.attachShader(pPathVertexShader);
  pathSegmentProgram.attachShader(pFragmentShader);
  pathSegmentProgram.link();
  return pathSegmentProgram;
}

ShaderProgram ShaderCollection::linkLineDrawProgram() const {
  const ShaderProgram lineProgram;
  lineProgram.attachShader(pLineVertexShader);
  lineProgram.attachShader(pFragmentShader);
  lineProgram.link();
  return lineProgram;
}
