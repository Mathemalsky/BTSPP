#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglerrors.hpp"

class ShaderProgram {
public:
  ShaderProgram() : pProgramID(glCreateProgram()) {}
  ~ShaderProgram() { GL_CALL(glDeleteProgram(pProgramID);) }

  GLuint id() const { return pProgramID; }

  void attachShader(const GLuint shader) { GL_CALL(glAttachShader(pProgramID, shader);) }

  void link() const;

  void setUniform(const char* name, const float value);
  void setUniform(const char* name, const int value);

  void use() { GL_CALL(glUseProgram(pProgramID);) }

private:
  GLuint pProgramID;
};

class ShaderCollection {
public:
  ShaderCollection();
  ~ShaderCollection();

  ShaderProgram linkCircleDrawProgram() const;

  // private:
  const GLuint pVertexShader;
  const GLuint pGeometryShader;
  const GLuint pFragmentShader;
};

GLuint linkShaders();
GLuint compileShader(const GLenum shaderType, const GLchar* shaderSource);
