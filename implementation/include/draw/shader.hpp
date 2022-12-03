#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

class ShaderProgram {
public:
  ShaderProgram() { pProgramID = glCreateProgram(); }
  ~ShaderProgram() { glDeleteProgram(pProgramID); }

  GLuint id() const { return pProgramID; }

  void attachShader(const GLuint shader) { glAttachShader(pProgramID, shader); }

  void link();

  void setUniform(const char* name, const float value);
  void setUniform(const char* name, const int value);

  void use() { glUseProgram(pProgramID); }

private:
  GLuint pProgramID;
};

class ShaderCollection {
public:
  ShaderCollection();
  ~ShaderCollection();

  ShaderProgram linkCircleDrawProgram() const;

private:
  const GLuint pVertexShader;
  const GLuint pGeometryShader;
  const GLuint pFragmentShader;
};

GLuint linkShaders();
