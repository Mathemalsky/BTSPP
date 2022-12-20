#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglerrors.hpp"

class ShaderProgram {
public:
  ShaderProgram() : pProgramID(glCreateProgram()) {}
  ~ShaderProgram() { GL_CALL(glDeleteProgram(pProgramID);) }

  GLuint id() const { return pProgramID; }

  void attachShader(const GLuint shader) const { GL_CALL(glAttachShader(pProgramID, shader);) }

  void link() const;

  void setUniform(const char* name, const double value) const;
  void setUniform(const char* name, const int value) const;
  void setUniform(const char* name, const double val1, const double val2) const;
  void setUniform(const char* name, const double val1, const double val2, const double val3, const double val4) const;

  void use() const { GL_CALL(glUseProgram(pProgramID);) }

private:
  const GLuint pProgramID;
};

class ShaderCollection {
public:
  ShaderCollection();
  ~ShaderCollection();

  ShaderProgram linkCircleDrawProgram() const;
  ShaderProgram linkLineSegementDrawProgram() const;

private:
  const GLuint pVertexShader;
  const GLuint pCircleShader;
  const GLuint pFragmentShader;
  const GLuint pLineVertexShader;
};

struct ShaderProgramCollection {
  ShaderProgramCollection(const ShaderProgram& drawCircles, const ShaderProgram& drawLineSegments)
    : drawCircles(drawCircles), drawLineSegments(drawLineSegments) {}
  const ShaderProgram drawCircles;
  const ShaderProgram drawLineSegments;
};
