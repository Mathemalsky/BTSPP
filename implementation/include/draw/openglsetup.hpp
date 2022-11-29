#pragma once

#include <string>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/shader.hpp"

struct VertexAttribute {
  std::string name;
  unsigned int size;
};

class OpenGLHandler {
public:
  OpenGLHandler() = default;

  ~OpenGLHandler();

  void addVertexArray();
  void addVertexBuffer();

  void insertAttribute(unsigned int position, std::string name, unsigned size);
  void pushBackAttribute(std::string name, unsigned int size) {
    pVertexAttributes.push_back(VertexAttribute{name, size});
  }

  void enableAllVertexAttribArrays();

  void linkShaderProgram() { pShaderProgramID = linkShaders(); }

  GLuint shaderProgramID() const { return pShaderProgramID; }
  GLuint vertexArrayID() const { return pVertexArrayID; }
  GLuint vertexBufferID() const { return pVertexBufferID; }

private:
  GLuint pShaderProgramID;
  GLuint pVertexArrayID;
  GLuint pVertexBufferID;

  std::vector<VertexAttribute> pVertexAttributes;
  static constexpr unsigned int pDataTypesSize = sizeof (float);
};
