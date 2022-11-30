#pragma once

#include <string>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/definitions.hpp"
#include "draw/shader.hpp"

#include "graph/geometry.hpp"

struct VertexAttribute {
  std::string name;
  unsigned int size;
};

class OpenGLHandler {
public:
  OpenGLHandler() : pVertAttrTotLen(0) {}

  ~OpenGLHandler();

  void addVertexArray();
  void addVertexBuffer();

  void insertAttribute(unsigned int position, std::string name, unsigned size);
  void pushBackAttribute(std::string name, unsigned int size);

  void enableAllVertexAttribArrays();

  void linkShaderProgram() { pShaderProgramID = linkShaders(); }

  GLuint shaderProgramID() const { return pShaderProgramID; }
  GLuint vertexArrayID() const { return pVertexArrayID; }
  GLuint vertexBufferID() const { return pVertexBufferID; }

  unsigned int vertAttrTotLen() const { return pVertAttrTotLen; }

  float* pointsToVertexBufferData(const Data<Point2D> points);

private:
  GLuint pShaderProgramID;
  GLuint pVertexArrayID;
  GLuint pVertexBufferID;

  std::vector<VertexAttribute> pVertexAttributes;
  unsigned int pVertAttrTotLen;
  static constexpr unsigned int pDataTypesSize = sizeof(float);
};
