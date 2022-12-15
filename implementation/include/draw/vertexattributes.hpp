#pragma once

#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglerrors.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

class VertexBuffer {
public:
  VertexBuffer() { GL_CALL(glGenBuffers(1, &pID);) }
  ~VertexBuffer() { GL_CALL(glDeleteBuffers(1, &pID);) }

  void bind() const { GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pID);) }

  template <typename Type>
  void bufferData(Data<Type>& dat, const GLuint componentsPerVertex);

  template <typename Type>
  void bufferSubData(Data<Type>& dat) const;

  GLenum type() const { return pType; }
  GLuint compPerVertex() const { return pComponentsPerVertex; }

private:
  GLuint pID;
  GLenum pType;
  GLuint pComponentsPerVertex;
};

template <typename Type>
void VertexBuffer::bufferData(Data<Type>& dat, const GLuint componentsPerVertex) {
  this->bind();
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, dat.byteSize(), dat.data(), GL_DYNAMIC_DRAW);)

  pComponentsPerVertex = componentsPerVertex;

  if (std::is_same<Type, float>{}) {
    pType = GL_FLOAT;
  }
  else {
    std::cerr << "[Buffer Data] Type not yet supported.\n";
  }
}

// template needs to be declared in header
template <typename Type>
void VertexBuffer::bufferSubData(Data<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, dat.byteSize(), dat.data());)
}

class ShaderBuffer {
public:
  ShaderBuffer() { GL_CALL(glGenBuffers(1, &pID);) }
  ~ShaderBuffer() { GL_CALL(glDeleteBuffers(1, &pID);) }

  GLuint id() const { return pID; }

  void bind() const { GL_CALL(glBindBuffer(GL_SHADER_STORAGE_BUFFER, pID);) }

  template <typename Type>
  void bufferData(std::vector<Type>& dat) const;

private:
  GLuint pID;
};

template <typename Type>
void ShaderBuffer::bufferData(std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferData(GL_SHADER_STORAGE_BUFFER, dat.size() * sizeof(Type), dat.data(), GL_DYNAMIC_DRAW);)
}

class VertexArray {
public:
  VertexArray() { GL_CALL(glGenVertexArrays(1, &pID);) }
  ~VertexArray() { GL_CALL(glDeleteVertexArrays(1, &pID);) }

  void bind() const { GL_CALL(glBindVertexArray(pID);) }

  void mapBufferToAttribute(const VertexBuffer& vbo, const GLuint shaderProgramID, const char* name);

  void enable(const GLuint shaderProgramID, const char* name) const;

  void bindBufferBase(const ShaderBuffer& shaderBuffer);

private:
  GLuint pID;
};
