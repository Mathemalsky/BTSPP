#include "draw/buffers.hpp"

template <>
VertexBuffer::VertexBuffer(const std::vector<float>& dat, const GLuint componentsPerVertex) :
  pComponentsPerVertex(componentsPerVertex), pType(GL_FLOAT) {
  GL_CALL(glGenBuffers(1, &pID);)
  this->bind();
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, bytes_of(dat), dat.data(), GL_DYNAMIC_DRAW);)
}

void VertexArray::mapBufferToAttribute(const VertexBuffer& vbo, const GLuint shaderProgramID, const char* name) const {
  this->bind();
  vbo.bind();
  GL_CALL(const GLint vertexAttribLocation = glGetAttribLocation(shaderProgramID, name);)
  GL_CALL(glVertexAttribPointer(vertexAttribLocation, vbo.compPerVertex(), vbo.type(), GL_FALSE, 0, nullptr);)
}

void VertexArray::enable(const GLuint shaderProgramID, const char* name) const {
  this->bind();
  GL_CALL(const GLint vertexAttribLocation = glGetAttribLocation(shaderProgramID, name);)
  GL_CALL(glEnableVertexAttribArray(vertexAttribLocation);)
}

void VertexArray::bindBufferBase(const ShaderBuffer& shaderBuffer, const GLuint bindingPoint) const {
  shaderBuffer.bind();
  this->bind();
  GL_CALL(glBindBufferBase(GL_SHADER_STORAGE_BUFFER, bindingPoint, shaderBuffer.id());)
}
