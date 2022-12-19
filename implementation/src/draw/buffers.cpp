#include "draw/buffers.hpp"

void VertexArray::mapBufferToAttribute(const VertexBuffer& vbo, const GLuint shaderProgramID, const char* name) {
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
