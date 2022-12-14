#include "draw/vertexattributes.hpp"

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
