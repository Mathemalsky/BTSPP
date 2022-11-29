#include "draw/openglsetup.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

OpenGLHandler::~OpenGLHandler() {
  // delete Opengl data
  glDeleteVertexArrays(1, &pVertexArrayID);
  glDeleteBuffers(1, &pVertexBufferID);
}

void OpenGLHandler::addVertexArray() {
  glGenVertexArrays(1, &pVertexArrayID);
  glBindVertexArray(pVertexArrayID);
}

void OpenGLHandler::addVertexBuffer() {
  glGenBuffers(1, &pVertexBufferID);
  glBindBuffer(GL_ARRAY_BUFFER, pVertexBufferID);
}

void OpenGLHandler::insertAttribute(unsigned int position, std::string name, unsigned int size) {
  const std::vector<VertexAttribute>::iterator it = pVertexAttributes.begin() + position;
  pVertexAttributes.insert(it, VertexAttribute{name, size});
}

void OpenGLHandler::enableAllVertexAttribArrays() {
  unsigned int offset      = 0;
  unsigned int totalLength = 0;
  for (const VertexAttribute& attribute : pVertexAttributes) {
    totalLength += attribute.size;
  }
  for (const VertexAttribute& attribute : pVertexAttributes) {
    const GLint vertexAttrib = glGetAttribLocation(shaderProgramID(), attribute.name.c_str());
    glEnableVertexAttribArray(vertexAttrib);
    glVertexAttribPointer(
      vertexAttrib, attribute.size, GL_FLOAT, GL_FALSE, totalLength * pDataTypesSize,
      (void*) ((long) (offset * pDataTypesSize)));
    offset += attribute.size;
  }
}
