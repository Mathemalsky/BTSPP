#include "draw/openglhandler.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "graph/graph.hpp"

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
  pVertAttrTotLen += size;
}

void OpenGLHandler::pushBackAttribute(std::string name, unsigned int size) {
  pVertexAttributes.push_back(VertexAttribute{name, size});
  pVertAttrTotLen += size;
}

void OpenGLHandler::enableAllVertexAttribArrays() {
  unsigned int offset = 0;
  for (const VertexAttribute& attribute : pVertexAttributes) {
    const GLint vertexAttrib = glGetAttribLocation(shaderProgramID(), attribute.name.c_str());
    glEnableVertexAttribArray(vertexAttrib);
    glVertexAttribPointer(
      vertexAttrib, attribute.size, GL_FLOAT, GL_FALSE, pVertAttrTotLen * pDataTypesSize,
      (void*) ((long) (offset * pDataTypesSize)));
    offset += attribute.size;
  }
}

// generalization NEEDED
void OpenGLHandler::pointsToVertexBufferData(Data<Point2D>& points) {
  const unsigned int size = points.size();
  float* vertexBufferData = new float[size * pVertAttrTotLen];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * pVertAttrTotLen + 0] = points[i].x;
    vertexBufferData[i * pVertAttrTotLen + 1] = points[i].y;
    vertexBufferData[i * pVertAttrTotLen + 2] = 0.1f;   // radius
    vertexBufferData[i * pVertAttrTotLen + 3] = 12.0f;  // steps
  }
  pVertexBufferData = Data(vertexBufferData, pVertAttrTotLen * size);
}
