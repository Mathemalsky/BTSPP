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
  pVertexAttributesTotalLength += size;
}

void OpenGLHandler::pushBackAttribute(std::string name, unsigned int size) {
  pVertexAttributes.push_back(VertexAttribute{name, size});
  pVertexAttributesTotalLength += size;
}

void OpenGLHandler::enableAllVertexAttribArrays() {
  unsigned int offset = 0;
  for (const VertexAttribute& attribute : pVertexAttributes) {
    const GLint vertexAttrib = glGetAttribLocation(shaderProgramID(), attribute.name.c_str());
    glEnableVertexAttribArray(vertexAttrib);
    glVertexAttribPointer(
      vertexAttrib, attribute.size, GL_FLOAT, GL_FALSE, pVertexAttributesTotalLength * pDataTypesSize,
      (void*) ((long) (offset * pDataTypesSize)));
    offset += attribute.size;
  }
}

// generalization NEEDED
float* OpenGLHandler::euclideanDistanceGraphToVertexBufferData(const Euclidean* graph) {
  const std::vector<Point2D>& positions = graph->allPositions();
  const unsigned int size               = positions.size();
  float* vertexBufferData               = new float[size * pVertexAttributesTotalLength];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * pVertexAttributesTotalLength + 0] = positions[i].x;
    vertexBufferData[i * pVertexAttributesTotalLength + 1] = positions[i].y;
    vertexBufferData[i * pVertexAttributesTotalLength + 2] = 0.1f;   // radius
    vertexBufferData[i * pVertexAttributesTotalLength + 3] = 12.0f;  // steps
  }
  return vertexBufferData;
}
