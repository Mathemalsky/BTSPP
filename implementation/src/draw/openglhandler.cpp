#include "draw/openglhandler.hpp"

#include "draw/openglerrors.hpp"

OpenGLVariables::~OpenGLVariables() {
  pVertexBufferData.~Data();
}

void OpenGLVariables::pointsToVertexBufferData(const Data<Point2D>& points) {
  const size_t size        = points.size();
  const size_t pointOffset = pVertAttr["vertexPosition"].offset;
  float* vertexBufferData = new float[size * pVertAttr.attrLen()];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * pVertAttr.attrLen() + pointOffset + 0] = points[i].x;  // cannot use memcpy here because we
    vertexBufferData[i * pVertAttr.attrLen() + pointOffset + 1] = points[i].y;  // need to cast double to float
  }
  pVertexBufferData = Data(vertexBufferData, pVertAttr.attrLen() * size /*, true */);
}

void OpenGLVariables::updatePointsInVertexBufferData(const Data<Point2D>& points) {
  const size_t size   = points.size();
  const size_t offset = pVertAttr["vertexPosition"].offset;
  for (unsigned int i = 0; i < size; ++i) {
    pVertexBufferData[i * pVertAttr.attrLen() + offset + 0] = points[i].x;  // cannot use memcpy here because we
    pVertexBufferData[i * pVertAttr.attrLen() + offset + 1] = points[i].y;  // need to cast double to float
  }
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, 20 * pVertAttr.attrLen() * pTypeSize, pVertexBufferData.data());)
}
