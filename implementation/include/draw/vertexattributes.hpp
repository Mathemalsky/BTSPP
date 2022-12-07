#pragma once

#include <cassert>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

struct VertexAttribute {
  VertexAttribute(const char* name, const size_t size, const size_t offset) : name(name), size(size), offset(offset) {}

  const char* name;    /*!< name of attribute */
  const size_t size;   /*!< number of variables */
  const size_t offset; /*!< offset in number of variables */
};

template <typename Type>
class VertexAttributes {
public:
  VertexAttributes() : pVertAttrTotLen(0), pTypeSize(sizeof(Type)) {}

  VertexAttribute& operator[](unsigned int index) {
    assert(index < pAttr.size() && "[Data] trying to access out of range");
    return pAttr[index];
  }
  const VertexAttribute& operator[](unsigned int index) const {
    assert(index < pAttr.size() && "[Data] trying to access out of range");
    return pAttr[index];
  }

  const VertexAttribute& operator[](const std::string& str) const {
    assert(pPosition.at(str) < pAttr.size() && "[Data] trying to access out of range");
    return pAttr[pPosition.at(str)];
  }

  size_t attrLen() const { return pVertAttrTotLen; }

  void emplaceBack(const char* name, const size_t size);

  void enableAllToShaderProgram(const GLuint shaderProgramID);

  void pointsToVertexBufferData(const Data<Point2D>& points);

  void vertexBufferDataToGL() {
    GL_CALL(glBufferData(
              GL_ARRAY_BUFFER, graph::POINTS.size() * pVertAttrTotLen * pTypeSize, pVertexBufferData.data(),
              GL_DYNAMIC_DRAW);)
  }

  void updatePointsInVertexBufferData(const Data<Point2D>& points);

private:
  std::vector<VertexAttribute> pAttr;
  std::unordered_map<std::string, size_t> pPosition;
  size_t pVertAttrTotLen;
  const size_t pTypeSize;

  Data<Type> pVertexBufferData;
};

/***********************************************************************************************************************
 *                          implementation of templates needs to be included in header
 **********************************************************************************************************************/

template <typename Type>
void VertexAttributes<Type>::emplaceBack(const char* name, const size_t size) {
  pPosition[name] = pAttr.size();
  pAttr.emplace_back(name, size, pVertAttrTotLen);
  pVertAttrTotLen += size;
}

template <typename Type>
void VertexAttributes<Type>::enableAllToShaderProgram(const GLuint shaderProgramID) {
  for (const VertexAttribute& attribute : pAttr) {
    const GLint vertexAttrib = glGetAttribLocation(shaderProgramID, attribute.name);
    GL_CALL(glEnableVertexAttribArray(vertexAttrib);)
    GL_CALL(glVertexAttribPointer(
              vertexAttrib, attribute.size, GL_FLOAT, GL_FALSE, pVertAttrTotLen * pTypeSize,
              (void*) (attribute.offset * pTypeSize));)
  }
}

template <typename Type>
void VertexAttributes<Type>::pointsToVertexBufferData(const Data<Point2D>& points) {
  const size_t size        = points.size();
  const size_t pointOffset = pAttr[pPosition["vertexPosition"]].offset;
  float* vertexBufferData  = new float[size * pVertAttrTotLen];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * pVertAttrTotLen + pointOffset + 0] = points[i].x;  // cannot use memcpy here because we
    vertexBufferData[i * pVertAttrTotLen + pointOffset + 1] = points[i].y;  // need to cast double to float
  }
  pVertexBufferData.set(vertexBufferData, pVertAttrTotLen * size, true);
}

template <typename Type>
void VertexAttributes<Type>::updatePointsInVertexBufferData(const Data<Point2D>& points) {
  const size_t size   = points.size();
  const size_t offset = pAttr[pPosition["vertexPosition"]].offset;
  for (unsigned int i = 0; i < size; ++i) {
    pVertexBufferData[i * pVertAttrTotLen + offset + 0] = points[i].x;  // cannot use memcpy here because we
    pVertexBufferData[i * pVertAttrTotLen + offset + 1] = points[i].y;  // need to cast double to float
  }
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, points.size() * pVertAttrTotLen * pTypeSize, pVertexBufferData.data());)
}

inline void vertexArray() {
  GLuint vertexArrayObject;
  GL_CALL(glGenVertexArrays(1, &vertexArrayObject);)
  GL_CALL(glBindVertexArray(vertexArrayObject);)
}

inline void vertexBuffer() {
  GLuint vertexBufferObject;
  GL_CALL(glGenBuffers(1, &vertexBufferObject);)
  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, vertexBufferObject);)
}
