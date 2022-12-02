#pragma once

#include <cassert>
#include <cstddef>
#include <cstring>  // needed for memcpy
#include <string>
#include <unordered_map>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/shader.hpp"

#include "graph/geometry.hpp"

#include "utility/datacontainer.hpp"

struct VertexAttribute {
  VertexAttribute(const char* Name, const size_t Size, const size_t Offset) : name(Name), size(Size), offset(Offset) {}

  const char* name; /*!< name of attribute */
  size_t size;      /*!< number of variables */
  size_t offset;    /*!< offset in number of variables */
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

private:
  std::vector<VertexAttribute> pAttr;
  std::unordered_map<std::string, size_t> pPosition;
  size_t pVertAttrTotLen;
  size_t pTypeSize;
};

template <typename Type>
class OpenGLHandler {
public:
  OpenGLHandler() = default;

  ~OpenGLHandler();

  void addVertexArray();
  void addVertexBuffer();

  void emplaceBackAttribute(const char* name, unsigned int size) { pVertAttr.emplaceBack(name, size); };

  void enableAllVertexAttribArrays() { pVertAttr.enableAllToShaderProgram(pShaderProgramID); }

  void linkShaderProgram() { pShaderProgramID = linkShaders(); }

  GLuint shaderProgramID() const { return pShaderProgramID; }
  GLuint vertexArrayID() const { return pVertexArrayID; }
  GLuint vertexBufferID() const { return pVertexBufferID; }

  void pointsToVertexBufferData(const Data<Point2D>& points);
  void updatePointsInVertexBufferData(const Data<Point2D>& points);
  void vertexBufferDataToGL() {
    glBufferData(GL_ARRAY_BUFFER, 20 * pVertAttr.attrLen() * pTypeSize, pVertexBufferData.data(), GL_DYNAMIC_DRAW);
  }

private:
  static constexpr size_t pTypeSize = sizeof(Type);

  GLuint pShaderProgramID;
  GLuint pVertexArrayID;
  GLuint pVertexBufferID;

  VertexAttributes<Type> pVertAttr;
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
    glEnableVertexAttribArray(vertexAttrib);
    glVertexAttribPointer(
      vertexAttrib, attribute.size, GL_FLOAT, GL_FALSE, pVertAttrTotLen * pTypeSize,
      (void*) (attribute.offset * pTypeSize));
  }
}

template <typename Type>
OpenGLHandler<Type>::~OpenGLHandler() {
  // delete Opengl data
  glDeleteVertexArrays(1, &pVertexArrayID);
  glDeleteBuffers(1, &pVertexBufferID);
  glDeleteProgram(pShaderProgramID);

  pVertexBufferData.~Data();
}

template <typename Type>
void OpenGLHandler<Type>::addVertexArray() {
  glGenVertexArrays(1, &pVertexArrayID);
  glBindVertexArray(pVertexArrayID);
}

template <typename Type>
void OpenGLHandler<Type>::addVertexBuffer() {
  glGenBuffers(1, &pVertexBufferID);
  glBindBuffer(GL_ARRAY_BUFFER, pVertexBufferID);
}

template <typename Type>
void OpenGLHandler<Type>::pointsToVertexBufferData(const Data<Point2D>& points) {
  const size_t size         = points.size();
  const size_t pointOffset  = pVertAttr["vertexPosition"].offset;
  const size_t radiusOffset = pVertAttr["size"].offset;
  const size_t stepsOffset  = pVertAttr["steps"].offset;
  float* vertexBufferData   = new float[size * pVertAttr.attrLen()];
  for (unsigned int i = 0; i < size; ++i) {
    vertexBufferData[i * pVertAttr.attrLen() + pointOffset + 0] = points[i].x;  // cannot use memcpy here because we
    vertexBufferData[i * pVertAttr.attrLen() + pointOffset + 1] = points[i].y;  // need to cast double to float
    vertexBufferData[i * pVertAttr.attrLen() + radiusOffset]    = 0.1f;         // radius
    vertexBufferData[i * pVertAttr.attrLen() + stepsOffset]     = 12.0f;        // steps
  }
  pVertexBufferData = Data(vertexBufferData, pVertAttr.attrLen() * size /*, true */);
}

template <typename Type>
void OpenGLHandler<Type>::updatePointsInVertexBufferData(const Data<Point2D>& points) {
  const size_t size   = points.size();
  const size_t offset = pVertAttr["vertexPosition"].offset;
  for (unsigned int i = 0; i < size; ++i) {
    pVertexBufferData[i * pVertAttr.attrLen() + offset + 0] = points[i].x;  // cannot use memcpy here because we
    pVertexBufferData[i * pVertAttr.attrLen() + offset + 1] = points[i].y;  // need to cast double to float
  }
  glBufferSubData(GL_ARRAY_BUFFER, 0, 20 * pVertAttr.attrLen() * pTypeSize, pVertexBufferData.data());
}
