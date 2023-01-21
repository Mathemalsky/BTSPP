#pragma once

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglerrors.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

#include "graph/geometry.hpp"

#include "utility/utils.hpp"

/*! \file buffers.hpp */

/***********************************************************************************************************************
 *                                               VertexBuffer class
 **********************************************************************************************************************/

/*!
 * \brief VertexBuffer manages an OpenGL vertex buffer object
 */
class VertexBuffer {
public:
  /*!
   * \brief constructor, invokes glGenBuffers
   */
  VertexBuffer() { GL_CALL(glGenBuffers(1, &pID);) }

  /*!
   * \brief desstructor, invokes glDeleteBuffers
   */
  ~VertexBuffer() { GL_CALL(glDeleteBuffers(1, &pID);) }

  /*!
   * \brief binds this buffer as GL_ARRAY_BUFFER
   */
  void bind() const { GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pID);) }

  /*!
   * \brief copies dat to OpenGL
   * \details calls bind(), copies dat whith hint GL_DYNAMIC_DRAW to a new memory block associated with this buffer
   * \param dat data to be copied
   */
  template <typename Type>
  void bufferData(std::vector<Type>& dat, const GLuint componentsPerVertex);

  template <typename Type>
  void bufferSubData(std::vector<Type>& dat) const;

  GLenum type() const { return pType; }
  GLuint compPerVertex() const { return pComponentsPerVertex; }

private:
  GLuint pID; /**< this buffer's OpenGL ID */
  GLenum pType;
  GLuint pComponentsPerVertex;
};

/***********************************************************************************************************************
 *                                               ShaderBuffer class
 **********************************************************************************************************************/

/*!
 * \brief ShaderBuffer manages an OpenGL shader buffer object
 * \details ShaderBuffer is used to control the way of passing data to shader programs
 */
class ShaderBuffer {
public:
  /*!
   * \brief constructor, invokes glGenBuffers
   */
  ShaderBuffer() { GL_CALL(glGenBuffers(1, &pID);) }

  /*!
   * \brief destructor, invokes glDeleteBuffers
   */
  ~ShaderBuffer() { GL_CALL(glDeleteBuffers(1, &pID);) }

  /*!
   * \brief returns OpenGL's internal id for the shader buffer object
   * \return OpenGL's internal id for the shader buffer object
   */
  GLuint id() const { return pID; }

  /*!
   * \brief binds this buffer as GL_SHADER_STORAGE_BUFFER
   */
  void bind() const { GL_CALL(glBindBuffer(GL_SHADER_STORAGE_BUFFER, pID);) }

  /*!
   * \brief copies dat to OpenGL
   * \details calls bind(), copies dat whith hint GL_DYNAMIC_DRAW to a new memory block associated with this buffer
   * \param dat data to be copied
   */
  template <typename Type>
  void bufferData(std::vector<Type>& dat) const;

  /*!
   * \brief replaces data in OpenGL with dat
   * \details calls bind(), copies dat into existing memory block associated with this buffer
   * \param dat data to be copied
   */
  template <typename Type>
  void bufferSubData(const std::vector<Type>& dat) const;

private:
  GLuint pID; /**< this buffer's OpenGL ID */
};

/***********************************************************************************************************************
 *                                               VertexArray class
 **********************************************************************************************************************/

class VertexArray {
public:
  VertexArray() { GL_CALL(glGenVertexArrays(1, &pID);) }
  ~VertexArray() { GL_CALL(glDeleteVertexArrays(1, &pID);) }

  void bind() const { GL_CALL(glBindVertexArray(pID);) }

  void mapBufferToAttribute(const VertexBuffer& vbo, const GLuint shaderProgramID, const char* name);

  void enable(const GLuint shaderProgramID, const char* name) const;

  void bindBufferBase(const ShaderBuffer& shaderBuffer, const GLuint bindingPoint) const;

private:
  GLuint pID;
};

/*!
 * \brief Buffers bundles various buffers
 */
struct Buffers {
  VertexBuffer coordinates;     /**< coordinates of graph vertices */
  ShaderBuffer tour;            /**< vertex indeces in order as they appear in the tour */
  ShaderBuffer tourCoordinates; /**< coordinates of graph vertices */
};

/***********************************************************************************************************************
 *                                          template implementation
 **********************************************************************************************************************/

template <typename Type>
void VertexBuffer::bufferData(std::vector<Type>& dat, const GLuint componentsPerVertex) {
  this->bind();
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, bytes_of(dat), dat.data(), GL_DYNAMIC_DRAW);)

  pComponentsPerVertex = componentsPerVertex;

  if (std::is_same<Type, float>{}) {
    pType = GL_FLOAT;
  }
  else {
    std::runtime_error("[Buffer Data] Type not yet supported.");
  }
}

template <typename Type>
void VertexBuffer::bufferSubData(std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, bytes_of(dat), dat.data());)  // 0 for no offset
}

template <typename Type>
void ShaderBuffer::bufferData(std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferData(GL_SHADER_STORAGE_BUFFER, bytes_of(dat), dat.data(), GL_DYNAMIC_DRAW);)
}

template <typename Type>
void ShaderBuffer::bufferSubData(const std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, bytes_of(dat), dat.data());)  // 0 for no offset
}
