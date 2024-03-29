/*
 * BTSPP is a tool to solve, approximate and draw instances of BTSVPP,
 * BTSPP, BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// graph library
#include "geometry.hpp"

#include "draw/openglerrors.hpp"
#include "draw/shader.hpp"
#include "draw/variables.hpp"

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
   * \brief constructor, copies dat to OpenGL
   * \details calls bind(), copies dat whith hint GL_DYNAMIC_DRAW to a new memory block associated with this buffer
   * \param dat data to be copied
   * \param componentsPerVertex stores how variables of type Type belong to each vertex. This is needed in
   * VertexArray::mapBufferToAttribute().
   */
  template <typename Type>
  VertexBuffer(const std::vector<Type>& dat, const GLuint componentsPerVertex);

  /*!
   * \brief desstructor, invokes glDeleteBuffers
   */
  ~VertexBuffer() { GL_CALL(glDeleteBuffers(1, &pID);) }

  /*!
   * \brief binds this buffer as GL_ARRAY_BUFFER
   */
  void bind() const { GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pID);) }

  /*!
   * \brief replaces data in OpenGL with dat
   * \details calls bind(), copies dat into existing memory block associated with this buffer
   * \param dat data to be copied
   */
  template <typename Type>
  void bufferSubData(const std::vector<Type>& dat) const;

  /*!
   * \brief compPerVertex returns the number of basic variables per vertex
   * \return number of basic variables per vertex
   */
  GLuint compPerVertex() const { return pComponentsPerVertex; }

  /*!
   * \brief type returns the buffers generic type
   * \return type of data stored in the buffer (usually float or unsigned int)
   */
  GLenum type() const { return pType; }

private:
  GLuint pID;                        /**< this buffer's OpenGL ID */
  const GLuint pComponentsPerVertex; /**< number of variables per vertex*/
  const GLenum pType;                /**< data type of variables in the buffer */
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
   * \brief constructor, invokes glGenBuffers, calls bufferData()
   */
  template <typename Type>
  ShaderBuffer(const std::vector<Type>& dat) {
    GL_CALL(glGenBuffers(1, &pID);)
    bufferData<Type>(dat);
  }

  /*!
   * \brief destructor, invokes glDeleteBuffers
   */
  ~ShaderBuffer() { GL_CALL(glDeleteBuffers(1, &pID);) }

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
  void bufferData(const std::vector<Type>& dat) const;

  /*!
   * \brief replaces data in OpenGL with dat
   * \details calls bind(), copies dat into existing memory block associated with this buffer
   * \param dat data to be copied
   */
  template <typename Type>
  void bufferSubData(const std::vector<Type>& dat) const;

  /*!
   * \brief returns OpenGL's internal id for the shader buffer object
   * \return OpenGL's internal id for the shader buffer object
   */
  GLuint id() const { return pID; }

private:
  GLuint pID; /**< this buffer's OpenGL ID */
};

/***********************************************************************************************************************
 *                                               VertexArray class
 **********************************************************************************************************************/

/*!
 * \brief The VertexArray class mamanages an OpenGL vertex array object
 */
class VertexArray {
public:
  /*!
   * \brief constructor, invokes glGenVertexArrays
   */
  VertexArray() { GL_CALL(glGenVertexArrays(1, &pID);) }

  /*!
   * \brief destructor, invokes glDeleteVertexArrays
   */
  ~VertexArray() { GL_CALL(glDeleteVertexArrays(1, &pID);) }

  /*!
   * \brief binds this VertexArray
   */
  void bind() const { GL_CALL(glBindVertexArray(pID);) }

  /*!
   * \brief mapBufferToAttribute maps data from vbo to an attribute of shaderprogram
   * \details calls VertexBuffer::bind() on the vertex buffer object vbo, finds location of attribute name from
   * shaderprogram with id shaderProgramID and maps data from vbo there \param vbo vertex buffer object that contains
   * the data \param shaderProgramID id of shaderprogram with attribute, which needs the data \param name the attributes
   * name
   */
  void mapBufferToAttribute(const std::shared_ptr<VertexBuffer> vbo, const GLuint shaderProgramID, const char* name) const;

  /*!
   * \brief enables the given attribute in given shaderprogram
   * \param shaderProgramID id of shaderprogram
   * \param name the attributes name
   */
  void enable(const GLuint shaderProgramID, const char* name) const;

  /*!
   * \brief bindBufferBase calls glBindBufferBase
   * \details calls ShaderBuffer::bind() on shaderBuffer, and the calls glBindBufferBase with the shaderBuffers ID and
   * bindingpoint as parameters \param shaderBuffer the ShaderBuffer to bind \param bindingPoint position to bind the
   * buffer (like an address)
   */
  void bindBufferBase(const std::shared_ptr<ShaderBuffer> shaderBuffer, const GLuint bindingPoint) const;

private:
  GLuint pID; /**< this buffers OpenGL ID */
};

/*!
 * \brief Buffers bundles various buffers
 */
struct Buffers {
  std::shared_ptr<VertexBuffer> coordinates;     /**< coordinates of graph vertices */
  std::shared_ptr<ShaderBuffer> tour;            /**< vertex indeces in order as they appear in the tour */
  std::shared_ptr<ShaderBuffer> tourCoordinates; /**< coordinates of graph vertices */
};

/***********************************************************************************************************************
 *                                          template implementation
 **********************************************************************************************************************/

template <typename Type>
void VertexBuffer::bufferSubData(const std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, bytes_of(dat), dat.data());)  // 0 for no offset
}

template <typename Type>
void ShaderBuffer::bufferData(const std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferData(GL_SHADER_STORAGE_BUFFER, bytes_of(dat), dat.data(), GL_DYNAMIC_DRAW);)
}

template <typename Type>
void ShaderBuffer::bufferSubData(const std::vector<Type>& dat) const {
  this->bind();
  GL_CALL(glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, bytes_of(dat), dat.data());)  // 0 for no offset
}
