/*
 * pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
 * BTSP and TSP. Drawing is limited to euclidean graphs.
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
#include "draw/buffers.hpp"

#include <memory>

template <>
VertexBuffer::VertexBuffer(const std::vector<float>& dat, const GLuint componentsPerVertex) :
  pComponentsPerVertex(componentsPerVertex),
  pType(GL_FLOAT) {
  GL_CALL(glGenBuffers(1, &pID);)
  this->bind();
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, bytes_of(dat), dat.data(), GL_DYNAMIC_DRAW);)
}

void VertexArray::mapBufferToAttribute(const std::shared_ptr<VertexBuffer> vbo, const GLuint shaderProgramID, const char* name) const {
  this->bind();
  vbo->bind();
  GL_CALL(const GLint vertexAttribLocation = glGetAttribLocation(shaderProgramID, name);)
  GL_CALL(glVertexAttribPointer(vertexAttribLocation, vbo->compPerVertex(), vbo->type(), GL_FALSE, 0, nullptr);)
}

void VertexArray::enable(const GLuint shaderProgramID, const char* name) const {
  this->bind();
  GL_CALL(const GLint vertexAttribLocation = glGetAttribLocation(shaderProgramID, name);)
  GL_CALL(glEnableVertexAttribArray(vertexAttribLocation);)
}

void VertexArray::bindBufferBase(const std::shared_ptr<ShaderBuffer> shaderBuffer, const GLuint bindingPoint) const {
  shaderBuffer->bind();
  this->bind();
  GL_CALL(glBindBufferBase(GL_SHADER_STORAGE_BUFFER, bindingPoint, shaderBuffer->id());)
}
