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
#pragma once

#include <array>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/openglerrors.hpp"

class ShaderProgram {
public:
  ShaderProgram() : pProgramID(glCreateProgram()) {}
  ~ShaderProgram() { GL_CALL(glDeleteProgram(pProgramID);) }

  GLuint id() const { return pProgramID; }

  void attachShader(const GLuint shader) const { GL_CALL(glAttachShader(pProgramID, shader);) }

  void link() const;

  void setUniform(const char* name, const float value) const;
  void setUniform(const char* name, const int value) const;
  void setUniform(const char* name, const float val1, const float val2) const;
  void setUniform(const char* name, const float val1, const float val2, const float val3, const float val4) const;
  void setUniform(const char* name, const std::array<float, 4>& value) const;

  void use() const { GL_CALL(glUseProgram(pProgramID);) }

private:
  const GLuint pProgramID;
};

class ShaderCollection {
public:
  ShaderCollection();
  ~ShaderCollection();

  ShaderProgram linkCircleDrawProgram() const;
  ShaderProgram linkPathSegementDrawProgram() const;
  ShaderProgram linkLineDrawProgram() const;

private:
  const GLuint pVertexShader;
  const GLuint pCircleShader;
  const GLuint pPathVertexShader;
  const GLuint pFragmentShader;
  const GLuint pLineVertexShader;
};

struct ShaderProgramCollection {
  ShaderProgramCollection(const ShaderProgram& drawCircles, const ShaderProgram& drawPathSegments, const ShaderProgram& drawLine) :
    drawCircles(drawCircles),
    drawPathSegments(drawPathSegments),
    drawLine(drawLine) {}
  const ShaderProgram drawCircles;
  const ShaderProgram drawPathSegments;
  const ShaderProgram drawLine;
};
