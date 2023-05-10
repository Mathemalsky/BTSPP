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
#include "draw/openglerrors.hpp"

#include <iostream>
#include <stdexcept>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

void getOpenglErrors(const char* file, const unsigned int line) {
  while (GLenum error = glGetError()) {
    std::cerr << "[OPEN_GL] ";
    switch (error) {
    case GL_INVALID_ENUM :
      std::cerr << "GL_INVALID_ENUM\n";
      break;
    case GL_INVALID_VALUE :
      std::cerr << "INVALID_VALUE\n";
      break;
    case GL_INVALID_OPERATION :
      std::cerr << "INVALID_OPERATION\n";
      break;
    default :
      std::cerr << error << std::endl;
    }
    std::cerr << "from " << file << ": line " << line << std::endl;
    throw std::runtime_error("[OPEN GL] OpenGL error: " + std::to_string(error));
  }
}
