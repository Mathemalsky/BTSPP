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

#include <memory>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/drawdata.hpp"
#include "draw/shader.hpp"

namespace drawing {
/*!
 * \brief draw draws the content on the screen
 * \param window pointer to the glfw window
 * \param programs containes all shaderprograms compiled in this project
 * \param drawData contains all data used for drawing including all buffers allocated in this project
 */
void draw(GLFWwindow* window, const ShaderProgramCollection& programs, const std::shared_ptr<DrawData> drawData);
}  // namespace drawing
