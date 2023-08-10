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

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/drawdata.hpp"

/*!
 * \brief imguiVersionHints sets version of opengl
 * \details specify opengl version and glsl version
 * \return glsl version
 */
const char* imguiVersionHints();

/*!
 * \brief setUpImgui performs checks and sets up the imgui elements
 * \param window asociated instance of glfw window
 * \param glsl_version
 */
void setUpImgui(GLFWwindow* window, const char* glsl_version);

/*!
 * \brief drawImgui draws the updated gui for every frame
 */
void drawImgui(drawing::Appearance& appearance);

/*!
 * \brief cleanUpImgui frees memory for gui and performs proper shutdown of windows
 */
void cleanUpImgui();
