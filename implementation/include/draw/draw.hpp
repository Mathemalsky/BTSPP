#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/shader.hpp"

/*!
 * \brief draw draws the content on the screen
 * \param window pointer to the glfw window
 */
void draw(GLFWwindow* window, const ShaderProgramCollection& programs, Buffers& buffers);
