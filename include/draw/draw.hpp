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
