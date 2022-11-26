#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "graph/graph.hpp"

struct DrawComponents {
  const Euclidean* euclidean;
};

/*!
 * \brief draw draws the content on the screen
 * \param window pointer to the glfw window
 */
void draw(GLFWwindow* window, const DrawComponents& components);

void testdraw(GLuint &vertexbuffer);
