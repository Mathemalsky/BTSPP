#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

/*!
 * \brief setUpImgui performs checks and sets up the imgui elements
 * \param window asociated instance of glfw window
 * \param glsl_version
 */
void setUpImgui(GLFWwindow* window, const char* glsl_version);

/*!
 * \brief initImGuiWindows asignes initial values to the elements
 */
void initImGuiWindows();

/*!
 * \brief drawImgui draws the updated gui for every frame
 */
void drawImgui();

/*!
 * \brief cleanUpImgui frees memory for gui and performs proper shutdown of windows
 */
void cleanUpImgui();
