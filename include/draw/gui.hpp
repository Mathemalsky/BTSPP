#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

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
void drawImgui();

/*!
 * \brief cleanUpImgui frees memory for gui and performs proper shutdown of windows
 */
void cleanUpImgui();
