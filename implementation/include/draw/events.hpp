#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

void keyCallback(
  [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);

void handleFastEvents(GLFWwindow* window);
