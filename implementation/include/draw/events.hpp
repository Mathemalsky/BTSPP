#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/variables.hpp"
#include "draw/buffers.hpp"

#include "graph/geometry.hpp"

void keyCallback(
    [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);

void handleEvents(GLFWwindow* window, const Buffers& coordinates);
