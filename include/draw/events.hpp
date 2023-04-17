#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/buffers.hpp"
#include "draw/drawdata.hpp"
#include "draw/variables.hpp"

#include "graph/geometry.hpp"

void keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);

void handleEvents(GLFWwindow* window, drawing::DrawData& drawData);
