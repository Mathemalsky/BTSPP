#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/variables.hpp"
#include "draw/buffers.hpp"

#include "graph/geometry.hpp"

void keyCallback(
    [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);

inline Point2D transformCoordinates(const double x, const double y) {
  return Point2D{2.0 * x / mainwindow::WIDTH - 1.0, -2.0 * y / mainwindow::HEIGHT + 1.0};
}

void handleEvents(GLFWwindow* window, const Buffers& coordinates);
