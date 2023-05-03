#pragma once
#include <memory>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// graph library
#include "geometry.hpp"

#include "draw/buffers.hpp"
#include "draw/drawdata.hpp"
#include "draw/variables.hpp"

void keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);

void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);

void handleEvents(GLFWwindow* window, std::shared_ptr<drawing::DrawData> drawData);
