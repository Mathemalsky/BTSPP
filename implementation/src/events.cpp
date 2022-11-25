#include "events.hpp"

#include <GLFW/glfw3.h>

#include "definitions.hpp"
#include "variables.hpp"

void toggleSettingsWindow() {
  if (imguiwindow::SHOW_SETTINGS_WINDOW == true) {
    imguiwindow::SHOW_SETTINGS_WINDOW = false;
  }
  else {
    imguiwindow::SHOW_SETTINGS_WINDOW = true;
  }
}

void keyCallback(
  [[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods) {
  if (key == GLFW_KEY_F3 && action == GLFW_PRESS) {
    toggleSettingsWindow();
  }
}
