#include <cstdio>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// error callback function which prints glfw errors in case they arise
static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  // set error colback function
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit()) {
    return -1;
  }
  return 0;
}
