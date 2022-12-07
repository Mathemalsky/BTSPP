#include "draw/openglerrors.hpp"

#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

void getOpenglErrors(const char* file, const unsigned int line) {
  while (GLenum error = glGetError()) {
    std::cerr << "[OPEN_GL] ";
    switch (error) {
    case GL_INVALID_VALUE:
      std::cerr << "INVALID_VALUE\n";
      break;
    case GL_INVALID_OPERATION:
      std::cerr << "INVALID_OPERATION\n";
      break;
    default:
      std::cerr << error << std::endl;
    }
    std::cerr << "from " << file << ": line " << line << std::endl;
  }
}
