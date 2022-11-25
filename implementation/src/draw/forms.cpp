#include "draw/forms.hpp"

#include <cmath>

#include <GLFW/glfw3.h>

#include "draw/variables.hpp"

#include "graph/geometry.hpp"

void drawCircle(const Point2D& position, const double outerRadius) {
  glBegin(GL_TRIANGLE_FAN);

  const unsigned int steps = 20;

  // set center
  glVertex2d(position.x, position.y);
  for (unsigned int i = 0; i <= steps; i++) {
    // calculate angle
    double angle = double(i) / double(steps) * 2.0 * M_PI;

    // set point
    glVertex2d(
      mainwindow::WIDTH * position.x + outerRadius * sin(angle),
      mainwindow::HEIGHT * position.y + outerRadius * cos(angle));
  }
  glEnd();
}
