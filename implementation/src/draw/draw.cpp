#include "draw/draw.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/forms.hpp"
#include "draw/variables.hpp"

#include "graph/graph.hpp"

// DEBUG
#include <iostream>

void drawEuclideanDistanceGraph(const Euclidean* graph) {
  const std::vector<Point2D>& positions = graph->allPositions();
  for (const Point2D& position : positions) {
    drawCircle(position, 3);
  }
}

void draw(GLFWwindow* window, const DrawComponents& components) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glColor3d(0.8, 0.8, 0.0);
  if (components.euclidean != nullptr) {
    // std::cerr << "drawing\n";
    drawEuclideanDistanceGraph(components.euclidean);
  }

  // DEBUG
  // drawCircle(Point2D{0.0, 0.0}, 0.3);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glColor3d(0.0, 1.0, 0.0);
  glBegin(GL_POLYGON);
  glVertex2d(100, 200);
  glVertex2d(700, 200);
  glVertex2d(700, 700);
  glVertex2d(100, 700);
  glEnd();
}
