#include "draw/draw.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "draw/forms.hpp"
#include "draw/variables.hpp"

#include "graph/graph.hpp"

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
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glColor3d(0.8, 0.8, 0.0);
  if (components.euclidean != nullptr) {
    drawEuclideanDistanceGraph(components.euclidean);
  }
}

void testdraw(GLuint& vertexbuffer) {
  // 1st attribute buffer : vertices
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  glVertexAttribPointer(
    0,         // attribute 0. No particular reason for 0, but must match the layout in the shader.
    2,         // size
    GL_FLOAT,  // type
    GL_FALSE,  // normalized?
    0,         // stride
    (void*) 0  // array buffer offset
  );
  // Draw the triangle !
  glDrawArrays(GL_TRIANGLES, 0, 3);  // Starting from vertex 0; 3 vertices total -> 1 triangle
  glDisableVertexAttribArray(0);
}
