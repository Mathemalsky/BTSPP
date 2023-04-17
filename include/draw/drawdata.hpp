#pragma once

#include <vector>

#include "draw/buffers.hpp"

#include "graph/graph.hpp"

class FloatVertices {
public:
  FloatVertices()  = default;
  ~FloatVertices() = default;

  /*!
   * @brief casts the doubles in coordinates of euclidean into floatVertices
   * @param euclidean
   */
  void updatePointsfFromEuclidean(graph::Euclidean& euclidean) {
    pFloatVertices.resize(2 * euclidean.numberOfNodes());
    const std::vector<graph::Point2D> points = euclidean.vertices();
    for (size_t i = 0; i < points.size(); ++i) {
      pFloatVertices[2 * i]     = static_cast<float>(points[i].x);
      pFloatVertices[2 * i + 1] = static_cast<float>(points[i].y);
    }
  }

  std::vector<float>& modify() { return pFloatVertices; }
  const std::vector<float>& read() const { return pFloatVertices; }

  const float& xCoord(const size_t u) const { return pFloatVertices[2 * u]; }
  const float& yCoord(const size_t u) const { return pFloatVertices[2 * u + 1]; }

private:
  std::vector<float> pFloatVertices;
};

struct Results {
  approximation::Result BTSP_APPROX_RESULT;
  approximation::Result BTSPP_APPROX_RESULT;
  exactsolver::Result BTSP_EXACT_RESULT;
  exactsolver::Result BTSPP_EXACT_RESULT;
};

struct DrawData {
  const Buffers& buffers;
  FloatVertices floatVertices;
};