#pragma once

#include <array>
#include <vector>

#include "draw/buffers.hpp"

#include "graph/graph.hpp"

#include "solve/definitions.hpp"

namespace drawing {
class FloatVertices {
public:
  FloatVertices()  = default;
  ~FloatVertices() = default;

  FloatVertices(const std::vector<float>& vertices) : pFloatVertices(vertices) {}

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

  /*!
   * @brief read only reference to the underlying float vector
   */
  const std::vector<float>& read() const { return pFloatVertices; }

  /*!
   * @brief access to x coordinate of vertex u
   * @param u vertex
   */
  const float& xCoord(const size_t u) const { return pFloatVertices[2 * u]; }

  /*!
   * @brief access to y coordinate of vertex u
   * @param u vertex
   */
  const float& yCoord(const size_t u) const { return pFloatVertices[2 * u + 1]; }

private:
  /*! vector of the vertices coordinates, order: x_0, y_0, x_1, y_1, ... */
  std::vector<float> pFloatVertices;
};

struct Results {
  approximation::Result BTSP_APPROX_RESULT;
  approximation::Result BTSPP_APPROX_RESULT;
  exactsolver::Result BTSP_EXACT_RESULT;
  exactsolver::Result BTSPP_EXACT_RESULT;
};

struct Appearance {
  std::array<RGBA_COLOUR, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> colour;
  std::array<float, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> thickness;
  RGBA_COLOUR clearColour;
  RGBA_COLOUR vertexColour;

  void init() {
    colour       = INITIAL_COLOUR;
    thickness    = INITIAL_THICKNESS;
    clearColour  = INITIAL_CLEAR_COLOUR;
    vertexColour = INITIAL_VERTEX_COLOUR;
  }
};

class VertexOrder {
public:
  VertexOrder() : pInitialized(INITIAL_ACTIVENESS) {}
  ~VertexOrder() = default;

  void updateOrder(const std::vector<unsigned int>& order, const ProblemType& type);
  const std::vector<unsigned int>& operator[](const ProblemType type) const { return pVertexOrder[std::to_underlying(type)]; }
  bool initialized(const ProblemType& type) const { return pInitialized[std::to_underlying(type)]; }

private:
  std::array<std::vector<unsigned int>, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> pVertexOrder;
  std::array<bool, std::to_underlying(ProblemType::NUMBER_OF_OPTIONS)> pInitialized;
};

struct DrawData {
  DrawData(const Buffers& buffers, const FloatVertices& floatVertices) : buffers(buffers), floatVertices(floatVertices) {
    appearance.init();
  }
  const Buffers& buffers;
  FloatVertices floatVertices;
  Results results;
  VertexOrder vertexOrder;
  Appearance appearance;
};
}  // namespace drawing