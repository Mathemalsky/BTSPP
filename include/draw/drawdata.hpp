#pragma once

#include <vector>

#include "draw/buffers.hpp"

#include "graph/graph.hpp"

class DrawGraph {
public:
  DrawGraph(const graph::Euclidean& euclidean) : pEuclidean(euclidean) {}
  ~DrawGraph() = default;

  const graph::Euclidean& euclidean() const { return pEuclidean; }
  graph::Euclidean& euclidean() { return pEuclidean; }

  void updatePointsfFromEuclidean() {
    pPoints_f.resize(2 * pEuclidean.numberOfNodes());
    const std::vector<graph::Point2D> points = pEuclidean.vertices();
    for (size_t i = 0; i < points.size(); ++i) {
      pPoints_f[2 * i]     = (float) points[i].x;
      pPoints_f[2 * i + 1] = (float) points[i].y;
    }
  }

private:
  graph::Euclidean pEuclidean;
  std::vector<float> pPoints_f;
};

struct Results {
  approximation::Result BTSP_APPROX_RESULT;
  approximation::Result BTSPP_APPROX_RESULT;
  exactsolver::Result BTSP_EXACT_RESULT;
  exactsolver::Result BTSPP_EXACT_RESULT;
};

struct DrawData {
  const Buffers& buffers;
  DrawGraph drawgraph;

  // calback functions as memberfunctions
  void keyCallback([[maybe_unused]] GLFWwindow* window, int key, [[maybe_unused]] int scancode, int action, [[maybe_unused]] int mods);
  void mouseButtonCallback([[maybe_unused]] GLFWwindow* window, int button, int action, [[maybe_unused]] int mods);
};