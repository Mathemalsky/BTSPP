#include "solve/commonfunctions.hpp"

#include <vector>

#include "graph/graph.hpp"

graph::Edge findBottleneck(const graph::Euclidean& euclidean, const std::vector<unsigned int>& tour, const bool cycle) {
  unsigned int bottleneckEdgeEnd = 0;
  double bottleneckWeight        = euclidean.weight(tour[0], tour[1]);
  for (unsigned int i = 1; i < euclidean.numberOfNodes() - 1; ++i) {
    if (euclidean.weight(tour[i], tour[i + 1]) > bottleneckWeight) {
      bottleneckEdgeEnd = i;
      bottleneckWeight  = euclidean.weight(tour[i], tour[i + 1]);
    }
  }
  if (cycle && euclidean.weight(tour.back(), tour[0]) > bottleneckWeight) {
    return graph::Edge{tour.back(), 0};
  }
  else {
    return graph::Edge{tour[bottleneckEdgeEnd], tour[bottleneckEdgeEnd + 1]};
  }
}
