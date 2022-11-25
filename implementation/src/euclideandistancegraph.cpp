#include "euclideandistancegraph.hpp"

#include <random>
#include <vector>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes) {
  std::vector<Point2D> positions(numOfNodes);
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distribution(0, 1.0);
  for (Point2D& point : positions) {
    point.x = distribution(eng);
    point.y = distribution(eng);
  }
  return Euclidean(positions);
}
