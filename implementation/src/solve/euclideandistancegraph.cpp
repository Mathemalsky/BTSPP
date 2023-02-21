#include "solve/euclideandistancegraph.hpp"

#include <array>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes) {
  std::array<uint_fast32_t, 2> random_data{2333342850, 259139102};
  // std::random_device src;
  // std::generate(random_data.begin(), random_data.end(), std::ref(src));

  std::seed_seq seed(random_data.begin(), random_data.end());
  std::default_random_engine generator;
  generator.seed(seed);

  std::cerr << "seed: ";
  seed.param(std::ostream_iterator<uint_fast32_t>(std::cerr, " "));
  std::cerr << "\n";

  std::vector<Point2D> positions(numOfNodes);
  std::uniform_real_distribution<double> distribution(-0.95, 0.95);
  for (Point2D& point : positions) {
    point.x = distribution(generator);
    point.y = distribution(generator);
  }
  return Euclidean(std::move(positions));
}
