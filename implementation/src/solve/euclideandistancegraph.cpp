#include "solve/euclideandistancegraph.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "graph/geometry.hpp"
#include "graph/graph.hpp"

graph::Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes) {
  std::array<uint_fast32_t, 2> random_data{2407893799, 27884848};
  std::random_device src;
  std::generate(random_data.begin(), random_data.end(), std::ref(src));

  std::seed_seq seed(random_data.begin(), random_data.end());
  std::default_random_engine generator;
  generator.seed(seed);

  std::cerr << "seed: ";
  seed.param(std::ostream_iterator<uint_fast32_t>(std::cerr, " "));
  std::cerr << "\n";

  std::vector<graph::Point2D> positions(numOfNodes);
  std::uniform_real_distribution<double> distribution(-0.95, 0.95);
  for (graph::Point2D& point : positions) {
    point.x = distribution(generator);
    point.y = distribution(generator);
  }
  return graph::Euclidean(std::move(positions));
}
