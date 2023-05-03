/*
 * pathBTSP is a tool to solve, approximate and draw instances of BTSPP,
 * BTSP and TSP. Drawing is limited to euclidean graphs.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "solve/euclideandistancegraph.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "geometry.hpp"
#include "graph.hpp"

#include "solve/definitions.hpp"

graph::Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes) {
  std::array<uint_fast32_t, SEED_LENGTH> randomData;
  std::random_device src;
  std::generate(randomData.begin(), randomData.end(), std::ref(src));
  return generateEuclideanDistanceGraph(numOfNodes, randomData);
}

graph::Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes, const std::array<uint_fast32_t, SEED_LENGTH>& randomData) {
  std::seed_seq seed(randomData.begin(), randomData.end());
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
