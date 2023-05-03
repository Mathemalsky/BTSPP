#pragma once

#include <array>

// graph library
#include "graph.hpp"

#include "solve/definitions.hpp"

graph::Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes);
graph::Euclidean generateEuclideanDistanceGraph(unsigned int numOfNodes, const std::array<uint_fast32_t, SEED_LENGTH>& randomData);
