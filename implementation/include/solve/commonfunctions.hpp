#pragma once

#include <vector>

#include "graph/graph.hpp"

graph::Edge findBottleneck(const graph::Euclidean& euclidean, const std::vector<unsigned int>& tour, const bool cycle);
