#pragma once

#include <vector>

#include "graph/graph.hpp"

Edge findBottleneck(const Euclidean& euclidean, const std::vector<unsigned int>& tour, const bool cycle);
