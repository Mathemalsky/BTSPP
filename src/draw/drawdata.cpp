#include "draw/drawdata.hpp"

#include <string.h>
#include <array>
#include <utility>
#include <vector>

#include "draw/definitions.hpp"

#include "solve/definitions.hpp"

namespace drawing {
void VertexOrder::updateOrder(const std::vector<unsigned int>& order, const ProblemType& type) {
  if (type == ProblemType::BTSP_approx || type == ProblemType::BTSP_exact || type == ProblemType::TSP_exact) {
    pVertexOrder[std::to_underlying(type)].resize(order.size() + PATH_OVERHEAD);
    std::memcpy(pVertexOrder[std::to_underlying(type)].data(), order.data(), bytes_of(order));
    std::memcpy(pVertexOrder[std::to_underlying(type)].data() + order.size(), order.data(), PATH_OVERHEAD * sizeof(unsigned int));
  }
  else if (type == ProblemType::BTSPP_approx || type == ProblemType::BTSPP_exact) {
    pVertexOrder[std::to_underlying(type)].resize(order.size() + PATH_OVERHEAD - 1);  // n-1 path segments to draw
    pVertexOrder[std::to_underlying(type)][0] = order[1];
    std::memcpy(pVertexOrder[std::to_underlying(type)].data() + 1, order.data(), bytes_of(order));
    pVertexOrder[std::to_underlying(type)].back() = order[order.size() - 2];
  }
  pInitialized[std::to_underlying(type)] = true;
}
}  // namespace drawing
