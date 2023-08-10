/*
 * BTSPP is a tool to solve, approximate and draw instances of BTSVPP,
 * BTSPP, BTSP and TSP. Drawing is limited to euclidean graphs.
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
#include "draw/drawdata.hpp"

#include <string.h>
#include <array>
#include <utility>
#include <vector>

#include "draw/definitions.hpp"

#include "solve/definitions.hpp"

namespace drawing {
void VertexOrder::updateOrder(const std::vector<size_t>& order, const ProblemType& type) {
  const std::vector<uint32_t> order32bit(order.begin(), order.end());  // cast to 32 bit integer
  if (type == ProblemType::BTSP_approx || type == ProblemType::BTSP_exact || type == ProblemType::TSP_exact) {
    pVertexOrder[std::to_underlying(type)].resize(order32bit.size() + PATH_OVERHEAD);
    std::memcpy(pVertexOrder[std::to_underlying(type)].data(), order32bit.data(), bytes_of(order32bit));
    std::memcpy(pVertexOrder[std::to_underlying(type)].data() + order32bit.size(), order32bit.data(), PATH_OVERHEAD * sizeof(uint32_t));
  }
  else if (type == ProblemType::BTSPP_approx || type == ProblemType::BTSPP_exact) {
    pVertexOrder[std::to_underlying(type)].resize(order32bit.size() + PATH_OVERHEAD - 1);  // n-1 path segments to draw
    pVertexOrder[std::to_underlying(type)][0] = order32bit[1];
    std::memcpy(pVertexOrder[std::to_underlying(type)].data() + 1, order32bit.data(), bytes_of(order32bit));
    pVertexOrder[std::to_underlying(type)].back() = order32bit[order32bit.size() - 2];
  }
  pInitialized[std::to_underlying(type)] = true;
}
}  // namespace drawing
