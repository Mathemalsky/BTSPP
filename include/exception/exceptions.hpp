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
#pragma once

#include <stdexcept>
#include <string>

class Exception : public std::exception {
public:
  explicit Exception(const char* message) : pMessage(message) {}
  explicit Exception(const std::string& message) : pMessage(message) {}

  virtual ~Exception() noexcept {}
  virtual const char* what() const noexcept { return pMessage.c_str(); }

protected:
  std::string pMessage;
};

class InvalidArgument : public Exception {
public:
  InvalidArgument(const std::string& msg) : Exception(msg) {}
};

class UnknownType : public Exception {
public:
  UnknownType(const std::string& msg) : Exception(msg) {}
};
