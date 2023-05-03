#pragma once

#include <stdexcept>
#include <string>

namespace graph {
class Exception : public std::exception {
public:
  explicit Exception(const char* message) : pMessage(message) {}
  explicit Exception(const std::string& message) : pMessage(message) {}

  virtual ~Exception() noexcept {}
  virtual const char* what() const noexcept { return pMessage.c_str(); }

protected:
  std::string pMessage;
};

class InfesableRequest : public Exception {
public:
  InfesableRequest(const std::string& msg) : Exception(msg) {}
};
}  // namespace graph
