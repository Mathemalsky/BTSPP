#pragma once

#include <cmath>
#include <fstream>

struct Point2D {
  double x, y;

  Point2D operator-(const Point2D& subtrahend) { return Point2D{x - subtrahend.x, y - subtrahend.y}; }

  void operator-=(const Point2D& subtrahend) {
    x -= subtrahend.x;
    y -= subtrahend.y;
  }

  void operator+=(const Point2D& subtrahend) {
    x += subtrahend.x;
    y += subtrahend.y;
  }
};

inline double norm2(const double x, const double y) {
  return std::sqrt(x * x + y * y);
}
inline double norm2(const Point2D& p) {
  return norm2(p.x, p.y);
}
inline double dist(const Point2D& p, const Point2D& q) {
  return norm2(p.x - q.x, p.y - q.y);
}

inline std::ostream& operator<<(std::ostream& os, const Point2D& point) {
  return os << "(" << point.x << ", " << point.y << ") ";
}
