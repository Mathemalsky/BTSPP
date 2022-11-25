#ifndef GRAPH_GEOMETRY_HPP
#define GRAPH_GEOMETRY_HPP

#include <cmath>

struct Point2D {
  double x, y;
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

#endif
