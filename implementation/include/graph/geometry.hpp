#pragma once

#include <cmath>
#include <fstream>

struct LineSegment;

struct Point2D {
  double x, y;

  Point2D operator-(const Point2D& subtrahend) const { return Point2D{x - subtrahend.x, y - subtrahend.y}; }

  void operator-=(const Point2D& subtrahend) {
    x -= subtrahend.x;
    y -= subtrahend.y;
  }

  void operator+=(const Point2D& subtrahend) {
    x += subtrahend.x;
    y += subtrahend.y;
  }

  /*!
   * \brief liesOnLeftSide checks if p lies in the left half space defined by directed vector seg.a -> seg.b
   * \details check is performed by checking if the euclidean scalar product of (seg.a, seg.b) and the vector obtained
   * by rotating (seg.b, p) by \pi /2 is nonpositiv \param seg half space defining vector \return bool if p lies in the
   * left closed half space
   */
  bool liesOnLeftSide(const LineSegment& seg) const;
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

struct LineSegment {
  Point2D a;
  Point2D b;
};

inline Point2D direction(const LineSegment& seg) {
  return seg.b - seg.a;
}

inline Point2D direction(const Point2D& a, const Point2D& b) {
  return b - a;
}

inline bool Point2D::liesOnLeftSide(const LineSegment& seg) const {
  const Point2D dir1 = direction(seg);
  const Point2D dir2 = direction(seg.b, *this);
  return (dir1.y * dir2.x - dir1.x * dir2.y <= 0.0);
}

inline bool intersect(const LineSegment& seg1, const LineSegment& seg2) {
  return (seg2.a.liesOnLeftSide(seg1) != seg2.b.liesOnLeftSide(seg1))
         && (seg1.a.liesOnLeftSide(seg2) != seg1.b.liesOnLeftSide(seg2));
}
