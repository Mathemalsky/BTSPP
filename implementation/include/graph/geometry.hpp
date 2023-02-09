#pragma once

#include <cmath>
#include <fstream>

struct LineSegment;
struct Vector2D;

using Point2D = Vector2D;

struct Vector2D {
  double x, y;

  Vector2D operator-(const Vector2D& subtrahend) const { return Vector2D{x - subtrahend.x, y - subtrahend.y}; }

  void operator-=(const Vector2D& subtrahend) {
    x -= subtrahend.x;
    y -= subtrahend.y;
  }

  void operator+=(const Vector2D& subtrahend) {
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

/*!
 * \brief norm2 computes the 2 norm
 * \param x
 * \param y
 * \return euclidean norm (2 norm) of the vector (x, y)
 */
inline double norm2(const double x, const double y) {
  return std::sqrt(x * x + y * y);
}

/*!
 * \brief norm2 computes the 2 norm of a vector
 * \varam v vector
 * \return euclidean norm (2 norm)
 */
inline double norm2(const Vector2D& v) {
  return norm2(v.x, v.y);
}

/*!
 * \brief dist messures the euclidean distance between p and q
 * \param p point
 * \param q point
 * \return euclidean norm (2 norm) of (p - q)
 */
inline double dist(const Point2D& p, const Point2D& q) {
  return norm2(p.x - q.x, p.y - q.y);
}

inline std::ostream& operator<<(std::ostream& os, const Point2D& point) {
  return os << "(" << point.x << ", " << point.y << ") ";
}

/*!
 * \brief The LineSegment struct represents a line segment between 2 points in euclidean plane
 */
struct LineSegment {
  Point2D a;
  Point2D b;
};

/*!
 * \brief direction computes the direction of a line segment from a to b
 * \param seg
 * \return direction (not normalized)
 */
inline Vector2D direction(const LineSegment& seg) {
  return seg.b - seg.a;
}

/*!
 * \brief direction computes the direction from a to b
 * \param a point
 * \param b point
 * \return direction (not normalized)
 */
inline Vector2D direction(const Point2D& a, const Point2D& b) {
  return b - a;
}

inline bool Point2D::liesOnLeftSide(const LineSegment& seg) const {
  const Vector2D dir1 = direction(seg);
  const Vector2D dir2 = direction(seg.b, *this);
  return (dir1.y * dir2.x - dir1.x * dir2.y <= 0.0);
}

/*!
 * \brief intersect checks if 2 line segments intersect
 * \param seg1 line segment
 * \param seg2 line segment
 * \return true if they intersect, else false
 */
inline bool intersect(const LineSegment& seg1, const LineSegment& seg2) {
  return (seg2.a.liesOnLeftSide(seg1) != seg2.b.liesOnLeftSide(seg1))
         && (seg1.a.liesOnLeftSide(seg2) != seg1.b.liesOnLeftSide(seg2));
}
