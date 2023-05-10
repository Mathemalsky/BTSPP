/*
 * GRAPH is a library to store and manipulate graphs as adjacency list or
 * as sparse eigen matrix. Different specialized types of graphs are
 * supported.
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

#include <cmath>
#include <fstream>

namespace graph {

struct LineSegment;
struct Vector2D;

using Point2D = Vector2D;

/*!
 * @brief Vector2D represents a vector or point in 2D
 */
struct Vector2D {
  //!@{
  double x, y; /**< coordinates */
  //!@}

  /*!
   * @brief subtract Vector2Ds componentwise
   * @param subtrahend will be subtracted
   */
  Vector2D operator-(const Vector2D& subtrahend) const { return Vector2D{x - subtrahend.x, y - subtrahend.y}; }

  /*!
   * @brief subtract Vector2D componentwise
   * @param subtrahend will be subtracted
   */
  void operator-=(const Vector2D& subtrahend) {
    x -= subtrahend.x;
    y -= subtrahend.y;
  }

  /*!
   * @brief subtract Vector2D componentwise
   * @param summand will be added
   */
  void operator+=(const Vector2D& summand) {
    x += summand.x;
    y += summand.y;
  }

  /*!
   * \brief liesOnLeftSide checks if p lies in the left half space defined by directed vector seg.a -> seg.b
   * \details check is performed by checking if the euclidean scalar product of (seg.a, seg.b) and the vector obtained
   * by rotating (seg.b, p) by \pi /2 is nonpositiv
   * \param seg half space defining vector
   * \return bool if p lies in the left closed half space
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
 * \param v vector
 * \return euclidean norm (2 norm)
 */
inline double norm2(const Vector2D& v) {
  return norm2(v.x, v.y);
}

/*!
 * \brief norm2squared computes the square of the 2 norm
 * \param x
 * \param y
 * \return square of euclidean norm (2 norm) of the vector (x, y)
 */
inline double norm2squared(const double x, const double y) {
  return x * x + y * y;
}

/*!
 * \brief norm2squared computes the square of the 2 norm
 * \param v vector
 * \return square of euclidean norm (2 norm)
 */
inline double norm2squared(const Vector2D& v) {
  return norm2squared(v.x, v.y);
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

/*!
 * \brief distSquared computes the square of the euclidean distance between p and q
 * \param p point
 * \param q point
 * \return square of euclidean norm (2 norm) of (p - q)
 */
inline double distSquared(const Point2D& p, const Point2D& q) {
  return norm2squared(p.x - q.x, p.y - q.y);
}

/*!
 * @brief add point in paranthesis to output stream
 * @param os output stream
 * @param point point or vector to be printed
 */
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

/*!
 * @brief checks if the line segment lies on the left side of the vector
 * @param seg is the line segment to check
 * @return true if seg lies on left side, false otherwise
 * @note inline keyword is needed for linker
 */
inline bool Vector2D::liesOnLeftSide(const LineSegment& seg) const {
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
  return (seg2.a.liesOnLeftSide(seg1) != seg2.b.liesOnLeftSide(seg1)) && (seg1.a.liesOnLeftSide(seg2) != seg1.b.liesOnLeftSide(seg2));
}
}  // namespace graph
