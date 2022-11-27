#pragma once

#include <vector>

#include "graph/geometry.hpp"

class Node {
public:
  Node() = default;
  Node(unsigned int v) : pIndex(v) {}

  ~Node() = default;

  void operator=(unsigned int v) { pIndex = v; }
  unsigned int& index() { return pIndex; }
  unsigned int index() const { return pIndex; }

private:
  unsigned int pIndex;
};

class Graph {
public:
  Graph()  = default;
  ~Graph() = default;

  virtual double distance(const Node& u, const Node& v) = 0;
  virtual bool adjacent(const Node& u, const Node& v)   = 0;

protected:
  unsigned int pNumberOfNodes;
};

class CompleteGraph : public Graph {
public:
  CompleteGraph()  = default;
  ~CompleteGraph() = default;

  bool adjacent([[maybe_unused]] const Node& u, [[maybe_unused]] const Node& v) override { return true; }
};

class Euclidean : public CompleteGraph {
public:
  Euclidean() = default;
  Euclidean(const std::vector<Point2D>& positions) : pPositions(positions) {}

  ~Euclidean() = default;

  double distance(const Node& u, const Node& v) override { return dist(pPositions[u.index()], pPositions[v.index()]); }
  Point2D position(const Node& v) const { return pPositions[v.index()]; }
  const std::vector<Point2D>& allPositions() const { return pPositions; }

private:
  std::vector<Point2D> pPositions;
};
