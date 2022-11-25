#pragma once

#include <vector>

#include "graph/geometry.hpp"

class Node {
public:
  virtual unsigned int& index()      = 0;
  virtual unsigned int index() const = 0;
};

class IndexNode : Node {
  IndexNode() = default;
  IndexNode(unsigned int v) : pIndex(v) {}

  ~IndexNode() = default;

  void operator=(unsigned int v) { pIndex = v; }
  unsigned int& index() override { return pIndex; }
  unsigned int index() const override { return pIndex; }

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
  std::vector<Node> pNodes;
};

class CompleteGraph : Graph {
public:
  CompleteGraph()  = default;
  ~CompleteGraph() = default;

  bool adjacent([[maybe_unused]] const Node& u, [[maybe_unused]] const Node& v) override { return true; }
};

class Euclidean : CompleteGraph {
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
