#pragma once

#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"

class Node {
public:
  Node() = default;
  Node(unsigned int v) : pIndex(v) {}

  ~Node() = default;

  void operator=(unsigned int v) { pIndex = v; }
  bool operator<(const Node other) const { return pIndex < other.pIndex; }
  bool operator>(const Node other) const { return pIndex > other.pIndex; }
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

  unsigned int numberOfNodes() { return pNumberOfNodes; }

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
  Euclidean(const std::vector<Point2D>& positions) : pPositions(positions) { pNumberOfNodes = positions.size(); }

  ~Euclidean() = default;

  double distance(const Node& u, const Node& v) override { return dist(pPositions[u.index()], pPositions[v.index()]); }
  Point2D position(const Node& v) const { return pPositions[v.index()]; }
  const std::vector<Point2D>& allPositions() const { return pPositions; }
  Point2D* pointer() { return &pPositions[0]; }

private:
  std::vector<Point2D> pPositions;
};

class EdgeCost {
public:
  EdgeCost() = default;
  EdgeCost(const double cost) : pCost(cost) {}

  double operator()() const { return pCost; }
  void operator=(const double cost) { pCost = cost; }

  EdgeCost operator+(const EdgeCost other) { return std::min(pCost, other.pCost); }
  EdgeCost operator*(const EdgeCost other) { return pCost + other.pCost; }

private:
  double pCost;
};

/*
 * Eigen needs some hints to deal with the custom type.
 */
namespace Eigen {
template <>
struct NumTraits<EdgeCost> : GenericNumTraits<EdgeCost> {
  typedef EdgeCost Real;
  typedef EdgeCost NonInteger;
  typedef EdgeCost Nested;

  enum {
    IsInteger             = 0,
    IsSigned              = 1,
    IsComplex             = 0,
    RequireInitialization = 0,
    ReadCost              = 1,
    AddCost               = 7,
    MulCost               = 3
  };
};
}  // namespace Eigen

class SimpleGraph : public Graph {
public:
  SimpleGraph() = default;
  SimpleGraph(const size_t numberOfNodes)
    : pAdjacencyMatrix(Eigen::SparseMatrix<EdgeCost>(numberOfNodes, numberOfNodes)) {}

  ~SimpleGraph();
  // constructor from 3 vectors

  // removal of a column and a row

  virtual void addEdge(const Node out, const Node in, const EdgeCost edge) = 0;

protected:
  Eigen::SparseMatrix<EdgeCost> pAdjacencyMatrix;
};

class UndirectedGraph : public SimpleGraph {
  void addEdge(const Node out, const Node in, const EdgeCost edge) override {
    if (out > in) {
      pAdjacencyMatrix.insert(out.index(), in.index()) = edge;
    }
    else if (out < in) {
    }
    else {
    }
  }
};

class Digraph : public SimpleGraph {};
