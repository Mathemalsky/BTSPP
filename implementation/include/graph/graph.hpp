#pragma once

#include <iostream>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"

struct Edge {
  size_t u;
  size_t v;
};

class Graph {
public:
  Graph()  = default;
  ~Graph() = default;

  Graph(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

  virtual double distance(const size_t u, const size_t v) const = 0;
  virtual double distance(const Edge e) const                   = 0;
  virtual bool adjacent(const size_t u, const size_t v) const   = 0;
  virtual bool connected() const                                = 0;

  size_t numberOfNodes() const { return pNumberOfNodes; }

protected:
  size_t pNumberOfNodes;
};

class CompleteGraph : public Graph {
public:
  CompleteGraph()  = default;
  ~CompleteGraph() = default;

  CompleteGraph(const size_t numberOfNodes) : Graph(numberOfNodes) {}

  bool adjacent([[maybe_unused]] const size_t u, [[maybe_unused]] const size_t v) const override { return true; }
  bool connected() const override { return true; }
};

class Euclidean : public CompleteGraph {
public:
  Euclidean() = default;
  Euclidean(const std::vector<Point2D>& positions) : CompleteGraph(positions.size()), pPositions(positions) {}
  Euclidean(std::vector<Point2D>&& positions) : CompleteGraph(positions.size()), pPositions(positions) {}

  ~Euclidean() = default;

  double distance(const size_t u, const size_t v) const override { return dist(pPositions[u], pPositions[v]); }
  double distance(const Edge e) const override { return dist(pPositions[e.u], pPositions[e.v]); }
  Point2D position(const size_t v) const { return pPositions[v]; }
  std::vector<Point2D>& verteces() { return this->pPositions; }

private:
  std::vector<Point2D> pPositions;
};

class EdgeCost {
public:
  EdgeCost() = default;
  EdgeCost(const double cost) : pCost(cost) {}

  double operator()() const { return pCost; }
  void operator=(const double cost) { pCost = cost; }

  EdgeCost operator+(const EdgeCost other) const { return std::min(pCost, other.pCost); }
  EdgeCost operator*(const EdgeCost other) const { return pCost + other.pCost; }

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
    ReadCost              = 5,
    AddCost               = 1,
    MulCost               = 1
  };
};
}  // namespace Eigen

class SimpleGraph : public Graph {
public:
  SimpleGraph() = default;
  SimpleGraph(const size_t numberOfNodes)
    : Graph(numberOfNodes), pAdjacencyMatrix(Eigen::SparseMatrix<EdgeCost>(numberOfNodes, numberOfNodes)) {}
  SimpleGraph(const std::vector<Eigen::Triplet<EdgeCost>>& tripletList) : Graph(tripletList.size()) {
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  ~SimpleGraph();

  // removal of a column and a row

  virtual void addEdge(const size_t out, const size_t in, const EdgeCost edge) = 0;

protected:
  Eigen::SparseMatrix<EdgeCost> pAdjacencyMatrix;
  void warnLoop(const size_t node) { std::cerr << "[Graph] Warning: inserting a loop on node " << node << std::endl; }
};

class UndirectedGraph : public SimpleGraph {
  void addEdge(const size_t out, const size_t in, const EdgeCost edge) override {
    out > in ? pAdjacencyMatrix.insert(out, in) = edge : pAdjacencyMatrix.insert(in, out);
  }

  bool connected() const override;
};

class Digraph : public SimpleGraph {
  void addEdge(size_t out, size_t in, const EdgeCost edge) override {
    if (out == in) {
      warnLoop(out);
    }
    pAdjacencyMatrix.insert(out, in) = edge;
  }
};
