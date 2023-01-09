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

  size_t numberOfNodes() const { return pNumberOfNodes; }

  virtual double distance(const size_t u, const size_t v) const = 0;
  virtual double distance(const Edge e) const                   = 0;
  virtual bool adjacent(const size_t u, const size_t v) const   = 0;
  virtual bool connected() const                                = 0;

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

  Point2D position(const size_t v) const { return pPositions[v]; }
  std::vector<Point2D>& verteces() { return this->pPositions; }

  double distance(const size_t u, const size_t v) const override { return dist(pPositions[u], pPositions[v]); }
  double distance(const Edge e) const override { return dist(pPositions[e.u], pPositions[e.v]); }

private:
  std::vector<Point2D> pPositions;
};

class EdgeCost {
public:
  EdgeCost() = default;
  EdgeCost(const double cost) : pCost(cost) {}

  double cost() const { return pCost; }

  double operator()() const { return pCost; }
  void operator=(const double cost) { pCost = cost; }
  bool operator==(const double compare) { return pCost == compare; }

  EdgeCost operator+(const EdgeCost other) const { return EdgeCost(std::min(pCost, other.pCost)); }
  EdgeCost operator*(const EdgeCost other) const { return EdgeCost(pCost + other.pCost); }
  EdgeCost& operator+=(const EdgeCost other) {
    this->pCost += other.pCost;
    return *this;
  }

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

  ~SimpleGraph() = default;

  void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }

  virtual void addEdge(const size_t out, const size_t in, const EdgeCost edgeCost) = 0;
  virtual void addEdge(const Edge& e, const EdgeCost edgeCost)                     = 0;

protected:
  Eigen::SparseMatrix<EdgeCost> pAdjacencyMatrix;
  void warnLoop(const size_t node) { std::cerr << "[Graph] Warning: inserting a loop on node " << node << std::endl; }
};

class UndirectedGraph : public SimpleGraph {
public:
  UndirectedGraph() = default;
  UndirectedGraph(const std::vector<Eigen::Triplet<EdgeCost>>& tripletList) : SimpleGraph(tripletList) {}

  ~UndirectedGraph() = default;

  bool adjacent(const size_t u, const size_t v) const override {
    return (u > v ? pAdjacencyMatrix.coeff(u, v) == 0.0 : pAdjacencyMatrix.coeff(v, u) == 0.0);
  }

  bool connected() const override;
  bool connectedWhithout(const size_t vertex) const;
  bool biconnected() const;

  double distance(const size_t u, const size_t v) const override {
    return (u > v ? pAdjacencyMatrix.coeff(u, v).cost() : pAdjacencyMatrix.coeff(v, u).cost());
  }
  double distance(const Edge e) const override {
    return (e.u > e.v ? pAdjacencyMatrix.coeff(e.u, e.v).cost() : pAdjacencyMatrix.coeff(e.v, e.u).cost());
  }

  void addEdge(const size_t out, const size_t in, const EdgeCost edgeCost) override {
    out > in ? pAdjacencyMatrix.insert(out, in) = edgeCost : pAdjacencyMatrix.insert(in, out) = edgeCost;
  }

  void addEdge(const Edge& e, const EdgeCost edgeCost) override {
    e.u > e.v ? pAdjacencyMatrix.insert(e.u, e.v) = edgeCost : pAdjacencyMatrix.insert(e.u, e.v) = edgeCost;
  }
};

class Digraph : public SimpleGraph {
  Digraph()  = default;
  ~Digraph() = default;

  void addEdge(size_t out, size_t in, const EdgeCost edge) override {
    if (out == in) {
      warnLoop(out);
    }
    pAdjacencyMatrix.insert(out, in) = edge;
  }
};
