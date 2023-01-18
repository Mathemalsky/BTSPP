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

  virtual bool adjacent(const size_t u, const size_t v) const = 0;
  virtual bool connected() const                              = 0;
  virtual double weight(const size_t u, const size_t v) const { return (double) adjacent(u, v); }
  virtual double weight(const Edge& e) const { return (double) adjacent(e.u, e.v); }

protected:
  size_t pNumberOfNodes;
};

class WeightedGraph : public virtual Graph {
public:
  WeightedGraph()  = default;
  ~WeightedGraph() = default;

  WeightedGraph(const size_t numberOfNodes) : Graph(numberOfNodes) {}

  virtual double weight(const size_t u, const size_t v) const = 0;
  virtual double weight(const Edge& e) const                  = 0;
};

class CompleteGraph : public virtual Graph {
public:
  CompleteGraph()  = default;
  ~CompleteGraph() = default;

  CompleteGraph(const size_t numberOfNodes) : Graph(numberOfNodes) {}

  bool adjacent([[maybe_unused]] const size_t u, [[maybe_unused]] const size_t v) const override { return true; }
  bool connected() const override { return true; }
};

class Euclidean : public CompleteGraph, WeightedGraph {
public:
  Euclidean()  = default;
  ~Euclidean() = default;

  Euclidean(const std::vector<Point2D>& positions) : CompleteGraph(positions.size()), pPositions(positions) {}
  Euclidean(std::vector<Point2D>&& positions) : CompleteGraph(positions.size()), pPositions(positions) {}

  double weight(const size_t u, const size_t v) const override { return dist(pPositions[u], pPositions[v]); }
  double weight(const Edge& e) const override { return dist(pPositions[e.u], pPositions[e.v]); }

  Point2D position(const size_t v) const { return pPositions[v]; }
  std::vector<Point2D>& verteces() { return this->pPositions; }

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
    this->pCost = std::min(this->pCost, other.pCost);
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

class SimpleGraph : public virtual Graph {
  virtual void addEdge(const size_t out, const size_t in, const EdgeCost edgeCost) = 0;
  virtual void addEdge(const Edge& e, const EdgeCost edgeCost)                     = 0;
};

class AdjacencyMatrixGraph : public SimpleGraph {
public:
  AdjacencyMatrixGraph()  = default;
  ~AdjacencyMatrixGraph() = default;

  AdjacencyMatrixGraph(const size_t numberOfNodes)
    : Graph(numberOfNodes), pAdjacencyMatrix(Eigen::SparseMatrix<EdgeCost>(numberOfNodes, numberOfNodes)) {}
  AdjacencyMatrixGraph(const std::vector<Eigen::Triplet<EdgeCost>>& tripletList) : Graph(tripletList.size()) {
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void addEdge(const size_t out, const size_t in, const EdgeCost edgeCost) override {
    pAdjacencyMatrix.insert(out, in) = edgeCost;
  }

  virtual void addEdge(const Edge& e, const EdgeCost edgeCost) override {
    pAdjacencyMatrix.insert(e.u, e.v) = edgeCost;
  }

  virtual bool adjacent(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v) == 0.0; }

  bool connected() const override;

  virtual double weight(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v).cost(); }
  virtual double weight(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v).cost(); }

  bool connectedWhithout(const size_t vertex) const;
  bool biconnected() const;

  const Eigen::SparseMatrix<EdgeCost>& matrix() const { return this->pAdjacencyMatrix; }

  void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }  // operator*= not provided by Eigen

protected:
  Eigen::SparseMatrix<EdgeCost> pAdjacencyMatrix;
  void warnLoop(const size_t node) { std::cerr << "[Graph] Warning: inserting a loop on node " << node << std::endl; }
};

template <class Storage>
class UndirectedGraph : public SimpleGraph {
public:
  UndirectedGraph()  = default;
  ~UndirectedGraph() = default;

  // UndirectedGraph(const std::vector<Eigen::Triplet<EdgeCost>>& tripletList) : AdjacencyMatrixGraph(tripletList) {}

  void addEdge(const size_t out, const size_t in, const EdgeCost edgeCost) override {
    out > in ? Storage::addEdge(out, in, edgeCost) : Storage::adEdge(in, out, edgeCost);
  }

  void addEdge(const Edge& e, const EdgeCost edgeCost) override {
    e.u > e.v ? Storage::addEdge(e, edgeCost) : Storage::adEdge(e.v, e.u, edgeCost);
  }

  bool adjacent(const size_t u, const size_t v) const override {
    return (u > v ? Storage::adjacent(u, v) : Storage::adjacent(v, u));
  }

  double weight(const size_t u, const size_t v) const override {
    return (u > v ? Storage::weight(u, v) : Storage::weight(v, u));
  }

  double weight(const Edge& e) const override { return (e.u > e.v ? Storage::weight(e) : Storage::weight(e.v, e.u)); }
};

template <class Storage>
class Digraph : public SimpleGraph {
  Digraph()  = default;
  ~Digraph() = default;

  void addEdge(const size_t out, const size_t in, const EdgeCost edgeCost) override {
    Storage::addEdge(out, in, edgeCost);
  }
};

class AdjacencyListGraph : public Graph {
  AdjacencyListGraph()  = default;
  ~AdjacencyListGraph() = default;

protected:
  struct arc {
    size_t neighboor;
    double dist;
  };
  std::vector<std::vector<arc>> pAdjacencyList;
};

class Tree : public AdjacencyListGraph {
public:
  bool connected() const override { return true; }

protected:
  size_t pRoot;
};

/*!
 * \brief class FixTree
 * \details This class is for trees where every node has exactly one parent. This allows to store the neighboors
 * more efficient. The node 0 is assumed to be the root node.
 */
class DfsTree : public Graph {
public:
  DfsTree()  = default;
  ~DfsTree() = default;

  DfsTree(const size_t numberOfNodes) : Graph(numberOfNodes) {
    pAdjacencyList.resize(numberOfNodes);
    pIndeces.resize(numberOfNodes);
  }

  bool adjacent(const size_t u, const size_t v) const override { return u == parent(v); }
  bool connected() const override { return true; }

  size_t index(const size_t u) const { return this->pIndeces[u]; }
  size_t& index(const size_t u) { return this->pIndeces[u]; }

  size_t parent(const size_t u) const { return this->pAdjacencyList[u]; }
  size_t& parent(const size_t u) { return this->pAdjacencyList[u]; }

private:
  std::vector<size_t> pAdjacencyList;
  std::vector<size_t> pIndeces;
};

struct OpenEarDecomposition {
  bool biconnected;
};

OpenEarDecomposition schmidt(const UndirectedGraph<AdjacencyMatrixGraph>& graph);
