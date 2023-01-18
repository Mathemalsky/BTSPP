#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"

struct Edge {
  size_t u;
  size_t v;
};

enum class Directionality { Undirected, Directed };

template <Directionality direct>
static Edge orientation(const size_t u, const size_t v) {
  if constexpr (direct == Directionality::Undirected) {
    return (u > v ? Edge{u, v} : Edge{v, u});
  }
  else {
    return Edge{u, v};
  }
}

template <Directionality direct>
Edge orientation(const Edge& e) {
  if constexpr (direct == Directionality::Undirected) {
    return (e.u > e.v ? e : Edge{e.v, e.u});
  }
  else {
    return e;
  }
}

class Graph {
public:
  Graph()  = default;
  ~Graph() = default;

  Graph(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

  virtual bool adjacent(const size_t u, const size_t v) const = 0;
  virtual bool adjacent(const Edge& e) const                  = 0;
  virtual bool connected() const                              = 0;

  size_t numberOfNodes() const { return pNumberOfNodes; }

protected:
  size_t pNumberOfNodes;
};

class WeightedGraph : public virtual Graph {
public:
  /*
  WeightedGraph()  = default;
  ~WeightedGraph() = default;

  WeightedGraph(const size_t numberOfNodes) : Graph(numberOfNodes) {}
  */
  virtual double weight(const size_t u, const size_t v) const = 0;
  virtual double weight(const Edge& e) const                  = 0;
};

class CompleteGraph : public virtual Graph {
public:
  /*
  CompleteGraph()  = default;
  ~CompleteGraph() = default;

  CompleteGraph(const size_t numberOfNodes) : Graph(numberOfNodes) {}
  */

  bool adjacent([[maybe_unused]] const size_t u, [[maybe_unused]] const size_t v) const override { return true; }
  bool adjacent([[maybe_unused]] const Edge& e) const override { return true; }
  bool connected() const override { return true; }
};

class Euclidean : public CompleteGraph, WeightedGraph {
public:
  Euclidean()  = default;
  ~Euclidean() = default;

  Euclidean(const std::vector<Point2D>& positions) : Graph(positions.size()), pPositions(positions) {}
  Euclidean(std::vector<Point2D>&& positions) : Graph(positions.size()), pPositions(positions) {}

  double weight(const size_t u, const size_t v) const override { return dist(pPositions[u], pPositions[v]); }
  double weight(const Edge& e) const override { return dist(pPositions[e.u], pPositions[e.v]); }

  Point2D position(const size_t v) const { return pPositions[v]; }
  std::vector<Point2D>& verteces() { return this->pPositions; }

private:
  std::vector<Point2D> pPositions;
};

class EdgeWeight {
public:
  EdgeWeight() = default;
  EdgeWeight(const double cost) : pCost(cost) {}

  double cost() const { return pCost; }

  double operator()() const { return pCost; }
  void operator=(const double cost) { pCost = cost; }
  bool operator==(const double compare) { return pCost == compare; }

  EdgeWeight operator+(const EdgeWeight other) const { return EdgeWeight(std::min(pCost, other.pCost)); }
  EdgeWeight operator*(const EdgeWeight other) const { return EdgeWeight(pCost + other.pCost); }
  EdgeWeight& operator+=(const EdgeWeight other) {
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
struct NumTraits<EdgeWeight> : GenericNumTraits<EdgeWeight> {
  typedef EdgeWeight Real;
  typedef EdgeWeight NonInteger;
  typedef EdgeWeight Nested;

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

class Modifyable : public virtual Graph {
protected:
  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) = 0;
  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight)                     = 0;
};

class AdjMatGraph : public virtual WeightedGraph, Modifyable {
public:
  AdjMatGraph()  = default;
  ~AdjMatGraph() = default;

  AdjMatGraph(const size_t numberOfNodes)
    : Graph(numberOfNodes), pAdjacencyMatrix(Eigen::SparseMatrix<EdgeWeight>(numberOfNodes, numberOfNodes)) {}
  AdjMatGraph(const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) : Graph(tripletList.size()) {
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  virtual bool adjacent(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v) == 0.0; }
  virtual bool adjacent(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v) == 0.0; }

  bool connected() const override;

  // double weight(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v).cost(); }
  virtual double weight(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v).cost(); }

  bool biconnected() const;
  bool connectedWhithout(const size_t vertex) const;

  const Eigen::SparseMatrix<EdgeWeight>& matrix() const { return this->pAdjacencyMatrix; }

  void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }  // operator*= not provided by Eigen

protected:
  Eigen::SparseMatrix<EdgeWeight> pAdjacencyMatrix;

  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    pAdjacencyMatrix.insert(out, in) = edgeWeight;
  }

  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight) override {
    pAdjacencyMatrix.insert(e.u, e.v) = edgeWeight;
  }
};

template <Directionality direct>
class AdjacencyMatrixGraph : public AdjMatGraph {
public:
  AdjacencyMatrixGraph()  = default;
  ~AdjacencyMatrixGraph() = default;

  AdjacencyMatrixGraph(const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) : AdjMatGraph(tripletList) {}

  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(orientation<direct>(out, in), edgeWeight);
  }

  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(orientation<direct>(e), edgeWeight);
  }

  bool adjacent(const size_t u, const size_t v) const override {
    return AdjMatGraph::adjacent(orientation<direct>(u, v));
  }

  double weight(const size_t u, const size_t v) const override {
    return AdjMatGraph::weight(orientation<direct>(u, v));
  }

  double weight(const Edge& e) const override { return AdjMatGraph::weight(orientation<direct>(e)); }
};

class AdjListGraph : public Modifyable {
public:
  virtual bool adjacent(const size_t u, const size_t v) const;
};

template <Directionality direct>
class AdjacencyListGraph : public AdjListGraph {
public:
  bool adjacent(const size_t u, const size_t v) const override {
    const Edge e                         = orientation<direct>(u, v);
    const std::vector<size_t>& neighbour = pAdjacencyList[e.u];
    const auto it                        = std::find(neighbour.begin(), neighbour.end(), e.v);
    return it != neighbour.end();
  }

private:
  std::vector<std::vector<size_t>> pAdjacencyList;
};

/*
template<class Directed>
class WeightedAdjacencyListGraph : public AdjListGraph, WeightedGraph {
public:
  bool adjacent(const size_t u, const size_t v) const override {
    const Edge e = Directed::orientation(u,v);

  }

private:
  struct arc {
    size_t neighboor;
    double dist;
  };
  std::vector<std::vector<arc>> pAdjacencyList;
}
*/

template <Directionality Directed>
class Tree : public AdjacencyListGraph<Directed> {
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
  bool adjacent(const Edge& e) const override { return e.u == parent(e.v); }
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

OpenEarDecomposition schmidt(const AdjacencyMatrixGraph<Directionality::Undirected>& graph);
