#pragma once

#include <algorithm>
#include <numeric>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"

struct Edge {
  size_t u;
  size_t v;
};

inline std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  return os << "(" << edge.u << ", " << edge.v << ") ";
}

class EdgeWeight {
public:
  EdgeWeight() = default;
  EdgeWeight(const double cost) : pCost(cost) {}

  double cost() const { return pCost; }

  double operator()() const { return pCost; }
  void operator=(const double cost) { pCost = cost; }
  bool operator==(const double compare) { return pCost == compare; }
  bool operator!=(const double compare) { return pCost != compare; }

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

enum class Directionality { Undirected, Directed };

template <Directionality DIRECT>
static Edge orientation(const size_t u, const size_t v) {
  if constexpr (DIRECT == Directionality::Undirected) {
    return (u > v ? Edge{u, v} : Edge{v, u});
  }
  else {
    return Edge{u, v};
  }
}

template <Directionality DIRECT>
static Edge orientation(const Edge& e) {
  if constexpr (DIRECT == Directionality::Undirected) {
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
  virtual size_t numberOfEdges() const                        = 0;

  size_t numberOfNodes() const { return pNumberOfNodes; }

protected:
  size_t pNumberOfNodes;
};

class WeightedGraph : public virtual Graph {
public:
  virtual double weight(const size_t u, const size_t v) const = 0;
  virtual double weight(const Edge& e) const                  = 0;
};

class CompleteGraph : public virtual Graph {
public:
  bool adjacent([[maybe_unused]] const size_t u, [[maybe_unused]] const size_t v) const override { return true; }
  bool adjacent([[maybe_unused]] const Edge& e) const override { return true; }
  bool connected() const override { return true; }
  size_t numberOfEdges() const override { return pNumberOfNodes * (pNumberOfNodes - 1) / 2; }
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

class Modifyable : public virtual Graph {
protected:
  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) = 0;
  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight)                     = 0;
};

// DBEUG
#include <iostream>

class AdjMatGraph : public Modifyable, WeightedGraph {
public:
  AdjMatGraph()  = default;
  ~AdjMatGraph() = default;

  AdjMatGraph(const size_t numberOfNodes) :
    Graph(numberOfNodes), pAdjacencyMatrix(Eigen::SparseMatrix<EdgeWeight>(numberOfNodes, numberOfNodes)) {}
  AdjMatGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) {
    // DEBUG
    pNumberOfNodes = numberOfNodes;

    pAdjacencyMatrix = Eigen::SparseMatrix<EdgeWeight>(numberOfNodes, numberOfNodes);
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

    // DEBUG
    std::cerr << "#nodes in AdjMatGraphClass: " << pNumberOfNodes << std::endl;
    std::cerr << "#matrix size in AdjMatGraphClass: " << pAdjacencyMatrix.rows() << "x" << pAdjacencyMatrix.cols()
              << std::endl;
  }

  virtual bool adjacent(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v) == 0.0; }
  virtual bool adjacent(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v) == 0.0; }

  virtual double weight(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v).cost(); }

  size_t numberOfEdges() const override { return pAdjacencyMatrix.nonZeros(); }

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

template <Directionality DIRECT>
class AdjacencyMatrixGraph : public AdjMatGraph {
public:
  AdjacencyMatrixGraph()  = default;
  ~AdjacencyMatrixGraph() = default;

  AdjacencyMatrixGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {
    // DEBUG
    std::cerr << "#nodes in AdjacencyMatrixGraphClass: " << pNumberOfNodes << std::endl;
    std::cerr << "#matrix size in AdjacencyMatrixGraphClass: " << pAdjacencyMatrix.rows() << "x"
              << pAdjacencyMatrix.cols() << std::endl;
  }

  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(orientation<DIRECT>(out, in), edgeWeight);
  }

  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(orientation<DIRECT>(e), edgeWeight);
  }

  bool adjacent(const size_t u, const size_t v) const override {
    return AdjMatGraph::adjacent(orientation<DIRECT>(u, v));
  }

  bool connected() const override;

  double weight(const size_t u, const size_t v) const override {
    return AdjMatGraph::weight(orientation<DIRECT>(u, v));
  }

  double weight(const Edge& e) const override { return AdjMatGraph::weight(orientation<DIRECT>(e)); }

  bool biconnected() const;
  bool connectedWhithout(const size_t vertex) const;
};

template <Directionality DIRECT>
class AdjacencyListGraph : public Modifyable {
public:
  void addEdge(const size_t out, const size_t in, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    const Edge e = orientation<DIRECT>(out, in);
    pAdjacencyList[e.u].push_back(e.v);
  }

  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    const Edge newEdge = orientation<DIRECT>(e);
    pAdjacencyList[newEdge.u].push_back(newEdge.v);
  }

  bool adjacent(const size_t u, const size_t v) const override {
    const Edge e                         = orientation<DIRECT>(u, v);
    const std::vector<size_t>& neighbour = pAdjacencyList[e.u];
    const auto it                        = std::find(neighbour.begin(), neighbour.end(), e.v);
    return it != neighbour.end();
  }

  bool adjacent(const Edge& e) const { return this->adjacent(e.u, e.v); }

  bool connected() const override;

  size_t numberOfEdges() const override {
    return std::accumulate(
        pAdjacencyList.begin(), pAdjacencyList.end(), 0,
        [](const unsigned int sum, const std::vector<size_t>& vec) { return sum + vec.size(); });
  }

  const std::vector<size_t>& neighbours(const size_t u) const { return pAdjacencyList[u]; }

private:
  std::vector<std::vector<size_t>> pAdjacencyList;
};

class Tree : public virtual Graph {
  bool connected() const override { return true; }
  size_t numberOfEdges() const override { return pNumberOfNodes - 1; }
};

/*!
 * \brief class FixTree
 * \details This class is for trees where every node has exactly one parent. This allows to store the neighboors
 * more efficient. The node 0 is assumed to be the root node.
 */
class DfsTree : public Tree {
public:
  DfsTree()  = default;
  ~DfsTree() = default;

  DfsTree(const size_t numberOfNodes) : Graph(numberOfNodes) {
    pAdjacencyList.resize(numberOfNodes);
    pIndeces.resize(numberOfNodes);
  }

  bool adjacent(const size_t u, const size_t v) const override { return v == parent(u); }
  bool adjacent(const Edge& e) const override { return e.v == parent(e.u); }

  size_t index(const size_t u) const { return this->pIndeces[u]; }
  size_t& index(const size_t u) { return this->pIndeces[u]; }

  size_t parent(const size_t u) const { return this->pAdjacencyList[u]; }
  size_t& parent(const size_t u) { return this->pAdjacencyList[u]; }

private:
  std::vector<size_t> pAdjacencyList;
  std::vector<size_t> pIndeces;
};

template <Directionality DIRECT>
DfsTree dfs(const AdjacencyMatrixGraph<DIRECT>& graph, const size_t rootNode = 0);

struct OpenEarDecomposition {
  std::vector<std::vector<size_t>> ears;
};

OpenEarDecomposition schmidt(const AdjacencyMatrixGraph<Directionality::Undirected>& graph);

/***********************************************************************************************************************
 *                                          template implementation
 **********************************************************************************************************************/

template <Directionality DIRECT>
DfsTree dfs(const AdjacencyMatrixGraph<DIRECT>& graph, const size_t rootNode) {
  const size_t numberOfNodes = graph.numberOfNodes();
  DfsTree tree(numberOfNodes);
  size_t indexCounter = 0;
  std::vector<bool> visited(numberOfNodes, false);
  std::stack<size_t> nodeStack;
  nodeStack.push(rootNode);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      visited[top] = true;
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(graph.matrix(), top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());     // do not push already visited nodes
          tree.parent(it.index()) = top;  // update the parent
        }
      }
      tree.index(top) = indexCounter;
      ++indexCounter;
    }
  }
  return tree;
}
