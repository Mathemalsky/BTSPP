#pragma once

#include <algorithm>
#include <numeric>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"
//#include "graph/iterators.hpp"

struct Edge {
  size_t u;
  size_t v;
};

inline std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  return os << "(" << edge.u << ", " << edge.v << ") ";
}

/***********************************************************************************************************************
 *                                          iterator declaration
 **********************************************************************************************************************/

class GraphIt {
  virtual Edge operator*() const = 0;
  virtual GraphIt& operator++()  = 0;
};

class AdjMatGraphIt;
class DfsTreeIt;

/***********************************************************************************************************************
 *                                             graph classes
 **********************************************************************************************************************/

class EdgeWeight {
public:
  EdgeWeight()  = default;
  ~EdgeWeight() = default;

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
  virtual bool adjacent(const size_t u, const size_t v) const = 0;
  virtual bool adjacent(const Edge& e) const                  = 0;
  virtual bool connected() const                              = 0;
  virtual size_t numberOfEdges() const                        = 0;
  virtual size_t numberOfNodes() const                        = 0;

protected:
};

template <class T>
class Iterator {
public:
  Iterator(const T& graph, size_t position) : pGraph(graph), pPosition(position) {}
  Edge operator*() const;
  Iterator<T>& operator++();
  bool operator!=(const Iterator<T>& other) const;

private:
  const T& pGraph;
  size_t pPosition;
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
  size_t numberOfEdges() const override { return numberOfNodes() * (numberOfNodes() - 1) / 2; }
};

class Euclidean : public CompleteGraph, public WeightedGraph {
public:
  Euclidean()          = default;
  virtual ~Euclidean() = default;

  Euclidean(const std::vector<Point2D>& positions) : pPositions(positions) {}
  Euclidean(std::vector<Point2D>&& positions) : pPositions(positions) {}

  size_t numberOfNodes() const override { return pPositions.size(); }

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

class AdjMatGraph : public Modifyable, public WeightedGraph {
public:
  AdjMatGraph()  = default;
  ~AdjMatGraph() = default;

  AdjMatGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    pAdjacencyMatrix(Eigen::SparseMatrix<EdgeWeight>(numberOfNodes, numberOfNodes)) {
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  virtual bool adjacent(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v) == 0.0; }
  virtual bool adjacent(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v) == 0.0; }

  virtual double weight(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v).cost(); }

  size_t numberOfEdges() const override { return pAdjacencyMatrix.nonZeros(); }
  size_t numberOfNodes() const override { return pAdjacencyMatrix.cols(); }

  void compressMatrix() { pAdjacencyMatrix.makeCompressed(); }
  const Eigen::SparseMatrix<EdgeWeight>& matrix() const { return this->pAdjacencyMatrix; }

  void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }  // operator*= not provided by Eigen

  AdjMatGraphIt begin() const;
  AdjMatGraphIt end() const;

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
class AdjacencyMatrixGraph : public virtual AdjMatGraph {
public:
  AdjacencyMatrixGraph()  = default;
  ~AdjacencyMatrixGraph() = default;

  AdjacencyMatrixGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {}

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
  AdjacencyListGraph()  = default;
  ~AdjacencyListGraph() = default;

  AdjacencyListGraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }

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
  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  const std::vector<size_t>& neighbours(const size_t u) const { return pAdjacencyList[u]; }

private:
  std::vector<std::vector<size_t>> pAdjacencyList;
};

class Tree : public virtual Graph {
  bool connected() const override { return true; }
  size_t numberOfEdges() const override { return numberOfNodes() - 1; }
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

  DfsTree(const size_t numberOfNodes) {
    pAdjacencyList.resize(numberOfNodes);
    pIndeces.resize(numberOfNodes);
  }

  bool adjacent(const size_t u, const size_t v) const override { return v == parent(u); }
  bool adjacent(const Edge& e) const override { return e.v == parent(e.u); }

  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  // size_t index(const size_t u) const { return this->pIndeces[u]; }
  // size_t& index(const size_t u) { return this->pIndeces[u]; }

  size_t parent(const size_t u) const { return this->pAdjacencyList[u]; }
  size_t& parent(const size_t u) { return this->pAdjacencyList[u]; }

  DfsTreeIt begin() const;
  DfsTreeIt end() const;

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
 *                                          iterator implementation
 **********************************************************************************************************************/

class AdjMatGraphIt : public GraphIt {
public:
  struct SparseMatrixPos {
    size_t innerIndex;
    size_t outerIndex;
    bool operator!=(const SparseMatrixPos& other) const { return innerIndex != other.innerIndex; }
  };

  AdjMatGraphIt(const AdjMatGraph& graph, SparseMatrixPos position) : pGraph(graph), pPosition(position) {}

  Edge operator*() const override {
    assert(pGraph.matrix().isCompressed() && "check is only valid on compressed matrix");
    const int* innerIndeces = pGraph.matrix().innerIndexPtr();
    return Edge{(size_t) (innerIndeces[pPosition.innerIndex]), pPosition.outerIndex};
  }

  AdjMatGraphIt& operator++() override {
    assert(pGraph.matrix().isCompressed() && "check is only valid on compressed matrix");
    const int* outerIndeces = pGraph.matrix().outerIndexPtr();
    ++pPosition.innerIndex;
    if (pPosition.innerIndex >= (size_t) outerIndeces[pPosition.outerIndex + 1]) {
      ++pPosition.outerIndex;
    }
    return *this;
  }

  bool operator!=(const AdjMatGraphIt& other) const { return pPosition != other.pPosition; }

private:
  const AdjMatGraph& pGraph;
  SparseMatrixPos pPosition;
};

inline AdjMatGraphIt AdjMatGraph::begin() const {
  return AdjMatGraphIt(*this, AdjMatGraphIt::SparseMatrixPos{0, 0});
}
inline AdjMatGraphIt AdjMatGraph::end() const {
  return AdjMatGraphIt(
      *this, AdjMatGraphIt::SparseMatrixPos{(size_t) pAdjacencyMatrix.nonZeros(), (size_t) pAdjacencyMatrix.cols()});
}

class DfsTreeIt : public GraphIt {
public:
  DfsTreeIt(const DfsTree& tree, const size_t position) : pTree(tree), pPosition(position) {}

  Edge operator*() const override { return Edge{pPosition, pTree.parent(pPosition)}; }

  DfsTreeIt& operator++() override {
    ++pPosition;
    return *this;
  }

  bool operator!=(const DfsTreeIt& other) const { return pPosition != other.pPosition; }

private:
  const DfsTree& pTree;
  size_t pPosition;
};

inline DfsTreeIt DfsTree::begin() const {
  return DfsTreeIt(*this, (size_t) 1);
}
inline DfsTreeIt DfsTree::end() const {
  return DfsTreeIt(*this, pAdjacencyList.size());
}

/***********************************************************************************************************************
 *                                          template implementation
 **********************************************************************************************************************/

// DEBUG
#include <iostream>
inline std::ostream& operator<<(std::ostream& os, const AdjMatGraph& graph) {
  os << "Number of nodes: " << graph.numberOfNodes() << std::endl;
  for (const Edge& e : graph) {
    os << e << std::endl;
  }
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const DfsTree& tree) {
  os << "Number of nodes: " << tree.numberOfNodes() << std::endl;
  for (const Edge& e : tree) {
    os << e << std::endl;
  }
  return os;
}

template <Directionality DIRECT>
DfsTree dfs(const AdjacencyMatrixGraph<DIRECT>& graph, const size_t rootNode) {
  const size_t numberOfNodes = graph.numberOfNodes();
  DfsTree tree(numberOfNodes);
  std::vector<bool> visited(numberOfNodes, false);
  std::stack<size_t> nodeStack;
  nodeStack.push(rootNode);
  while (!nodeStack.empty()) {
    const size_t top = nodeStack.top();
    nodeStack.pop();
    if (!visited[top]) {
      visited[top] = true;
      for (unsigned int i = 0; i < top; ++i) {
        if (graph.matrix().coeff(top, i) != 0 && !visited[i]) {
          nodeStack.push(i);
          tree.parent(i) = top;
        }
      }
      for (Eigen::SparseMatrix<EdgeWeight>::InnerIterator it(graph.matrix(), top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());     // do not push already visited nodes
          tree.parent(it.index()) = top;  // update the parent
        }
      }
    }
  }
  return tree;
}
