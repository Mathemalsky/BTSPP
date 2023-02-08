#pragma once

#include <algorithm>
#include <numeric>
#include <stack>
#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"

/***********************************************************************************************************************
 *                                                Edges
 **********************************************************************************************************************/

struct Edge {
  size_t u;
  size_t v;
  Edge reverse() const { return Edge{v, u}; }
  void invert() { std::swap(u, v); }
};

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
    return (e.u > e.v ? e : e.reverse());
  }
  else {
    return e;
  }
}

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

  bool operator<=(const EdgeWeight other) { return pCost <= other.pCost; }
  EdgeWeight operator+(const EdgeWeight other) const { return EdgeWeight(std::min(pCost, other.pCost)); }
  EdgeWeight operator*(const EdgeWeight other) const { return EdgeWeight(pCost + other.pCost); }
  EdgeWeight& operator+=(const EdgeWeight other) {
    pCost = std::min(pCost, other.pCost);
    return *this;
  }

  friend EdgeWeight abs(const EdgeWeight weight) { return EdgeWeight(std::abs(weight.pCost)); }

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

/***********************************************************************************************************************
 *                                          iterator declaration
 **********************************************************************************************************************/

class GraphIt {
  virtual Edge operator*() const = 0;
  virtual GraphIt& operator++()  = 0;
};

template <Directionality DIRECT>
class AdjacencyListGraphIt;
class AdjMatGraphIt;
class DfsTreeIt;

/***********************************************************************************************************************
 *                                             graph classes
 **********************************************************************************************************************/

class Graph {
public:
  virtual bool adjacent(const size_t u, const size_t v) const = 0;
  virtual bool adjacent(const Edge& e) const                  = 0;
  virtual bool connected() const                              = 0;
  virtual size_t numberOfEdges() const                        = 0;
  virtual size_t numberOfNodes() const                        = 0;

protected:
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
  std::vector<Point2D>& verteces() { return pPositions; }

private:
  std::vector<Point2D> pPositions;
};

class Modifyable : public virtual Graph {
protected:
  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) = 0;
  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight)                     = 0;
};

/*!
 * \brief The AdjacencyListGraph class implements a modifyable graph with adjacency list as internal storage
 * \details If the edge is directed, it is directed from vertex associated with outer storage index to vertex associated
 * with inner storage index
 */
template <Directionality DIRECT>
class AdjacencyListGraph : public Modifyable {
public:
  AdjacencyListGraph()  = default;
  ~AdjacencyListGraph() = default;

  AdjacencyListGraph(const AdjacencyListGraph& graph) = default;
  AdjacencyListGraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }
  AdjacencyListGraph(const std::vector<std::vector<size_t>>& adjacencyList) : pAdjacencyList(adjacencyList) {}

  void addEdge(const size_t out, const size_t in, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    pAdjacencyList[out].push_back(in);
    if constexpr (DIRECT == Directionality::Undirected) {
      pAdjacencyList[in].push_back(out);
    }
  }

  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    pAdjacencyList[e.u].push_back(e.v);
    if constexpr (DIRECT == Directionality::Undirected) {
      pAdjacencyList[e.v].push_back(e.u);
    }
  }

  bool adjacent(const size_t u, const size_t v) const override {
    const Edge e                         = orientation<DIRECT>(u, v);
    const std::vector<size_t>& neighbour = pAdjacencyList[e.u];
    const auto it                        = std::find(neighbour.begin(), neighbour.end(), e.v);
    return it != neighbour.end();
  }

  bool adjacent(const Edge& e) const { return adjacent(e.u, e.v); }

  /*!
   * \brief AdjacencyListGraph::connected checks the graph is connected
   * \details Check if the connected componend containing vertex 0 is the whole graph, by performing a dfs
   * \return true if the graph is connected, else false
   */
  bool connected() const override;

  size_t numberOfEdges() const override {
    return std::accumulate(
        pAdjacencyList.begin(), pAdjacencyList.end(), 0,
        [](const unsigned int sum, const std::vector<size_t>& vec) { return sum + vec.size(); });
  }
  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  AdjacencyListGraphIt<DIRECT> begin() const;
  AdjacencyListGraphIt<DIRECT> end() const;

  const std::vector<std::vector<size_t>>& adjacencyList() const { return pAdjacencyList; }

  AdjacencyListGraph<Directionality::Undirected> undirected() const;

  const std::vector<size_t>& neighbours(const size_t u) const { return pAdjacencyList[u]; }
  size_t numberOfNeighbours(const size_t u) const { return pAdjacencyList[u].size(); }

private:
  std::vector<std::vector<size_t>> pAdjacencyList;
};

class AdjMatGraph : public Modifyable, public WeightedGraph {
public:
  AdjMatGraph()  = default;
  ~AdjMatGraph() = default;

  AdjMatGraph(const size_t numberOfNodes) : pAdjacencyMatrix(numberOfNodes, numberOfNodes) {}
  AdjMatGraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : pAdjacencyMatrix(mat) {}
  AdjMatGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    pAdjacencyMatrix(Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>(numberOfNodes, numberOfNodes)) {
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  virtual bool adjacent(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v) == 0.0; }
  virtual bool adjacent(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v) == 0.0; }

  double weight(const Edge& e) const override {
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != 0 && "edgeweight 0 cann also mean the edge does not exists");
    return pAdjacencyMatrix.coeff(e.u, e.v).cost();
  }
  double weight(const size_t u, const size_t v) const override {
    assert(pAdjacencyMatrix.coeff(u, v) != 0 && "edgeweight 0 cann also mean the edge does not exists");
    return pAdjacencyMatrix.coeff(u, v).cost();
  }

  size_t numberOfEdges() const override { return pAdjacencyMatrix.nonZeros(); }
  size_t numberOfNodes() const override { return pAdjacencyMatrix.cols(); }

  AdjMatGraphIt begin() const;
  AdjMatGraphIt end() const;

  void compressMatrix() { pAdjacencyMatrix.makeCompressed(); }

  const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& matrix() const { return pAdjacencyMatrix; }

  void prune() { pAdjacencyMatrix.prune(0.0, Eigen::NumTraits<EdgeWeight>::dummy_precision()); }

  void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }  // operator*= not provided by Eigen

protected:
  Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor> pAdjacencyMatrix;

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

  AdjacencyMatrixGraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : AdjMatGraph(mat) {}
  AdjacencyMatrixGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {}

  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(out, in, edgeWeight);
    if constexpr (DIRECT == Directionality::Undirected) {
      AdjMatGraph::addEdge(in, out, edgeWeight);
    }
  }

  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { this->addEdge(e.u, e.v, edgeWeight); }

  bool connected() const override;

  bool biconnected() const;
  bool connectedWhithout(const size_t vertex) const;
  void removeEdge(const Edge& e) {
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != 0 && "edge to be removed does not exist in graph");
    pAdjacencyMatrix.coeffRef(e.u, e.v) = 0.0;
    if constexpr (DIRECT == Directionality::Undirected) {
      assert(pAdjacencyMatrix.coeff(e.v, e.u) != 0 && "edge to be removed does not exist in graph");
      pAdjacencyMatrix.coeffRef(e.v, e.u) = 0.0;
    }
  }
  AdjacencyMatrixGraph<Directionality::Undirected> removeUncriticalEdges() const;
  AdjacencyMatrixGraph<Directionality::Undirected> undirected() const;
};

class Tree : public virtual Graph {
  bool connected() const override { return true; }
  size_t numberOfEdges() const override { return numberOfNodes() - 1; }
};

/*!
 * \brief class FixTree
 * \details This class is for trees where every node has exactly one parent. This allows to store the neighboors
 * more efficient.
 */
class DfsTree : public Tree {
public:
  DfsTree()  = default;
  ~DfsTree() = default;

  DfsTree(const size_t numberOfNodes) {
    pAdjacencyList.resize(numberOfNodes);
    pExplorationOrder.resize(numberOfNodes);
  }

  bool adjacent(const size_t u, const size_t v) const override { return u != 0 && v == parent(u); }
  bool adjacent(const Edge& e) const override { return e.u != 0 && e.v == parent(e.u); }

  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  DfsTreeIt begin() const;
  DfsTreeIt end() const;

  std::vector<size_t> explorationOrder() const { return pExplorationOrder; }
  std::vector<size_t>& explorationOrder() { return pExplorationOrder; }

  size_t parent(const size_t u) const {
    assert(
        u != pExplorationOrder[0] && "You are trying to access the root nodes parent, which is uninitialized memory!");
    return pAdjacencyList[u];
  }

  size_t& parent(const size_t u) {
    assert(
        u != pExplorationOrder[0] && "You are trying to access the root nodes parent, which is uninitialized memory!");
    return pAdjacencyList[u];
  }

private:
  std::vector<size_t> pAdjacencyList;
  std::vector<size_t> pExplorationOrder;
};

/***********************************************************************************************************************
 *                                          iterator implementation
 **********************************************************************************************************************/

struct AdjListPos {
  size_t outerIndex;
  size_t innerIndex;
  bool operator!=(const AdjListPos& other) const { return outerIndex != other.outerIndex; }
};

template <Directionality DIRECT>
class AdjacencyListGraphIt : public GraphIt {
public:
  AdjacencyListGraphIt(const AdjacencyListGraph<DIRECT>& graph, AdjListPos position) :
    pGraph(graph), pPosition(position) {}

  Edge operator*() const override {
    return Edge{pPosition.outerIndex, pGraph.neighbours(pPosition.outerIndex)[pPosition.innerIndex]};
  }

  AdjacencyListGraphIt& operator++() override {
    ++pPosition.innerIndex;
    while (pPosition.innerIndex == pGraph.numberOfNeighbours(pPosition.outerIndex)
           && pPosition.outerIndex < pGraph.numberOfNodes()) {
      pPosition.innerIndex = 0;
      ++pPosition.outerIndex;
    }
    return *this;
  }

  bool operator!=(const AdjacencyListGraphIt<DIRECT>& other) const { return pPosition != other.pPosition; }

private:
  const AdjacencyListGraph<DIRECT>& pGraph;
  AdjListPos pPosition;
};

template <Directionality DIRECT>
inline AdjacencyListGraphIt<DIRECT> AdjacencyListGraph<DIRECT>::begin() const {
  return AdjacencyListGraphIt<DIRECT>(*this, AdjListPos{0, 0});
}

template <Directionality DIRECT>
inline AdjacencyListGraphIt<DIRECT> AdjacencyListGraph<DIRECT>::end() const {
  return AdjacencyListGraphIt(*this, AdjListPos{numberOfNodes(), 0});
}

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
    return Edge{pPosition.outerIndex, (size_t) (innerIndeces[pPosition.innerIndex])};
  }

  AdjMatGraphIt& operator++() override {
    assert(pGraph.matrix().isCompressed() && "check is only valid on compressed matrix");
    const int* outerIndeces = pGraph.matrix().outerIndexPtr();
    ++pPosition.innerIndex;
    while (pPosition.innerIndex >= (size_t) outerIndeces[pPosition.outerIndex + 1]) {
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
 *                                           graph algorithms
 **********************************************************************************************************************/

template <Directionality DIRECT>
DfsTree dfs(const AdjacencyMatrixGraph<DIRECT>& graph, const size_t rootNode = 0) {
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
      tree.explorationOrder().push_back(top);  // store order of node exploration
      for (Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>::InnerIterator it(graph.matrix(), top); it; ++it) {
        if (!visited[it.index()]) {
          nodeStack.push(it.index());     // do not push already visited nodes
          tree.parent(it.index()) = top;  // update the parent
        }
      }
    }
  }
  return tree;
}

struct EarDecomposition {
  std::vector<std::vector<size_t>> ears;
  std::vector<size_t> articulationPoints;

  bool open() const { return articulationPoints.empty(); }
};

EarDecomposition schmidt(const AdjacencyMatrixGraph<Directionality::Undirected>& graph);
AdjacencyMatrixGraph<Directionality::Undirected> earDecompToGraph(const EarDecomposition& earDecomposition);
AdjacencyMatrixGraph<Directionality::Undirected> biconnectedSpanningGraph(const Euclidean& euclidean);
