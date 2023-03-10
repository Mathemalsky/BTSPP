#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>  // can be removed when error handling is implemented
#include <vector>

#include <Eigen/SparseCore>

#include "graph/geometry.hpp"

#include "utility/utils.hpp"

/***********************************************************************************************************************
 *                                  forward declare functions from algorithm header
 **********************************************************************************************************************/

template <typename G>
bool checkBiconnectivity(const G& graph);

/***********************************************************************************************************************
 *                                                        Edges
 **********************************************************************************************************************/

/*!
 * \brief The Edge struct represents an edge from u to v in directed graphs or just between u and v in undirected graphs
 */
struct Edge {
  size_t u;
  size_t v;

  /*!
   * @brief edges with different outgoing node or different ingoing node are different
   * @details parallel edges are not different
   * @param other edge to compare
   * @return true if different, else otherwise
   */
  bool operator!=(const Edge& other) { return u != other.u || v != other.v; }

  /*!
   * \brief reverse creates a new add directed from v to u
   * \return edge in opposite direction
   */
  Edge reverse() const { return Edge{v, u}; }

  /*!
   * \brief invert swaps the internal storage of u and v
   */
  void invert() { std::swap(u, v); }
};

/*!
 * \brief The Directionality enum can be Directed or Undirected
 * \details This enum can be passed as template argument to some graph classes.
 */
enum class Directionality {
  Undirected,
  Directed
};

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
 *                                                  graph classes
 **********************************************************************************************************************/

/*!
 * \brief The Graph class is the base class of all more specialist types of graph. It's designed to contain the
 * functions graphs off all types have in common.
 */
class Graph {
public:
  virtual bool adjacent(const size_t u, const size_t v) const = 0;
  virtual bool adjacent(const Edge& e) const                  = 0;
  virtual bool connected() const                              = 0;
  virtual size_t numberOfEdges() const                        = 0;
  virtual size_t numberOfNodes() const                        = 0;
};

/*!
 * \brief The WeightedGraph class is an interface all graphs with weighted edges inherit from.
 */
class WeightedGraph : public virtual Graph {
public:
  virtual double weight(const size_t u, const size_t v) const = 0;
  virtual double weight(const Edge& e) const                  = 0;
};

/*!
 * \brief The CompleteGraph class is a inclomplete class for complete graphs, undirected simple graphs.
 * \details Complete Graph uses the property that every vertex is connected to every other vertex and hence edges do not
 * need to be stored explicitly. This class does not support directed graphs or multigraphs.
 */
class CompleteGraph : public virtual Graph {
public:
  bool adjacent([[maybe_unused]] const size_t u, [[maybe_unused]] const size_t v) const override { return true; }
  bool adjacent([[maybe_unused]] const Edge& e) const override { return true; }
  bool connected() const override { return true; }
  size_t numberOfEdges() const override { return numberOfNodes() * (numberOfNodes() - 1) / 2; }
};

class DirectedGraph : public virtual Graph {};

class UndirectedGraph : public virtual Graph {};

/*!
 * \brief The Euclidean class
 * \details An euclidean graph is a graph where each vertex has a position in 2 dimensional euclidean plane and the
 * weights of the edges are the euclidean distances between them.
 */
class Euclidean : public CompleteGraph, public WeightedGraph, public UndirectedGraph {
private:
  class Edges {
  private:
    class Iterator {
    public:
      struct Position {
        size_t index;
        const size_t numberOfNodes;
      };
      Iterator(const Position& pos) : pPosition(pos) {}

      Edge operator*() const {
        return Edge{pPosition.index / pPosition.numberOfNodes, pPosition.index % pPosition.numberOfNodes};
      }

      Iterator& operator++() {
        ++pPosition.index;
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition.index != other.pPosition.index; }

    private:
      Position pPosition;
    };  // end Iterator class

  public:
    Edges(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

    Iterator begin() const { return Iterator(Iterator::Position{0, pNumberOfNodes}); }
    Iterator end() const { return Iterator(Iterator::Position{pNumberOfNodes * pNumberOfNodes, pNumberOfNodes}); }

  private:
    const size_t pNumberOfNodes;
  };  // end Edges class

public:
  Euclidean()          = default;
  virtual ~Euclidean() = default;

  Euclidean(const std::vector<Point2D>& positions) : pPositions(positions) {}
  Euclidean(std::vector<Point2D>&& positions) : pPositions(positions) {}

  size_t numberOfNodes() const override { return pPositions.size(); }

  double weight(const size_t u, const size_t v) const override { return dist(pPositions[u], pPositions[v]); }
  double weight(const Edge& e) const override { return dist(pPositions[e.u], pPositions[e.v]); }

  Edges edges() const { return Edges(numberOfNodes()); }

  Point2D position(const size_t v) const { return pPositions[v]; }
  std::vector<Point2D>& verteces() { return pPositions; }

private:
  std::vector<Point2D> pPositions;
};

/*!
 * \brief The Modifyable class is an interface all modifyable graph classes inherit from.
 */
class Modifyable : public virtual Graph {
protected:
  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) = 0;
  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight)                     = 0;
};

/*!
 * \brief The AdjListGraph class is an abstract class for graphs which store adjacency as list.
 */
class AdjListGraph : public Modifyable {
private:
  class Nodes {
  private:
    class Iterator {
    public:
      Iterator(const size_t pos) : pPosition(pos) {}

      size_t operator*() const { return pPosition; }

      Iterator& operator++() {
        ++pPosition;
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition != other.pPosition; }

    private:
      size_t pPosition;
    };  // end Iterator class

  public:
    Nodes(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

    Iterator begin() const { return Iterator(0); }
    Iterator end() const { return Iterator(pNumberOfNodes); }

  private:
    const size_t pNumberOfNodes;
  };  // end Nodes class
public:
  AdjListGraph()  = default;
  ~AdjListGraph() = default;

  AdjListGraph(const AdjListGraph& graph) = default;
  AdjListGraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }
  AdjListGraph(const std::vector<std::vector<size_t>>& adjacencyList) : pAdjacencyList(adjacencyList) {}

  bool adjacent(const size_t u, const size_t v) const override {
    const std::vector<size_t>& neighbour         = pAdjacencyList[u];
    const std::vector<size_t>::const_iterator it = std::find(neighbour.begin(), neighbour.end(), v);
    return it != neighbour.end();
  }

  bool adjacent(const Edge& e) const override { return adjacent(e.u, e.v); }

  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  const std::vector<std::vector<size_t>>& adjacencyList() const { return pAdjacencyList; }

  /*!
   * \brief degree returns (out-) degree of node
   * \details in case of undirected graph, this function returns the outdegree
   * \param u node
   * \return u's degree
   */
  size_t degree(const size_t u) const { return pAdjacencyList[u].size(); }
  std::vector<size_t> degrees() const {
    std::vector<size_t> vec(numberOfNodes());
    for (size_t i = 0; i < numberOfNodes(); ++i) {
      vec[i] = degree(i);
    }
    return vec;
  }

  size_t neighbourAny(const size_t u) {
    assert(pAdjacencyList[u].size() > 0 && "There is no neighbour!");
    return pAdjacencyList[u][0];
  }

  template <typename func>
  size_t neighbourAnyExcept(const size_t u, func&& criteria) {
    for (const size_t v : pAdjacencyList[u]) {
      if (!criteria(v)) {
        return v;
      }
    }
    throw std::runtime_error("There is no node adjacent to " + std::to_string(u) + "not matching the criteria!");
  }

  template <typename func>
  size_t neighbourAnyPrefer(const size_t u, func&& criteria) {
    for (const size_t v : pAdjacencyList[u]) {
      if (criteria(v)) {
        return v;
      }
    }
    return neighbourAny(u);
  }

  const std::vector<size_t>& neighbours(const size_t u) const { return pAdjacencyList[u]; }

  Nodes nodes() const { return Nodes(pAdjacencyList.size()); }

  void removeAllEdges() {
    for (std::vector<size_t>& vec : pAdjacencyList) {
      vec.resize(0);
    }
  }

protected:
  std::vector<std::vector<size_t>> pAdjacencyList;
};

/*!
 * \brief The AdjacencyListGraph class implemnts undirected graphs based on adjacency lists.
 * \details For the sake of faster iteration over all neighbours, each edge is stored twice: once for each possible
 * direction.
 */
class AdjacencyListGraph : public AdjListGraph, public UndirectedGraph {
private:
  class Edges {
  private:
    class Iterator {
    public:
      struct AdjListPos {
        size_t outerIndex;
        size_t innerIndex;
      };

      Iterator(const std::vector<std::vector<size_t>>& adjacencyList, const AdjListPos& pos) :
        pAdjacencyList(adjacencyList),
        pPosition(pos) {}

      Edge operator*() const {
        return Edge{pPosition.outerIndex, pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex]};
      }

      Iterator& operator++() {
        assert(pPosition.outerIndex < pAdjacencyList.size() && "Iterator is already behind end!");

        ++pPosition.innerIndex;
        while (pPosition.outerIndex < pAdjacencyList.size() && (outOfNeighbours() || !toLowerIndex())) {
          if (outOfNeighbours()) {
            pPosition.innerIndex = 0;
            ++pPosition.outerIndex;
          }
          else {
            ++pPosition.innerIndex;
          }
        }
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition.outerIndex != other.pPosition.outerIndex; }

      bool toLowerIndex() const {
        return pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex] < pPosition.outerIndex;
      }

    private:
      bool outOfNeighbours() const { return pPosition.innerIndex >= pAdjacencyList[pPosition.outerIndex].size(); }

      const std::vector<std::vector<size_t>>& pAdjacencyList;
      AdjListPos pPosition;
    };  // end Iterator class

  public:
    Edges(const std::vector<std::vector<size_t>>& adjacencyList) : pAdjacencyList(adjacencyList) {}

    /*!
     * \brief begin returns the begin iterator for iteration
     * \details We start with 1 because for node 0 there is no node with strictly smaller index. If this inital iterator
     * does not point to an edge to lower index, it is increased to the next edge satisfying this condition.
     * \return Iterator infront of the first element.
     */
    Iterator begin() const {
      Iterator it(pAdjacencyList, Iterator::AdjListPos{1, 0});
      return it.toLowerIndex() ? it : ++it;
    }

    /*!
     * \brief end returns the end Iterator for comparison.
     * \details The 0 is just an arbitrary number because the outerIndex of AdjListPos is not taken into
     * account for comparison.
     * \return Iterator behind the last element.
     */
    Iterator end() const { return Iterator(pAdjacencyList, Iterator::AdjListPos{pAdjacencyList.size(), 0}); }

  private:
    const std::vector<std::vector<size_t>>& pAdjacencyList;
  };  // end Edges class

public:
  AdjacencyListGraph()  = default;
  ~AdjacencyListGraph() = default;

  AdjacencyListGraph(const AdjacencyListGraph& graph) = default;
  AdjacencyListGraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }
  AdjacencyListGraph(const std::vector<std::vector<size_t>>& adjacencyList) : AdjListGraph(adjacencyList) {}

  void addEdge(const size_t out, const size_t in, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    pAdjacencyList[out].push_back(in);
    pAdjacencyList[in].push_back(out);
  }

  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    pAdjacencyList[e.u].push_back(e.v);
    pAdjacencyList[e.v].push_back(e.u);
  }

  /*!
   * \brief AdjacencyListGraph::connected checks the graph is connected
   * \details Check if the connected componend containing vertex 0 is the whole graph, by performing a dfs.
   * \return true if the graph is connected, else false
   */
  bool connected() const override;

  size_t numberOfEdges() const override {
    return std::accumulate(pAdjacencyList.begin(), pAdjacencyList.end(), 0,
                           [](const unsigned int sum, const std::vector<size_t>& vec) { return sum + vec.size(); })
           / 2;
  }

  bool biconnected() const { return checkBiconnectivity(*this); }

  Edges edges() const { return Edges(pAdjacencyList); }

  void removeEdge(const Edge& e) {
    [[maybe_unused]] const bool removed = removeAnyElementByValue(pAdjacencyList[e.u], e.v);
    assert(removed && "Edge to be removed does not exist in graph!");
    [[maybe_unused]] const bool removed2 = removeAnyElementByValue(pAdjacencyList[e.v], e.u);
    assert(removed2 && "Edge to be removed does not exist in graph!");
  }

  void removeEdge(const size_t u, const size_t v) {
    [[maybe_unused]] const bool removed = removeAnyElementByValue(pAdjacencyList[u], v);
    assert(removed && "Edge to be removed does not exist in graph!");
    [[maybe_unused]] const bool removed2 = removeAnyElementByValue(pAdjacencyList[v], u);
    assert(removed2 && "Edge to be removed does not exist in graph!");
  }

  AdjacencyListGraph removeUncriticalEdges() const;
};

/*!
 * \brief The AdjacencyListDigraph class implements directed graphs based on adjacency lists.
 * \details The functions checking for connectivity are checking for connectivity in the sense of weak connectivity.
 */
class AdjacencyListDigraph : public AdjListGraph, public DirectedGraph {
private:
  class Edges {
  private:
    class Iterator {
    public:
      struct AdjListPos {
        size_t outerIndex;
        size_t innerIndex;
      };
      Iterator(const std::vector<std::vector<size_t>>& adjacencyList, const AdjListPos& pos) :
        pAdjacencyList(adjacencyList),
        pPosition(pos) {}

      Edge operator*() const {
        return Edge{pPosition.outerIndex, pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex]};
      }

      Iterator& operator++() {
        ++pPosition.innerIndex;
        while (pPosition.innerIndex >= pAdjacencyList[pPosition.outerIndex].size()
               && pPosition.outerIndex < pAdjacencyList.size()) {
          pPosition.innerIndex = 0;
          ++pPosition.outerIndex;
        }
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition.outerIndex != other.pPosition.outerIndex; }

    private:
      const std::vector<std::vector<size_t>>& pAdjacencyList;
      AdjListPos pPosition;
    };  // end Iterator class

  public:
    Edges(const std::vector<std::vector<size_t>>& adjacencyList) : pAdjacencyList(adjacencyList) {}

    Iterator begin() const { return Iterator(pAdjacencyList, Iterator::AdjListPos{0, 0}); }

    /*!
     * \brief end returns the end Iterator for comparison.
     * \details The 0 is just an arbitrary number because the outerIndex of AdjListPos is not taken into
     * account for comparison.
     * \return Iterator behind the last element.
     */
    Iterator end() const { return Iterator(pAdjacencyList, Iterator::AdjListPos{pAdjacencyList.size(), 0}); }

  private:
    const std::vector<std::vector<size_t>>& pAdjacencyList;
  };  // end Edges class

public:
  AdjacencyListDigraph()  = default;
  ~AdjacencyListDigraph() = default;

  AdjacencyListDigraph(const AdjacencyListDigraph& graph) = default;
  AdjacencyListDigraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }
  AdjacencyListDigraph(const std::vector<std::vector<size_t>>& adjacencyList) : AdjListGraph(adjacencyList) {}

  void addEdge(const size_t out, const size_t in, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    pAdjacencyList[out].push_back(in);
  }

  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override {
    pAdjacencyList[e.u].push_back(e.v);
  }

  /*!
   * \brief AdjacencyListDiGraph::connected checks the graph is connected
   * \details Check if the connected componend containing vertex 0 is the whole graph, by performing a dfs. Directed
   * graphs are considered as connected if the undirected graph obtained by adding an anti parallel edge for each edge
   * is connected.
   * \return true if the graph is connected, else false
   */
  bool connected() const override { return undirected().connected(); }

  size_t numberOfEdges() const override {
    return std::accumulate(pAdjacencyList.begin(), pAdjacencyList.end(), 0,
                           [](const unsigned int sum, const std::vector<size_t>& vec) { return sum + vec.size(); });
  }

  bool biconnected() const { return checkBiconnectivity(this->undirected()); }

  Edges edges() const { return Edges(pAdjacencyList); }

  void removeEdge(const Edge& e) {
    [[maybe_unused]] const bool removed = removeAnyElementByValue(pAdjacencyList[e.u], e.v);
    assert(removed && "Edge to be removed does not exist in graph!");
  }

  void removeEdge(const size_t u, const size_t v) {
    [[maybe_unused]] const bool removed = removeAnyElementByValue(pAdjacencyList[u], v);
    assert(removed && "Edge to be removed does not exist in graph!");
  }

  AdjacencyListDigraph removeUncriticalEdges() const;
  AdjacencyListGraph undirected() const;
};

/*!
 * \brief The AdjListGraph class is an abstract class for graphs which store adjacency as sparse matrix.
 * \details Graphs with adjacency matrix storage are generally assumed to be simple (without parallel edges) weighted
 * graphs.
 */
class AdjMatGraph : public Modifyable, public WeightedGraph {
private:
  class Neighbours {
  private:
    class Iterator {
    public:
      Iterator(const int* ptr) : pInnerIndeces(ptr) {}

      size_t operator*() const { return *pInnerIndeces; }

      Iterator operator++() {
        ++pInnerIndeces;
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pInnerIndeces != other.pInnerIndeces; }

    private:
      const int* pInnerIndeces;
    };  // end Iterator class

  public:
    Neighbours(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const size_t node) :
      pAdjacencyMatrix(adjacencyMatrix),
      pNode(node) {
      assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
    }

    Iterator begin() const {
      return Iterator(pAdjacencyMatrix.innerIndexPtr() + pAdjacencyMatrix.outerIndexPtr()[pNode]);
    }

    Iterator end() const {
      return Iterator(pAdjacencyMatrix.innerIndexPtr() + pAdjacencyMatrix.outerIndexPtr()[pNode + 1]);
    }

  private:
    const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix;
    const size_t pNode;
  };  // end Neighbours class

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
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != 0 && "Edgeweight 0 cann also mean the edge does not exist!");
    return pAdjacencyMatrix.coeff(e.u, e.v).cost();
  }
  double weight(const size_t u, const size_t v) const override {
    assert(pAdjacencyMatrix.coeff(u, v) != 0 && "Edgeweight 0 cann also mean the edge does not exist!");
    return pAdjacencyMatrix.coeff(u, v).cost();
  }

  size_t numberOfEdges() const override { return pAdjacencyMatrix.nonZeros(); }
  size_t numberOfNodes() const override { return pAdjacencyMatrix.cols(); }

  void compressMatrix() { pAdjacencyMatrix.makeCompressed(); }

  const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& matrix() const { return pAdjacencyMatrix; }

  Neighbours neighbours(const size_t node) const { return Neighbours(pAdjacencyMatrix, node); }

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

/*!
 * \brief The AdjacencyMatrixGraph class implements the concrete classes for undirected graphs.
 * \details For the sake of faster iteration over all neighbours, each edge is stored twice: once for each possible
 * direction.
 */
class AdjacencyMatrixGraph : public AdjMatGraph, public UndirectedGraph {
private:
  class Edges {
  private:
    class Iterator {
    public:
      struct SparseMatrixPos {
        size_t innerIndex;
        size_t outerIndex;
      };

      Iterator(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const SparseMatrixPos& pos) :
        pAdjacencyMatrix(adjacencyMatrix),
        pPosition(pos) {
        assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
      }

      Edge operator*() const {
        const int* const innerIndeces = pAdjacencyMatrix.innerIndexPtr();
        return Edge{pPosition.outerIndex, static_cast<size_t>(innerIndeces[pPosition.innerIndex])};
      }

      Iterator& operator++() {
        const int* const outerIndeces = pAdjacencyMatrix.outerIndexPtr();
        const int* const innerIndeces = pAdjacencyMatrix.innerIndexPtr();
        ++pPosition.innerIndex;
        while (pPosition.innerIndex < static_cast<size_t>(pAdjacencyMatrix.nonZeros()) && !valid()) {
          if (static_cast<size_t>(innerIndeces[pPosition.innerIndex]) >= pPosition.outerIndex) {
            pPosition.innerIndex = outerIndeces[pPosition.outerIndex + 1];  // skip rest of the row
          }
          ++pPosition.outerIndex;  // goes to next row
        }
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition.innerIndex != other.pPosition.innerIndex; }

      bool valid() const {
        const int* const outerIndeces = pAdjacencyMatrix.outerIndexPtr();
        const int* const innerIndeces = pAdjacencyMatrix.innerIndexPtr();
        return pPosition.innerIndex < static_cast<size_t>(outerIndeces[pPosition.outerIndex + 1])  // correct row
               && static_cast<size_t>(innerIndeces[pPosition.innerIndex]) < pPosition.outerIndex;  // lower triangular
      }

    private:
      const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix;
      SparseMatrixPos pPosition;
    };  // end Iterator class

  public:
    Edges(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix) :
      pAdjacencyMatrix(adjacencyMatrix) {}

    /*!
     * \brief begin returns a an Iterator to the begin of the strictly lower triangle in adjacency matrix.
     * \details The Iterator starts in second row, because the intersect of first row and strictly lower triangle is
     * empty.
     * \return Iterator to the begin of the strictly lower triangle in adjacency matrix
     */
    Iterator begin() const {
      assert(pAdjacencyMatrix.rows() > 1 && "Adjacency matrix with less than 2 rows has no entries below diagonal!");
      Iterator it(pAdjacencyMatrix,
                  Iterator::SparseMatrixPos{static_cast<size_t>(pAdjacencyMatrix.outerIndexPtr()[1]), 0});
      return it.valid() ? it : ++it;
    }

    /*!
     * \brief end returns the end Iterator for comparison.
     * \details The 0 is just an arbitrary number because the outerIndex of SparseMatrixPos is not taken into
     * account for comparison.
     * \return Iterator behind the last element.
     */
    Iterator end() const {
      return Iterator(pAdjacencyMatrix, Iterator::SparseMatrixPos{static_cast<size_t>(pAdjacencyMatrix.nonZeros()), 0});
    }

  private:
    const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix;
  };  // end Edges class

public:
  AdjacencyMatrixGraph()  = default;
  ~AdjacencyMatrixGraph() = default;

  AdjacencyMatrixGraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : AdjMatGraph(mat) {}
  AdjacencyMatrixGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {}

  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(out, in, edgeWeight);
    AdjMatGraph::addEdge(in, out, edgeWeight);
  }

  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { addEdge(e.u, e.v, edgeWeight); }

  bool connected() const override;

  bool biconnected() const { return checkBiconnectivity(*this); }

  Edges edges() const { return Edges(pAdjacencyMatrix); }

  void removeEdge(const Edge& e) {
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != 0 && "Edge to be removed does not exist in graph!");
    pAdjacencyMatrix.coeffRef(e.u, e.v) = 0.0;
    assert(pAdjacencyMatrix.coeff(e.v, e.u) != 0 && "Edge to be removed does not exist in graph!");
    pAdjacencyMatrix.coeffRef(e.v, e.u) = 0.0;
  }
  AdjacencyMatrixGraph removeUncriticalEdges() const;
};

/*!
 * \brief The AdjacencyMatrixGraph class implements the concrete class for directed graphs.
 * \details The functions checking for connectivity are checking for connectivity in the sense of weak connectivity.
 * Digraphs may have antiparallel edges.
 */
class AdjacencyMatrixDigraph : public AdjMatGraph, public DirectedGraph {
private:
  class Edges {
  private:
    class Iterator {
    public:
      struct SparseMatrixPos {
        size_t innerIndex;
        size_t outerIndex;
      };

      Iterator(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const SparseMatrixPos& pos) :
        pAdjacencyMatrix(adjacencyMatrix),
        pPosition(pos) {
        assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
      }

      Edge operator*() const {
        const int* innerIndeces = pAdjacencyMatrix.innerIndexPtr();
        return Edge{pPosition.outerIndex, static_cast<size_t>(innerIndeces[pPosition.innerIndex])};
      }

      Iterator& operator++() {
        const int* outerIndeces = pAdjacencyMatrix.outerIndexPtr();
        ++pPosition.innerIndex;
        while (pPosition.innerIndex >= static_cast<size_t>(outerIndeces[pPosition.outerIndex + 1])) {
          ++pPosition.outerIndex;
        }
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition.innerIndex != other.pPosition.innerIndex; }

    private:
      const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix;
      SparseMatrixPos pPosition;
    };  // end Iterator class

  public:
    Edges(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix) :
      pAdjacencyMatrix(adjacencyMatrix) {}

    Iterator begin() const { return Iterator(pAdjacencyMatrix, Iterator::SparseMatrixPos{0, 0}); }

    /*!
     * \brief end returns the end Iterator for comparison.
     * \details The 0 is just an arbitrary number because the outerIndex of SparseMatrixPos is not taken into
     * account for comparison.
     * \return Iterator behind the last element.
     */
    Iterator end() const {
      return Iterator(pAdjacencyMatrix, Iterator::SparseMatrixPos{static_cast<size_t>(pAdjacencyMatrix.nonZeros()), 0});
    }

  private:
    const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix;
  };  // end Edges class

public:
  AdjacencyMatrixDigraph()  = default;
  ~AdjacencyMatrixDigraph() = default;

  AdjacencyMatrixDigraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : AdjMatGraph(mat) {}
  AdjacencyMatrixDigraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {}

  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(out, in, edgeWeight);
  }

  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { addEdge(e.u, e.v, edgeWeight); }

  bool connected() const override { return undirected().connected(); };

  bool biconnected() const { return checkBiconnectivity(undirected()); }

  Edges edges() const { return Edges(pAdjacencyMatrix); }

  void removeEdge(const Edge& e) {
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != 0 && "Edge to be removed does not exist in graph!");
    pAdjacencyMatrix.coeffRef(e.u, e.v) = 0.0;
  }
  AdjacencyMatrixDigraph removeUncriticalEdges() const;
  AdjacencyMatrixGraph undirected() const;
};

/*!
 * \brief The Tree class is an abstract which implements the methods all trees have in common.
 */
class Tree : public virtual Graph {
  bool connected() const override { return true; }
  size_t numberOfEdges() const override { return numberOfNodes() - 1; }
};

/*!
 * \brief class DfsTree
 * \details This class is for trees where every node except for the root has exactly one parent and an edge directed to
 * it's parent. This allows to store the neighboors more efficient.
 */
class DfsTree : public Tree, public DirectedGraph {
private:
  class Edges {
  private:
    class Iterator {
    public:
      Iterator(const std::vector<size_t>& adjacencyList, const size_t root, const size_t pos) :
        pAdjacencyList(adjacencyList),
        pRoot(root),
        pPosition(pos) {}

      Edge operator*() const {
        assert(pPosition != pRoot && "No outging edge from root node!");
        return Edge{pPosition, pAdjacencyList[pPosition]};
      }

      Iterator& operator++() {
        ++pPosition;
        if (pPosition == pRoot) {
          ++pPosition;
        }
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition != other.pPosition; }

    private:
      const std::vector<size_t>& pAdjacencyList;
      const size_t pRoot;
      size_t pPosition;
    };  // end Iterator class

  public:
    Edges(const std::vector<size_t>& adjacencyList, const size_t root) : pAdjacencyList(adjacencyList), pRoot(root) {}

    Iterator begin() const {
      if (pRoot != 0) {
        return Iterator(pAdjacencyList, pRoot, 0);  // skip the root node
      }
      else {
        return Iterator(pAdjacencyList, pRoot, 1);
      }
    }

    Iterator end() const { return Iterator(pAdjacencyList, pRoot, pAdjacencyList.size()); }

  private:
    const std::vector<size_t>& pAdjacencyList;
    const size_t pRoot;
  };  // end Edges class

public:
  DfsTree()  = default;
  ~DfsTree() = default;

  DfsTree(const size_t numberOfNodes) : pAdjacencyList(numberOfNodes) {
    pExplorationOrder.reserve(numberOfNodes);  // just reserve, because dfs performs push_backs
  }

  bool adjacent(const size_t u, const size_t v) const override { return u != 0 && v == parent(u); }
  bool adjacent(const Edge& e) const override { return e.u != 0 && e.v == parent(e.u); }

  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  Edges edges() const { return Edges(pAdjacencyList, root()); }

  std::vector<size_t> explorationOrder() const { return pExplorationOrder; }
  std::vector<size_t>& explorationOrder() { return pExplorationOrder; }

  size_t parent(const size_t u) const {
    assert(u != root() && "You are trying to access the root nodes parent, which is uninitialized memory!");
    return pAdjacencyList[u];
  }

  size_t& parent(const size_t u) {
    assert(u != root() && "You are trying to access the root nodes parent, which is uninitialized memory!");
    return pAdjacencyList[u];
  }

  size_t root() const { return pExplorationOrder[0]; }

private:
  std::vector<size_t> pAdjacencyList;
  std::vector<size_t> pExplorationOrder;
};
