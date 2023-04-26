#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>  // can be removed when error handling is implemented
#include <vector>

#include <Eigen/SparseCore>

#include "exception/exceptions.hpp"

#include "graph/geometry.hpp"
#include "graph/utils.hpp"

namespace graph {

/***********************************************************************************************************************
 *                                  forward declare functions from algorithm header
 **********************************************************************************************************************/

/*!
 * @brief checks if the graph is biconnected
 * @tparam G type of graph
 * @param graph
 * @return true if the graph is biconnected
 */
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
  bool operator!=(const Edge& other) const { return u != other.u || v != other.v; }

  /*!
   * @brief edges with same outgoing node and same ingoing node are identic
   * @param other edge to compare
   * @return bool
   */
  bool operator==(const Edge& other) const { return u == other.u && v == other.v; }

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
 */
enum class Directionality {
  Undirected,
  Directed
};

class EdgeWeight {
public:
  EdgeWeight()  = default;
  ~EdgeWeight() = default;

  /*!
   * @brief constructs EdgeWeight with value weight
   */
  EdgeWeight(const double weight) : pWeight(weight) {}

  /*!
   * @brief returns the weight as double
   */
  double weight() const { return pWeight; }

  /*!
   * @brief assignment operator from double
   */
  void operator=(const double weight) { pWeight = weight; }

  /*!
   * @brief compares for inquality with double
   */
  bool operator==(const double compare) { return pWeight == compare; }

  /*!
   * @brief compares for inequality with double
   */
  bool operator!=(const double compare) { return pWeight != compare; }

  /*!
   * @brief compares for \leq with double
   * @details This function is needed in Eigen.
   */
  bool operator<=(const EdgeWeight other) { return pWeight <= other.pWeight; }

  /*!
   * @brief operator+ implements the minimum of the edgeweights
   * @details This function is needed in Eigen.
   */
  EdgeWeight operator+(const EdgeWeight other) const { return EdgeWeight(std::min(pWeight, other.pWeight)); }

  /*!
   * @brief operator* implements the sum of the edgeweights
   * @details This function is needed in Eigen.
   */
  EdgeWeight operator*(const EdgeWeight other) const { return EdgeWeight(pWeight + other.pWeight); }

  /*!
   * @brief operator+= setsto  the minimum of the edgeweights
   * @details This function is needed in Eigen.
   */
  EdgeWeight& operator+=(const EdgeWeight other) {
    pWeight = std::min(pWeight, other.pWeight);
    return *this;
  }

  /*!
   * @brief returns absolut value
   * @details This function is needed in Eigen.
   */
  friend EdgeWeight abs(const EdgeWeight weight) { return EdgeWeight(std::abs(weight.pWeight)); }

private:
  double pWeight; /*!< stores the weight of the edge*/
};
}  // namespace graph

/*
 * Eigen needs some hints to deal with the custom type.
 */
namespace Eigen {
template <>
struct NumTraits<graph::EdgeWeight> : GenericNumTraits<graph::EdgeWeight> {
  typedef graph::EdgeWeight Real;
  typedef graph::EdgeWeight NonInteger;
  typedef graph::EdgeWeight Nested;

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

namespace graph {

/***********************************************************************************************************************
 *                                                  graph classes
 **********************************************************************************************************************/

/*!
 * \brief The Graph class is the base class of all more specialist types of graph. It's designed to contain the
 * functions graphs off all types have in common.
 */
class Graph {
public:
  /*!
   * @brief checks if there is an edge (u,v) in the graph
   * @param u outgoing vertex
   * @param v incomming vertex
   */
  virtual bool adjacent(const size_t u, const size_t v) const = 0;

  /*!
   * @brief checks if there is an edge e in the graph
   * @param e edge to check
   */
  virtual bool adjacent(const Edge& e) const = 0;

  /*!
   * @brief checks if the graph is connected
   */
  virtual bool connected() const = 0;

  /*!
   * @brief returns the number of edges in the graph
   */
  virtual size_t numberOfEdges() const = 0;

  /*!
   * @brief returns the number of nodes in the graph
   */
  virtual size_t numberOfNodes() const = 0;
};

/*!
 * \brief The WeightedGraph class is an interface all graphs with weighted edges inherit from.
 */
class WeightedGraph : public virtual Graph {
public:
  /*!
   * @brief returns the edge's weight
   * @param u outgoing vertex
   * @param v incomming vertex
   */
  virtual double weight(const size_t u, const size_t v) const = 0;

  /*!
   * @brief returns the edge's weight
   * @param e edge
   */
  virtual double weight(const Edge& e) const = 0;
};

/*!
 * \brief The CompleteGraph class is a inclomplete class for complete graphs, undirected simple graphs.
 * \details Complete Graph uses the property that every vertex is connected to every other vertex and hence edges do not
 * need to be stored explicitly. This class does not support directed graphs or multigraphs.
 */
class CompleteGraph : public virtual Graph {
public:
  /*!
   * @brief returns true because it's a complete graph and every node is connected to every node
   */
  bool adjacent([[maybe_unused]] const size_t u, [[maybe_unused]] const size_t v) const override { return true; }

  /*!
   * @brief returns true because it's a complete graph and every node is connected to every node
   */
  bool adjacent([[maybe_unused]] const Edge& e) const override { return true; }

  /*!
   * @brief returns true because it's a complete graph and every node is connected to every node
   */
  bool connected() const override { return true; }

  /*!
   * @brief returns true because it's a complete graph and every node is connected to every node
   */
  size_t numberOfEdges() const override { return numberOfNodes() * (numberOfNodes() - 1) / 2; }
};

/*!
 * @brief class DirectedGraph is a virtual base class for all directed graphs
 */
class DirectedGraph : public virtual Graph {};

/*!
 * @brief class UndirectedGraph is a virtual base class for all directed graphs
 */
class UndirectedGraph : public virtual Graph {};

/*!
 * \brief The Euclidean class
 * \details An euclidean graph is a graph where each vertex has a position in 2 dimensional euclidean plane and the
 * weights of the edges are the euclidean distances between them.
 */
class Euclidean : public CompleteGraph, public WeightedGraph, public UndirectedGraph {
private:
  /*!
   * @brief Edges is a facade class to iterate over all edges in euclidean graph.
   */
  class Edges {
  private:
    /*!
     * @brief Iterator on edges of an Euclidean graph
     */
    class Iterator {
    public:
      /*!
       * @brief wraps position variables
       */
      struct Position {
        size_t index;               /**< index of edge, edge is (index / numberOfNodes, index % numberOfNodes)*/
        const size_t numberOfNodes; /**< number of nodes in euclidean graph*/
      };

      /*!
       * @brief constructs iterator at position pos
       */
      Iterator(const Position& pos) : pPosition(pos) {}

      /*!
       * @brief constructs edge from current position as explicit object
       * @return edge from current position as explicit object
       */
      Edge operator*() const { return Edge{pPosition.index / pPosition.numberOfNodes, pPosition.index % pPosition.numberOfNodes}; }

      /*!
       * @brief increases index by one
       * @return reference to this iterator after increasing index
       */
      Iterator& operator++() {
        ++pPosition.index;
        return *this;
      }

      /*!
       * @brief compares for inequality
       */
      bool operator!=(const Iterator& other) const { return pPosition.index != other.pPosition.index; }

    private:
      Position pPosition; /**< position of index in graph*/
    };                    // end Iterator class

  public:
    /*!
     * @brief creates edges instance of edges
     * @param numberOfNodes needed to compute the index correct
     */
    Edges(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

    /*!
     * @brief creates an iterator pointing in front of the first edge
     * @return begin iterator
     */
    Iterator begin() const { return Iterator(Iterator::Position{0, pNumberOfNodes}); }

    /*!
     * @brief creates an iterator pointing in behind the last edge
     * @return end iterator
     */
    Iterator end() const { return Iterator(Iterator::Position{pNumberOfNodes * pNumberOfNodes, pNumberOfNodes}); }

  private:
    const size_t pNumberOfNodes; /**< number of nodes in euclidean graph */
  };                             // end Edges class

public:
  Euclidean()          = default;
  virtual ~Euclidean() = default;

  /*!
   * @brief constructor from vector of Point2D
   */
  Euclidean(const std::vector<Point2D>& positions) : pPositions(positions) {}

  /*!
   * @brief move constructor from vector of Point2D
   */
  Euclidean(std::vector<Point2D>&& positions) : pPositions(positions) {}

  /*!
   * @brief returns number of points in internal vector as number of nodes
   */
  size_t numberOfNodes() const override { return pPositions.size(); }

  /*!
   * @brief returns euclidean distance between nodes
   * @param u start node
   * @param v end node
   */
  double weight(const size_t u, const size_t v) const override { return dist(pPositions[u], pPositions[v]); }

  /*!
   * @brief returns euclidean distance between nodes
   * @param e edge for which weight is requested
   */
  double weight(const Edge& e) const override { return dist(pPositions[e.u], pPositions[e.v]); }

  /*!
   * @brief creates edges object to iterate of the edges
   */
  Edges edges() const { return Edges(numberOfNodes()); }

  /*!
   * @brief position of node in plane
   * @param v node
   * @return position of node v
   */
  Point2D position(const size_t v) const { return pPositions[v]; }

  /*!
   * @brief returns reference to internal Point2D vector
   */
  std::vector<Point2D>& vertices() { return pPositions; }

private:
  std::vector<Point2D> pPositions; /**< vector of Point2D holding the nodes position in plane*/
};

/*!
 * \brief The Modifyable class is an interface all modifyable graph classes inherit from.
 */
class Modifyable : public virtual Graph {
protected:
  /*!
   * @brief adds an edge from out to in with weight edgeWeight
   * @param out start vertex of edge
   * @param end vertex of edge
   * @param edgeWeight of new edge
   */
  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) = 0;

  /*!
   * @brief adds an edge from out to in with weight edgeWeight
   * @param e edge to be added
   * @param edgeWeight of new edge
   */
  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight) = 0;
};

/*!
 * @brief The AdjListGraph class is an abstract class for graphs which store adjacency as list.
 * @details This class can handle multigraphs, but does not differentiate between parallel edges.
 */
class AdjListGraph : public Modifyable {
private:
  /*!
   * @brief Nodes is a facade class to iterate over all nodes in AdjListGraph.
   */
  class Nodes {
  private:
    /*!
     * @brief Iterator on nodes of an AdjListGraph graph
     */
    class Iterator {
    public:
      /*!
       * @brief creates Iterator at position pos.
       */
      Iterator(const size_t pos) : pPosition(pos) {}

      /*!
       * @brief returns node at current position
       */
      size_t operator*() const { return pPosition; }

      /*!
       * @brief moves iterator one node forward
       * @return reference to this iterator after incrementation
       */
      Iterator& operator++() {
        ++pPosition;
        return *this;
      }

      /*!
       * @brief compares for inequality
       */
      bool operator!=(const Iterator& other) const { return pPosition != other.pPosition; }

    private:
      size_t pPosition; /**< position of index in nodes*/
    };                  // end Iterator class

  public:
    /*!
     * constructor for for nodes facade
     */
    Nodes(const size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

    /*!
     * @brief creates an iterator pointing in front of the first node
     * @return begin iterator
     */
    Iterator begin() const { return Iterator(0); }

    /*!
     * @brief creates an iterator pointing in behind the last node
     * @return end iterator
     */
    Iterator end() const { return Iterator(pNumberOfNodes); }

  private:
    const size_t pNumberOfNodes; /**< number of nodes in adjacency list graph */
  };                             // end Nodes class
public:
  AdjListGraph()  = default;
  ~AdjListGraph() = default;

  AdjListGraph(const AdjListGraph& graph) = default;

  /*!
   * @brief constructor, resizes the internal adjacency list to numberOfNodes
   * @param numberOfNodes number of nodes in constructed graph
   */
  AdjListGraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }

  /*!
   * @brief constructor, sets the internal adjacency list to adjacencyList
   * @param adjacencyList adjacency list
   */
  AdjListGraph(const std::vector<std::vector<size_t>>& adjacencyList) : pAdjacencyList(adjacencyList) {}

  /*!
   * @brief checks if v appears in the list of neighbours of u
   * @param u node
   * @param v node
   */
  bool adjacent(const size_t u, const size_t v) const override {
    const std::vector<size_t>& neighbour         = pAdjacencyList[u];
    const std::vector<size_t>::const_iterator it = std::find(neighbour.begin(), neighbour.end(), v);
    return it != neighbour.end();
  }

  /*!
   * @brief checks if there is an edge e in the graph
   * @param e edge to check
   */
  bool adjacent(const Edge& e) const override { return adjacent(e.u, e.v); }

  /*!
   * @brief returns the number of nodes in adjacency list
   */
  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  /*!
   * @brief gives complete adjacency list
   * @details The i^th entry of the outer vectors contains the indices of all nodes reached with edges outgoing from node i.
   * @return vector of vector of size_t (indices)
   */
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

  /*!
   * @brief returns a neighbour of u
   * @param u node
   * @return first neighbour of u in the adjacency list
   */
  size_t neighbourAny(const size_t u) {
    assert(pAdjacencyList[u].size() > 0 && "There is no neighbour!");
    return pAdjacencyList[u][0];
  }

  /*!
   * @brief returns a neighbour of u not satisfying the criteria
   * @tparam func type of lambda function
   * @param u node
   * @param criteria lambda function implemnting the criteria
   * @return first neighbour not matching the criteria
   */
  template <typename func>
  size_t neighbourAnyExcept(const size_t u, func&& criteria) {
    for (const size_t v : pAdjacencyList[u]) {
      if (!criteria(v)) {
        return v;
      }
    }
    throw InfesableRequest("There is no node adjacent to " + std::to_string(u) + "not matching the criteria!");
  }

  /*!
   * @brief returns a neighbour, if possible satisfying the criteria
   * @tparam func type of lambda function
   * @param u node
   * @param criteria lambda function implemnting the criteria
   * @return first neighbour matching the criteria, if there is no then the fist neigbour
   */
  template <typename func>
  size_t neighbourAnyPrefer(const size_t u, func&& criteria) {
    for (const size_t v : pAdjacencyList[u]) {
      if (criteria(v)) {
        return v;
      }
    }
    return neighbourAny(u);
  }

  /*!
   * @brief returns list of neighbours
   * @param u node
   * @return const refference to vector in adjacency list
   */
  const std::vector<size_t>& neighbours(const size_t u) const { return pAdjacencyList[u]; }

  /*!
   * @brief creates iteratable instance of Nodes
   * @return instance of Nodes with the number of nodes in adjacency list
   */
  Nodes nodes() const { return Nodes(pAdjacencyList.size()); }

  /*!
   * @brief deletes all edges in the graph
   */
  void removeAllEdges() {
    for (std::vector<size_t>& vec : pAdjacencyList) {
      vec.resize(0);
    }
  }

protected:
  std::vector<std::vector<size_t>> pAdjacencyList; /**< i-th vector contains neighbours of i */
};

/*!
 * \brief The AdjacencyListGraph class implemnts undirected graphs based on adjacency lists.
 * \details For the sake of faster iteration over all neighbours, each edge is stored twice: once for each possible
 * direction.
 */
class AdjacencyListGraph : public AdjListGraph, public UndirectedGraph {
private:
  /*!
   * @brief Edges is a facade class to iterate over all edges in AdjacencyListGraph graph.
   */
  class Edges {
  private:
    /*!
     * @brief Iterator on edges of an AdjacencyListGraph graph
     */
    class Iterator {
    public:
      /*!
       * @brief AdjListPos bundles inner and outer index of the position
       */
      struct AdjListPos {
        size_t outerIndex; /**< outer index (out node) */
        size_t innerIndex; /**< inner index (in node) */
      };

      /*!
       * @brief creates an iterator pointing at pos
       * @param adjacencyList reference to graphs adjacency list
       * @param pos contains position as pair (outer, inner)
       */
      Iterator(const std::vector<std::vector<size_t>>& adjacencyList, const AdjListPos& pos) :
        pAdjacencyList(adjacencyList),
        pPosition(pos) {}

      /*!
       * @brief constructs edge from current position as explicit object
       * @return edge from current position as explicit object
       */
      Edge operator*() const { return Edge{pPosition.outerIndex, pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex]}; }

      /*!
       * @brief moves the operator one position forward
       * @details Increases inner index by 1. Repeats that until a valid (outer index, inner index) pair is found, or iterartor is out of
       * bounds. Valid index pair has outer index > inner index and points to an element in adjacencyList
       * @return reference to passed iterator
       */
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

      /*!
       * @brief compares for inequality
       */
      bool operator!=(const Iterator& other) const { return pPosition.outerIndex != other.pPosition.outerIndex; }

      /*!
       * @brief checks if inner index < outer index
       */
      bool toLowerIndex() const { return pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex] < pPosition.outerIndex; }

    private:
      /*!
       * @brief checks if the inner index is still in the range of neighbours vector
       */
      bool outOfNeighbours() const { return pPosition.innerIndex >= pAdjacencyList[pPosition.outerIndex].size(); }

      const std::vector<std::vector<size_t>>& pAdjacencyList; /**< reference to graphs adjacency list */
      AdjListPos pPosition;                                   /**< position in the adjacency list */
    };                                                        // end Iterator class

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

  /*!
   * @brief counts the number of edges
   * @return size_t
   */
  size_t numberOfEdges() const override {
    return std::accumulate(pAdjacencyList.begin(),
                           pAdjacencyList.end(),
                           0,
                           [](const unsigned int sum, const std::vector<size_t>& vec) { return sum + vec.size(); })
           / 2;
  }

  /*!
   * @brief checks if the graph is biconnected
   * @return bool
   */
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

      Edge operator*() const { return Edge{pPosition.outerIndex, pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex]}; }

      Iterator& operator++() {
        ++pPosition.innerIndex;
        while (pPosition.innerIndex >= pAdjacencyList[pPosition.outerIndex].size() && pPosition.outerIndex < pAdjacencyList.size()) {
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

  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = 1.0) override { pAdjacencyList[e.u].push_back(e.v); }

  /*!
   * \brief AdjacencyListDiGraph::connected checks the graph is connected
   * \details Check if the connected componend containing vertex 0 is the whole graph, by performing a dfs. Directed
   * graphs are considered as connected if the undirected graph obtained by adding an anti parallel edge for each edge
   * is connected.
   * \return true if the graph is connected, else false
   */
  bool connected() const override { return undirected().connected(); }

  size_t numberOfEdges() const override {
    return std::accumulate(pAdjacencyList.begin(), pAdjacencyList.end(), 0, [](const unsigned int sum, const std::vector<size_t>& vec) {
      return sum + vec.size();
    });
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

  /*!
   * @brief creates an undirected graph with the same nodes and edges
   * @details antiparallel edge are wedded
   * @return Adjacency list graph
   */
  AdjacencyListGraph undirected() const;
};

/*!
 * \brief.The AdjListGraph class is an abstract class for graphs which store adjacency as sparse matrix.
 * \details Graphs with adjacency matrix storage are generally assumed to be simple (without parallel edges) weighted
 * graphs.
 */
class AdjMatGraph : public Modifyable, public WeightedGraph {
private:
  class Neighbours {
  private:
    class Iterator {
    public:
      Iterator(const int* ptr) : pInnerIndices(ptr) {}

      size_t operator*() const { return *pInnerIndices; }

      Iterator operator++() {
        ++pInnerIndices;
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pInnerIndices != other.pInnerIndices; }

    private:
      const int* pInnerIndices;
    };  // end Iterator class

  public:
    Neighbours(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const size_t node) :
      pAdjacencyMatrix(adjacencyMatrix),
      pNode(node) {
      assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
    }

    Iterator begin() const { return Iterator(pAdjacencyMatrix.innerIndexPtr() + pAdjacencyMatrix.outerIndexPtr()[pNode]); }

    Iterator end() const { return Iterator(pAdjacencyMatrix.innerIndexPtr() + pAdjacencyMatrix.outerIndexPtr()[pNode + 1]); }

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
    return pAdjacencyMatrix.coeff(e.u, e.v).weight();
  }
  double weight(const size_t u, const size_t v) const override {
    assert(pAdjacencyMatrix.coeff(u, v) != 0 && "Edgeweight 0 cann also mean the edge does not exist!");
    return pAdjacencyMatrix.coeff(u, v).weight();
  }

  size_t numberOfEdges() const override { return pAdjacencyMatrix.nonZeros(); }
  size_t numberOfNodes() const override { return pAdjacencyMatrix.cols(); }

  void compressMatrix() { pAdjacencyMatrix.makeCompressed(); }

  const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& matrix() const { return pAdjacencyMatrix; }

  Neighbours neighbours(const size_t node) const { return Neighbours(pAdjacencyMatrix, node); }

  void pruneMatrix() { pAdjacencyMatrix.prune(0.0, Eigen::NumTraits<EdgeWeight>::dummy_precision()); }

  void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }  // operator*= not provided by Eigen

protected:
  Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor> pAdjacencyMatrix;

  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    pAdjacencyMatrix.insert(out, in) = edgeWeight;
  }

  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { pAdjacencyMatrix.insert(e.u, e.v) = edgeWeight; }
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
        const int* const innerIndices = pAdjacencyMatrix.innerIndexPtr();
        return Edge{pPosition.outerIndex, static_cast<size_t>(innerIndices[pPosition.innerIndex])};
      }

      Iterator& operator++() {
        const int* const outerIndices = pAdjacencyMatrix.outerIndexPtr();
        const int* const innerIndices = pAdjacencyMatrix.innerIndexPtr();
        ++pPosition.innerIndex;
        while (pPosition.innerIndex < static_cast<size_t>(pAdjacencyMatrix.nonZeros()) && !valid()) {
          if (static_cast<size_t>(innerIndices[pPosition.innerIndex]) >= pPosition.outerIndex) {
            pPosition.innerIndex = outerIndices[pPosition.outerIndex + 1];  // skip rest of the row
          }
          ++pPosition.outerIndex;  // goes to next row
        }
        return *this;
      }

      bool operator!=(const Iterator& other) const { return pPosition.innerIndex != other.pPosition.innerIndex; }

      bool valid() const {
        const int* const outerIndices = pAdjacencyMatrix.outerIndexPtr();
        const int* const innerIndices = pAdjacencyMatrix.innerIndexPtr();
        return pPosition.innerIndex < static_cast<size_t>(outerIndices[pPosition.outerIndex + 1])  // correct row
               && static_cast<size_t>(innerIndices[pPosition.innerIndex]) < pPosition.outerIndex;  // lower triangular
      }

    private:
      const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix;
      SparseMatrixPos pPosition;
    };  // end Iterator class

  public:
    Edges(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix) : pAdjacencyMatrix(adjacencyMatrix) {}

    /*!
     * \brief begin returns a an Iterator to the begin of the strictly lower triangle in adjacency matrix.
     * \details The Iterator starts in second row, because the intersect of first row and strictly lower triangle is
     * empty.
     * \return Iterator to the begin of the strictly lower triangle in adjacency matrix
     */
    Iterator begin() const {
      assert(pAdjacencyMatrix.rows() > 1 && "Adjacency matrix with less than 2 rows has no entries below diagonal!");
      Iterator it(pAdjacencyMatrix, Iterator::SparseMatrixPos{static_cast<size_t>(pAdjacencyMatrix.outerIndexPtr()[1]), 0});
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
};

/*!
 * \brief The AdjacencyMatrixGraph class implements the concrete class for directed graphs.
 * \details The functions checking for connectivity are checking for connectivity in the sense of weak connectivity.
 * AdjacencyMatrixDigraph may have antiparallel edges.
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
        const int* innerIndices = pAdjacencyMatrix.innerIndexPtr();
        return Edge{pPosition.outerIndex, static_cast<size_t>(innerIndices[pPosition.innerIndex])};
      }

      Iterator& operator++() {
        const int* outerIndices = pAdjacencyMatrix.outerIndexPtr();
        ++pPosition.innerIndex;
        while (pPosition.innerIndex >= static_cast<size_t>(outerIndices[pPosition.outerIndex + 1])) {
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
    Edges(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix) : pAdjacencyMatrix(adjacencyMatrix) {}

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

  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override { AdjMatGraph::addEdge(out, in, edgeWeight); }

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

}  // namespace graph
