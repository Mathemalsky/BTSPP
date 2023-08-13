/*
 * GRAPH is a library to store and manipulate graphs as adjacency list or
 * as sparse eigen matrix. Different specialized types of graphs are
 * supported.
 * Copyright (C) 2023 Jurek Rostalsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <Eigen/SparseCore>

#include "exceptions.hpp"
#include "geometry.hpp"
#include "utils.hpp"

namespace graph {

/***********************************************************************************************************************
 *                                            forward declare graph class
 **********************************************************************************************************************/

class Graph;

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
  requires(std::is_base_of_v<Graph, G>)
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
   * @details EdgeWeight(0) sets weight to infinity. This is because Eigen calls EdgeWeight(0) when
   * implicit entries of adjacency matrix are accessed. By the definition of + and * operator infinity
   * is the neutral element of + operation and absorbing element of * operation. To set an edge weight
   * to zero. EdgeWeight(0, true) muss be called.
   */
  explicit EdgeWeight(const double weight, const bool trueZero = false) {
    if (trueZero || weight != 0) {
      pWeight = weight;
    }
    else {
      pWeight = std::numeric_limits<double>::infinity();
    }
  }

  /*!
   * @brief returns the weight as double
   */
  double weight() const { return pWeight; }

  /*!
   * @brief compares for inquality
   */
  bool operator==(const EdgeWeight compare) { return pWeight == compare.pWeight; }

  /*!
   * @brief compares for inequality
   */
  bool operator!=(const EdgeWeight compare) { return pWeight != compare.pWeight; }

  /*!
   * @brief compares for \leq
   * @details This function is needed in Eigen.
   * @attention this function is used for comparison to prune elements in sparse matrix an more.
   * Not all are well defined for this + and * operator
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
 * @brief The Graph class is the base class of all more specialist types of graph. It's designed to contain the
 * functions graphs off all types have in common.
 * @details Currently no subclass implements adding nodes after calling the constructor.
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

  /*!
   * @brief fastWeight is a cheap function monotone in the weight of the edge
   * @details fastWeight is often cheaper than weight and it's monotone (for fixed subclass of graph) in the weight. So ist can be used for
   * faster sorting a graphs edges by weigth.
   * If no fastweight function is provided fall back to weight function.
   * @return fast weight
   */
  virtual double fastWeight(const Edge& e) const { return weight(e); };
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
   * @details
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
       * @return this iterator after increasing index
       */
      Iterator& operator++() {
        ++pPosition.index;
        makeValid();
        return *this;
      }

      /*!
       * @brief compares for inequality
       */
      bool operator!=(const Iterator& other) const { return pPosition.index != other.pPosition.index; }

      /*!
       * @brief increases position until edge is directed from higher to lower index
       */
      void makeValid() {
        const size_t outerIndex = pPosition.index / pPosition.numberOfNodes;
        if (pPosition.index >= outerIndex * (pPosition.numberOfNodes + 1)) {
          pPosition.index = (outerIndex + 1) * pPosition.numberOfNodes;
        }
      }

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
    Iterator begin() const {
      Iterator it(Iterator::Position{0, pNumberOfNodes});
      it.makeValid();
      return it;
    }

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
  Euclidean(std::vector<Point2D>&& positions) : pPositions(std::move(positions)) {}

  /*!
   * @brief number of Edges in graph. Depends only on the number of nodes since it is a complete graph
   */
  size_t numberOfEdges() const override { return pPositions.size() * (pPositions.size() - 1) / 2; }

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
   * @return square of distance
   */
  double weight(const Edge& e) const override { return dist(pPositions[e.u], pPositions[e.v]); }

  /*!
   * @brief returns square of euclidean distance between nodes
   * @details this is cheaper, because no root needs to be compted
   * @param e edge for which fast weight is requested
   * @return square of distance
   */
  double fastWeight(const Edge& e) const override { return distSquared(pPositions[e.u], pPositions[e.v]); }

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
   * @param in end vertex of edge
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
       * @return this iterator after incrementation
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
   * @brief appends empty entries to adjacency list
   * @param numberOfNodesToAdd number of isolated nodes that are added to the graph
   */
  void addIsolatedNodes(const size_t numberOfNodesToAdd) { pAdjacencyList.resize(numberOfNodes() + numberOfNodesToAdd); }

  /*!
   * @brief gives complete adjacency list
   * @details The i^th entry of the outer vectors contains the indices of all nodes reached with edges outgoing from node i.
   * @return read only reference to adjacency list
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
  template <typename Function>
  size_t neighbourAnyExcept(const size_t u, Function&& criteria) {
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
  template <typename Function>
  size_t neighbourAnyPrefer(const size_t u, Function&& criteria) {
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
       * @return this iterator after incrementation
       */
      Iterator& operator++() {
        assert(pPosition.outerIndex < pAdjacencyList.size() && "Iterator is already behind end!");

        ++pPosition.innerIndex;
        makeValid();
        return *this;
      }

      /*!
       * @brief compares for inequality
       */
      bool operator!=(const Iterator& other) const {
        return pPosition.outerIndex != other.pPosition.outerIndex || pPosition.innerIndex != other.pPosition.innerIndex;
      }

      /*!
       * @brief if invalid, set iterator to next valid position or to end iterator
       */
      void makeValid() {
        while (pPosition.outerIndex < pAdjacencyList.size() && (outOfNeighbours() || !toLowerIndex())) {
          if (outOfNeighbours()) {
            pPosition.innerIndex = 0;
            ++pPosition.outerIndex;
          }
          else {
            ++pPosition.innerIndex;
          }
        }
      }

    private:
      /*!
       * @brief checks if the inner index is still in the range of neighbours vector
       */
      bool outOfNeighbours() const { return pPosition.innerIndex >= pAdjacencyList[pPosition.outerIndex].size(); }

      /*!
       * @brief checks if inner index < outer index
       */
      bool toLowerIndex() const { return pAdjacencyList[pPosition.outerIndex][pPosition.innerIndex] < pPosition.outerIndex; }

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
      it.makeValid();
      return it;
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

  /*!
   * @brief constructor, resizes the internal adjacency list to numberOfNodes
   * @param numberOfNodes number of nodes in constructed graph
   */
  AdjacencyListGraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }

  /*!
   * @brief constructor, resizes the internal adjacency list, and reserves space for neighbours
   * @param numberOfNodes number of nodes in constructed graph
   * @param reserveNeighbours number of neighbours that can be added without reallocating
   */
  AdjacencyListGraph(const size_t numberOfNodes, const size_t reserveNeighbours) {
    pAdjacencyList.resize(numberOfNodes);
    for (std::vector<size_t>& vec : pAdjacencyList) {
      vec.reserve(reserveNeighbours);
    }
  }

  /*!
   * @brief constructor, sets the internal adjacency list to adjacencyList
   * @param adjacencyList adjacency list
   */
  AdjacencyListGraph(const std::vector<std::vector<size_t>>& adjacencyList) : AdjListGraph(adjacencyList) {}

  /*!
   * @brief adds an edge to the graph
   * @details adds also an reverse directed edge
   * @param out node at the edge's end
   * @param in node at the edge's end
   * @param edgeWeight not used
   */
  void addEdge(const size_t out, const size_t in, [[maybe_unused]] const EdgeWeight edgeWeight = EdgeWeight(1.0)) override {
    pAdjacencyList[out].push_back(in);
    pAdjacencyList[in].push_back(out);
  }

  /*!
   * @brief adds an edge to the graph
   * @details adds also an reverse directed edge
   * @param e edge to add
   * @param edgeWeight
   */
  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = EdgeWeight(1.0)) override {
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

  /*!
   * @brief creates iteratable instance of Edges
   * @return instance of Edges with adjacency list of this graph
   */
  Edges edges() const { return Edges(pAdjacencyList); }

  /*!
   * @brief removes an edge from the graph
   * @details removes copy for revrse directed edge as well, removes only one per direction
   * @param e edge to be removed
   */
  void removeEdge(const Edge& e) {
    [[maybe_unused]] const bool removed = removeAnyElementByValue(pAdjacencyList[e.u], e.v);
    assert(removed && "Edge to be removed does not exist in graph!");
    [[maybe_unused]] const bool removed2 = removeAnyElementByValue(pAdjacencyList[e.v], e.u);
    assert(removed2 && "Edge to be removed does not exist in graph!");
  }

  /*!
   * @brief removes an edge from the graph
   * @details removes copy for revrse directed edge as well, removes only one per direction
   * @param u one end of the edge
   * @param v other end of the edge
   */
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
  /*!
   * @brief Edges is a facade class to iterate over all edges in AdjacencyListDigraph graph.
   */
  class Edges {
  private:
    /*!
     * @brief Iterator on edges of an AdjacencyListDigraph graph
     */
    class Iterator {
    public:
      /*!
       * @brief Position in adjacency list containing outer index and inner index
       */
      struct AdjListPos {
        size_t outerIndex;
        size_t innerIndex;
      };

      /*!
       * @brief creates an Iterator poining to pos in adjacencyList
       * @param adjacencyList adjacency list of graph
       * @param pos position
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
       * @brief moves iterator forward by one position
       * @details if at the end of neighbours of a node set Iterator to next nodes neighbours
       * @return this iterator after incrementation
       */
      Iterator& operator++() {
        ++pPosition.innerIndex;
        makeValid();
        return *this;
      }

      /*!
       * @brief compares iterators for inequality
       * @param other iterator to compare with
       * @return if the iterators are different
       */
      bool operator!=(const Iterator& other) const {
        return pPosition.outerIndex != other.pPosition.outerIndex || pPosition.innerIndex != other.pPosition.innerIndex;
      }

      /*!
       * @brief if the iterator is not pointing to valid element, set the iterator to next valid element
       */
      void makeValid() {
        while (pPosition.outerIndex < pAdjacencyList.size() && pPosition.innerIndex >= pAdjacencyList[pPosition.outerIndex].size()) {
          pPosition.innerIndex = 0;
          ++pPosition.outerIndex;
        }
      }

    private:
      const std::vector<std::vector<size_t>>& pAdjacencyList; /**< i-th vector contains neighbours of i */
      AdjListPos pPosition;                                   /**< outer and inner index of current position */
    };                                                        // end Iterator class

  public:
    /*!
     * @brief creates edges instance of edges
     * @param adjacencyList adjacency list of graph
     */
    Edges(const std::vector<std::vector<size_t>>& adjacencyList) : pAdjacencyList(adjacencyList) {}

    /*!
     * @brief creates an iterator pointing in front of the first edge
     * @return begin iterator
     */
    Iterator begin() const {
      Iterator it(pAdjacencyList, Iterator::AdjListPos{0, 0});
      it.makeValid();
      return it;
    }

    /*!
     * \brief end returns the end Iterator for comparison.
     * \details The 0 is just an arbitrary number because the outerIndex of AdjListPos is not taken into
     * account for comparison.
     * \return Iterator behind the last element.
     */
    Iterator end() const { return Iterator(pAdjacencyList, Iterator::AdjListPos{pAdjacencyList.size(), 0}); }

  private:
    const std::vector<std::vector<size_t>>& pAdjacencyList; /**< i-th vector contains neighbours of i */
  };                                                        // end Edges class

public:
  AdjacencyListDigraph()  = default;
  ~AdjacencyListDigraph() = default;

  AdjacencyListDigraph(const AdjacencyListDigraph& graph) = default;

  /*!
   * @brief constructor, resizes the internal adjacency list to numberOfNodes
   * @param numberOfNodes number of nodes in constructed graph
   */
  AdjacencyListDigraph(const size_t numberOfNodes) { pAdjacencyList.resize(numberOfNodes); }

  /*!
   * @brief constructor, sets the internal adjacency list to adjacencyList
   * @param adjacencyList adjacency list
   */
  AdjacencyListDigraph(const std::vector<std::vector<size_t>>& adjacencyList) : AdjListGraph(adjacencyList) {}

  /*!
   * @brief adds an edge to the graph
   * @param out node at the edge's begin
   * @param in node at the edge is directed to
   * @param edgeWeight not used
   */
  void addEdge(const size_t out, const size_t in, [[maybe_unused]] const EdgeWeight edgeWeight = EdgeWeight(1.0)) override {
    pAdjacencyList[out].push_back(in);
  }

  /*!
   * @brief adds an edge to the graph
   * @param e edge to be added
   * @param edgeWeight
   */
  void addEdge(const Edge& e, [[maybe_unused]] const EdgeWeight edgeWeight = EdgeWeight(1.0)) override {
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

  /*!
   * @brief returns number of edges in the graph
   * @details counts the entries in adjacency list
   * @return number of edges in the graph
   */
  size_t numberOfEdges() const override {
    return std::accumulate(pAdjacencyList.begin(), pAdjacencyList.end(), 0, [](const unsigned int sum, const std::vector<size_t>& vec) {
      return sum + vec.size();
    });
  }

  /*!
   * @brief checks if the graph is biconnected
   * @return true if the graph is biconnected
   */
  bool biconnected() const { return checkBiconnectivity(this->undirected()); }

  /*!
   * @brief creates iteratable instance of Edges
   * @return instance of Edges with adjacency list of this graph
   */
  Edges edges() const { return Edges(pAdjacencyList); }

  /*!
   * @brief removes an edge from the graph
   * @param e edge to be removed
   */
  void removeEdge(const Edge& e) {
    [[maybe_unused]] const bool removed = removeAnyElementByValue(pAdjacencyList[e.u], e.v);
    assert(removed && "Edge to be removed does not exist in graph!");
  }

  /*!
   * @brief removes an edge from the graph
   * @param u outgoing node of the edge to be removed
   * @param v ingoing node of the edge to be removed
   */
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
  /*!
   * @brief Neighbours is a facade class to iterate over all neighbours of an node
   */
  class Neighbours {
  private:
    /*!
     * @brief Iterator on neighbours of a node of the graph
     */
    class Iterator {
    public:
      /*!
       * @brief creates an iterator from given pointer
       * @param ptr pointer inner index vector
       */
      Iterator(const int* ptr) : pInnerIndices(ptr) {}

      /*!
       * @brief dereferences the pointer to the current inner index
       * @return convert int to size_t
       */
      size_t operator*() const { return static_cast<size_t>(*pInnerIndices); }

      /*!
       * @brief increments the iterator
       * @return the iterator after incrementation
       */
      Iterator operator++() {
        ++pInnerIndices;
        return *this;
      }

      /*!
       * @brief compares iterators for inequality
       * @param other Iterator to compare with
       * @return true if iterators are different
       */
      bool operator!=(const Iterator& other) const { return pInnerIndices != other.pInnerIndices; }

    private:
      const int* pInnerIndices; /**< pointer to inner index */
    };                          // end Iterator class

  public:
    /*!
     * @brief
     * @param adjacencyMatrix
     * @param node
     */
    Neighbours(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const size_t node) :
      pAdjacencyMatrix(adjacencyMatrix),
      pNode(node) {
      assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
    }

    /*!
     * @brief creates begin iterator
     * @return Iterator pointing to the first inner index
     */
    Iterator begin() const { return Iterator(pAdjacencyMatrix.innerIndexPtr() + pAdjacencyMatrix.outerIndexPtr()[pNode]); }

    /*!
     * @brief creates end iterator
     * @return Iterator pointing behind the last inner index
     */
    Iterator end() const { return Iterator(pAdjacencyMatrix.innerIndexPtr() + pAdjacencyMatrix.outerIndexPtr()[pNode + 1]); }

  private:
    const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix; /**< the sparse adjacency matrix */
    const size_t pNode;                                                       /**< node index, iterate over neighbours of this node */
  };                                                                          // end Neighbours class

public:
  AdjMatGraph()  = default;
  ~AdjMatGraph() = default;

  /*!
   * @brief constructs adjacency matrix graph with quadratic adjacency matrix of given size
   * @param numberOfNodes dimension of adjacency matrix
   */
  AdjMatGraph(const size_t numberOfNodes) : pAdjacencyMatrix(numberOfNodes, numberOfNodes) {}

  /*!
   * @brief constructs adjacency matrix graph from sparse eigen matrix
   * @param mat
   */
  AdjMatGraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : pAdjacencyMatrix(mat) {
    assert(mat.rows() == mat.cols() && "Adjacency matrix should be quadratic!");
  }

  /*!
   * @brief constructs sparse matrix from triplet list and adjacency matrix graph from that matrix
   * @param numberOfNodes number of nodes in the graph
   * @param tripletList list of triplets (row index, column index, entry)
   */
  AdjMatGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    pAdjacencyMatrix(Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>(numberOfNodes, numberOfNodes)) {
    pAdjacencyMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  /*!
   * @brief checks if there is an explicit entry for that edge
   * @param u node at start of edge
   * @param v node at end of edge
   * @return true if edge exists
   */
  virtual bool adjacent(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v) != EdgeWeight(0.0); }

  /*!
   * @brief checks if there is an explicit entry for that edge
   * @param e edge to check
   * @return true if edge exists
   */
  virtual bool adjacent(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v) != EdgeWeight(0.0); }

  /*!
   * @brief accesses the weight of an edge
   * @param u node at start of edge
   * @param v node at end of edge
   * @return weight of edge
   */
  double weight(const size_t u, const size_t v) const override { return pAdjacencyMatrix.coeff(u, v).weight(); }

  /*!
   * @brief accesses the weight of an edge
   * @param e edge whose weight is to check
   * @return weight of edge
   */
  double weight(const Edge& e) const override { return pAdjacencyMatrix.coeff(e.u, e.v).weight(); }

  /*!
   * @brief number of edges is number of nonzeros
   * @return number of edges
   */
  size_t numberOfEdges() const override { return pAdjacencyMatrix.nonZeros(); }

  /*!
   * @brief number of nodes is number of columns in matrix which equals number of rows in matrix
   * @return number of columns in matrix
   */
  size_t numberOfNodes() const override { return pAdjacencyMatrix.cols(); }

  /*!
   * @brief calls makeCompressed() on the underlying sparse eigen matrix
   */
  void compressMatrix() { pAdjacencyMatrix.makeCompressed(); }

  /*!
   * @brief provides read only accessto the underlying sparse matrix
   * @return const reference to sparse eigen matrix
   */
  const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& matrix() const { return pAdjacencyMatrix; }

  /*!
   * @brief creates neighbours facade object to iterate over all neighbours of this graph
   * @param node whose neighbours are to iterate over
   * @return neighbours facade
   */
  Neighbours neighbours(const size_t node) const { return Neighbours(pAdjacencyMatrix, node); }

  /*!
   * @brief makes explicit zeros implicit zeros
   * @attention eigen and EdgeWeight are not yet sufficiently adjusted to each other to use this function
   */
  // void pruneMatrix() { pAdjacencyMatrix.prune(0.0, Eigen::NumTraits<EdgeWeight>::dummy_precision()); }

  /*!
   * @brief computes the square of a graph
   * @attention eigen and EdgeWeight are not yet sufficiently adjusted to each other to use this function
   */
  // void square() { pAdjacencyMatrix = pAdjacencyMatrix * pAdjacencyMatrix; }  // operator*= not provided by Eigen

protected:
  Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor> pAdjacencyMatrix; /**< sparse adjacency matrix */

  /*!
   * @brief adds an egde to the graph by inserting an entry in adjacency matrix
   * @param out outgoing node
   * @param in ingoing node
   * @param edgeWeight weight of the edge
   */
  virtual void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    pAdjacencyMatrix.insert(out, in) = edgeWeight;
  }

  /*!
   * @brief adds an egde to the graph by inserting an entry in adjacency matrix
   * @param e edge to insert
   * @param edgeWeight weight of the edge
   */
  virtual void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { pAdjacencyMatrix.insert(e.u, e.v) = edgeWeight; }
};

/*!
 * \brief The AdjacencyMatrixGraph class implements the concrete classes for undirected graphs.
 * \details For the sake of faster iteration over all neighbours, each edge is stored twice: once for each possible
 * direction.
 */
class AdjacencyMatrixGraph : public AdjMatGraph, public UndirectedGraph {
private:
  /*!
   * @brief Edges is a facade class to iterate over all edges of a AdjacencyMatrixGraph
   */
  class Edges {
  private:
    /*!
     * @brief Iterator to iterator of a AdjacencyMatrixGraphs edges
     */
    class Iterator {
    public:
      /*!
       * @brief bundle inner and outer index
       */
      struct SparseMatrixPos {
        size_t innerIndex; /**< inner index is column index */
        size_t outerIndex; /**< outer index is row index */
      };

      /*!
       * @brief creates an Iterator on the AdjacencyMatrixGraph
       * @param adjacencyMatrix of the AdjacencyMatrixGraph
       * @param pos position in sparse matrix
       */
      Iterator(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const SparseMatrixPos& pos) :
        pAdjacencyMatrix(adjacencyMatrix),
        pPosition(pos) {
        assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
      }

      /*!
       * @brief creates Edge for current edge
       * @return Edge
       */
      Edge operator*() const {
        const int* const innerIndices = pAdjacencyMatrix.innerIndexPtr();
        return Edge{pPosition.outerIndex, static_cast<size_t>(innerIndices[pPosition.innerIndex])};
      }

      /*!
       * @brief increments iterator
       * @return Iterator after incrementation
       */
      Iterator& operator++() {
        const int* const outerIndices = pAdjacencyMatrix.outerIndexPtr();
        const int* const innerIndices = pAdjacencyMatrix.innerIndexPtr();
        ++pPosition.innerIndex;
        while (pPosition.innerIndex < static_cast<size_t>(pAdjacencyMatrix.nonZeros()) && !valid()) {
          if (static_cast<size_t>(innerIndices[pPosition.innerIndex]) >= pPosition.outerIndex) {
            pPosition.innerIndex = outerIndices[pPosition.outerIndex + 1];  // skip rest of the row
          }
          ++pPosition.outerIndex;                                           // goes to next row
        }
        return *this;
      }

      /*!
       * @brief compares iterators for inequality
       * @param other iterator t compare with
       * @return true if iterators are diffrent
       */
      bool operator!=(const Iterator& other) const { return pPosition.innerIndex != other.pPosition.innerIndex; }

      /*!
       * @brief checks if the iterator points to a valid element and this element is in strictly lower left triangle of matrix
       * @return true if iterator points to fitting entry
       */
      bool valid() const {
        const int* const outerIndices = pAdjacencyMatrix.outerIndexPtr();
        const int* const innerIndices = pAdjacencyMatrix.innerIndexPtr();
        return pPosition.innerIndex < static_cast<size_t>(outerIndices[pPosition.outerIndex + 1])  // correct row
               && static_cast<size_t>(innerIndices[pPosition.innerIndex]) < pPosition.outerIndex;  // lower triangular
      }

    private:
      const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix; /**< seigen sparse matrix as adjacency matrix of graph */
      SparseMatrixPos pPosition;                                                /**< position of iterator in inner and outer indices */
    };                                                                          // end Iterator class

  public:
    /*!
     * @brief constructs facade class from adjacency matrix of this graph
     * @param adjacencyMatrix
     */
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
    const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix; /**< adjacency matrix from the asociated graph */
  };                                                                          // end Edges class

public:
  AdjacencyMatrixGraph()  = default;
  ~AdjacencyMatrixGraph() = default;

  /*!
   * @brief constructs, AdjacencyMatrixGraph from adjacency matrix
   * @param mat adjacency matrix
   */
  AdjacencyMatrixGraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : AdjMatGraph(mat) {}

  /*!
   * @brief constructs AdjacencyMatrixGraph from triplet list
   * @param numberOfNodes number of nodes, which is the dimension of the quadratic adjacency matrix
   * @param tripletList list of triplets (row index, column index, edgeweight)
   */
  AdjacencyMatrixGraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {}

  /*!
   * @brief adds an edge as entry in the adjacency matrix
   * @attention In wortst case this takes O(numberOfNonZeros), when ever possible add all needed edges before creating the adjacency matrix
   * @param out outgoing node
   * @param in ingoing node
   * @param edgeWeight weight of edge
   */
  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override {
    AdjMatGraph::addEdge(out, in, edgeWeight);
    AdjMatGraph::addEdge(in, out, edgeWeight);
  }

  /*!
   * @brief adds an edge as entry in the adjacency matrix
   * @attention In wortst case this takes O(numberOfNonZeros), when ever possible add all needed edges before creating the adjacency matrix
   * @param e edge to be added
   * @param edgeWeight weight of edge
   */
  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { addEdge(e.u, e.v, edgeWeight); }

  /*!
   * @brief checks if the graph is connected
   * @details this is accomplished by a dfs
   * @return true if graph is connected
   */
  bool connected() const override;

  /*!
   * @brief checks if the graph is biconnected
   * @details this is accomplished by schmidts algorithm
   * @return true if graph is biconnected
   */
  bool biconnected() const { return checkBiconnectivity(*this); }

  /*!
   * @brief creates a facade to iterate over the edges in the graph
   * @return facade class Edges
   */
  Edges edges() const { return Edges(pAdjacencyMatrix); }

  /*!
   * @brief removes an edge from the graph (and the reverse directed copy too)
   * @param e edge to be removed
   * @attention sets weight to zero, does not make the entry implicit
   * @attention eigen and EdgeWeight are not yet sufficiently adjusted to each other to use this function
   */
  /*
  void removeEdge(const Edge& e) {
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != EdgeWeight(0.0) && "Edge to be removed does not exist in graph!");
    pAdjacencyMatrix.coeffRef(e.u, e.v) = EdgeWeight(0.0);
    assert(pAdjacencyMatrix.coeff(e.v, e.u) != EdgeWeight(0.0) && "Edge to be removed does not exist in graph!");
    pAdjacencyMatrix.coeffRef(e.v, e.u) = EdgeWeight(0.0);
  }*/
};

/*!
 * \brief The AdjacencyMatrixGraph class implements the concrete class for directed graphs.
 * \details The functions checking for connectivity are checking for connectivity in the sense of weak connectivity.
 * AdjacencyMatrixDigraph may have antiparallel edges.
 */
class AdjacencyMatrixDigraph : public AdjMatGraph, public DirectedGraph {
private:
  /*!
   * @brief Edges is a facade class to iterate over all edges of a AdjacencyMatrixGraph
   */
  class Edges {
  private:
    /*!
     * @brief Iterator to iterator of a AdjacencyMatrixGraphs edges
     */
    class Iterator {
    public:
      /*!
       * @brief bundle inner and outer index
       */
      struct SparseMatrixPos {
        size_t innerIndex; /**< inner index is column index */
        size_t outerIndex; /**< outer index is row index */
      };

      /*!
       * @brief creates an Iterator on the AdjacencyMatrixGraph
       * @param adjacencyMatrix of the AdjacencyMatrixGraph
       * @param pos position in sparse matrix
       */
      Iterator(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix, const SparseMatrixPos& pos) :
        pAdjacencyMatrix(adjacencyMatrix),
        pPosition(pos) {
        assert(pAdjacencyMatrix.isCompressed() && "Iterating over uncompressed matrix results in undefined behavior!");
      }

      /*!
       * @brief creates Edge for current edge
       * @return Edge
       */
      Edge operator*() const {
        const int* innerIndices = pAdjacencyMatrix.innerIndexPtr();
        return Edge{pPosition.outerIndex, static_cast<size_t>(innerIndices[pPosition.innerIndex])};
      }

      /*!
       * @brief increments iterator
       * @return Iterator after incrementation
       */
      Iterator& operator++() {
        const int* outerIndices = pAdjacencyMatrix.outerIndexPtr();
        ++pPosition.innerIndex;
        while (pPosition.innerIndex >= static_cast<size_t>(outerIndices[pPosition.outerIndex + 1])) {
          ++pPosition.outerIndex;
        }
        return *this;
      }

      /*!
       * @brief compares iterators for inequality
       * @param other iterator t compare with
       * @return true if iterators are diffrent
       */
      bool operator!=(const Iterator& other) const { return pPosition.innerIndex != other.pPosition.innerIndex; }

    private:
      const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix; /**< seigen sparse matrix as adjacency matrix of graph */
      SparseMatrixPos pPosition;                                                /**< position of iterator in inner and outer indices */
    };                                                                          // end Iterator class

  public:
    /*!
     * @brief constructs facade class from adjacency matrix of this graph
     * @param adjacencyMatrix
     */
    Edges(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& adjacencyMatrix) : pAdjacencyMatrix(adjacencyMatrix) {}

    /*!
     * \brief begin returns a an Iterator to the begin of the strictly lower triangle in adjacency matrix.
     * \details The Iterator starts in second row, because the intersect of first row and strictly lower triangle is
     * empty.
     * \return Iterator to the begin of the strictly lower triangle in adjacency matrix
     */
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
    const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& pAdjacencyMatrix; /**< adjacency matrix from the asociated graph */
  };                                                                          // end Edges class

public:
  AdjacencyMatrixDigraph()  = default;
  ~AdjacencyMatrixDigraph() = default;

  /*!
   * @brief constructs digraph from sparse eigen row major matrix
   * @param mat sparse eigen matrix as adjacency matrix
   */
  AdjacencyMatrixDigraph(const Eigen::SparseMatrix<EdgeWeight, Eigen::RowMajor>& mat) : AdjMatGraph(mat) {}

  /*!
   * @brief constructs digraph from triplet list
   * @param numberOfNodes
   * @param tripletList list of triplets (row index, column index, edgeweight)
   */
  AdjacencyMatrixDigraph(const size_t numberOfNodes, const std::vector<Eigen::Triplet<EdgeWeight>>& tripletList) :
    AdjMatGraph(numberOfNodes, tripletList) {}

  /*!
   * @brief adds an edge to the graph
   * @attention In wortst case this takes O(numberOfNonZeros), when ever possible add all needed edges before creating the adjacency matrix
   * @param out outgoing node
   * @param in ingoing node
   * @param edgeWeight weight of edge
   */
  void addEdge(const size_t out, const size_t in, const EdgeWeight edgeWeight) override { AdjMatGraph::addEdge(out, in, edgeWeight); }

  /*!
   * @brief adds an edge to the graph
   * @attention In wortst case this takes O(numberOfNonZeros), when ever possible add all needed edges before creating the adjacency matrix
   * @param e edge to be added
   * @param edgeWeight weight of edge
   */
  void addEdge(const Edge& e, const EdgeWeight edgeWeight) override { addEdge(e.u, e.v, edgeWeight); }

  /*!
   * @brief checks if the graph is weakly connected
   * @details weakly connected means that the undirected graph obtained by ignoring the edges direction is connected
   * @return true if graph is weakly connected
   */
  bool connected() const override { return undirected().connected(); };

  /*!
   * @brief checks if the graph is weakly  biconnected
   * @details weakly biconnected means that the undirected graph obtained by ignoring the edges direction is biconnected
   * @return
   */
  bool biconnected() const { return checkBiconnectivity(undirected()); }

  /*!
   * @brief returns instance of edges facade class to iterator over all edges
   * @return edges facade
   */
  Edges edges() const { return Edges(pAdjacencyMatrix); }

  /*!
   * @brief removes an edge from the graph (and the reverse directed copy too)
   * @param e edge to be removed
   * @attention sets weight to zero, does not make the entry implicit
   * @attention eigen and EdgeWeight are not yet sufficiently adjusted to each other to use this function
   */
  /* void removeEdge(const Edge& e) {
    assert(pAdjacencyMatrix.coeff(e.u, e.v) != EdgeWeight(0.0) && "Edge to be removed does not exist in graph!");
    pAdjacencyMatrix.coeffRef(e.u, e.v) = EdgeWeight(0.0);
  } */

  AdjacencyMatrixGraph undirected() const;
};

/*!
 * \brief The Tree class is an abstract which implements the methods all trees have in common.
 */
class Tree : public virtual Graph {
  /*!
   * @brief trees are always connected and hence the return is always true
   * @return true
   */
  bool connected() const override { return true; }

  /*!
   * @brief a tree has one edge less than nodes
   * @return number of nodes -1
   */
  size_t numberOfEdges() const override { return numberOfNodes() - 1; }
};

/*!
 * \brief class DfsTree
 * \details This class is for trees where every node except for the root has exactly one parent and an edge directed to
 * it's parent. This allows to store the neighboors more efficient.
 */
class DfsTree : public Tree, public DirectedGraph {
private:
  /*!
   * @brief Edges is a facade class to iterate over all edges in the DfsTree
   */
  class Edges {
  private:
    /*!
     * @brief iterator on a DfsTrees edges
     */
    class Iterator {
    public:
      /*!
       * @brief creates iterator object
       * @param adjacencyList the trees adjacency list
       * @param root index of root node
       * @param pos iterator position
       */
      Iterator(const std::vector<size_t>& adjacencyList, const size_t root, const size_t pos) :
        pAdjacencyList(adjacencyList),
        pRoot(root),
        pPosition(pos) {}

      /*!
       * @brief creates an edge object for the edge of current iteration
       * @return Edge starting at current position
       */
      Edge operator*() const {
        assert(pPosition != pRoot && "No outging edge from root node!");
        return Edge{pPosition, pAdjacencyList[pPosition]};
      }

      /*!
       * @brief increments the iterator, the root is skipped
       * @return this iterator after incrementation
       */
      Iterator& operator++() {
        ++pPosition;
        if (pPosition == pRoot) {
          ++pPosition;
        }
        return *this;
      }

      /*!
       * @brief compares for inequality
       * @param other Iterator to compare
       * @return true if iterators are different
       */
      bool operator!=(const Iterator& other) const { return pPosition != other.pPosition; }

    private:
      const std::vector<size_t>& pAdjacencyList; /**< the tree's adjacency list */
      const size_t pRoot;                        /**< root index */
      size_t pPosition;                          /**< position in adjacenecy list */
    };                                           // end Iterator class

  public:
    /*!
     * @brief constructs facade class from adjacency list of this graph
     * @param adjacencyList list of parent nodes
     * @param root of tree
     */
    Edges(const std::vector<size_t>& adjacencyList, const size_t root) : pAdjacencyList(adjacencyList), pRoot(root) {}

    /*!
     * @brief creates iterator to first node
     * @details if first node is root, set iterator to second
     * @return begin iterator
     */
    Iterator begin() const {
      if (pRoot != 0) {
        return Iterator(pAdjacencyList, pRoot, 0);  // skip the root node
      }
      else {
        return Iterator(pAdjacencyList, pRoot, 1);
      }
    }

    /*!
     * @brief creates iterator behind last element
     * @return end iterator
     */
    Iterator end() const { return Iterator(pAdjacencyList, pRoot, pAdjacencyList.size()); }

  private:
    const std::vector<size_t>& pAdjacencyList; /**< adjacency list */
    const size_t pRoot;                        /**< root node index */
  };                                           // end Edges class

public:
  DfsTree()  = default;
  ~DfsTree() = default;

  /*!
   * @brief constructs DfsTree from number of nodes
   * @details resizes the adjacency list to numberOfNodes, reserve numberOfNodes entries in exploration Order
   * @param numberOfNodes number of nodes in dfsTree
   */
  DfsTree(const size_t numberOfNodes) : pAdjacencyList(numberOfNodes) {
    pExplorationOrder.reserve(numberOfNodes);  // just reserve, because dfs performs push_backs
  }

  /*!
   * @brief checks if node v is the parent of u
   * @param u
   * @param v
   * @return true if v is parent of u
   */
  bool adjacent(const size_t u, const size_t v) const override { return u != 0 && v == parent(u); }

  /*!
   * @brief checks if the edge is in the graph
   * @param e edge
   * @return true if e is in the graph
   */
  bool adjacent(const Edge& e) const override { return e.u != 0 && e.v == parent(e.u); }

  /*!
   * @brief returns number of nodes
   * @return size of pAdjacencyList
   */
  size_t numberOfNodes() const override { return pAdjacencyList.size(); }

  /*!
   * @brief creates iteratable instance of Edges
   * @return instance of Edges with adjacency list and root of this graph
   */
  Edges edges() const { return Edges(pAdjacencyList, root()); }

  /*!
   * @brief read only access to exploration order
   * @return const reference to exploration order
   */
  const std::vector<size_t>& explorationOrder() const { return pExplorationOrder; }

  /*!
   * @brief writable access to exploration order
   * @return reference to exploration order
   */
  std::vector<size_t>& explorationOrder() { return pExplorationOrder; }

  /*!
   * @brief returns the parent of u
   * @param u node
   * @return index of u's parent
   */
  size_t parent(const size_t u) const {
    assert(u != root() && "You are trying to access the root nodes parent, which is uninitialized memory!");
    return pAdjacencyList[u];
  }

  /*!
   * @brief writableaccess to the parent of u
   * @param u node
   * @return reference index of u's parent
   */
  size_t& parent(const size_t u) {
    assert(u != root() && "You are trying to access the root nodes parent, which is uninitialized memory!");
    return pAdjacencyList[u];
  }

  /*!
   * @brief returns the root node
   * @return the first explored node is the root node
   */
  size_t root() const { return pExplorationOrder.front(); }

private:
  std::vector<size_t> pAdjacencyList;    /**< adjancy list, all nodes except root have exatly one parent */
  std::vector<size_t> pExplorationOrder; /**< order of node exploration, starts with root */
};

}  // namespace graph
