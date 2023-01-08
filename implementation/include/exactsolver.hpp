#pragma once

#include <vector>

#include <draw/definitions.hpp>

#include "graph/graph.hpp"

/*!
 * \brief solveExact solve the tsp to optimality using Miller-Tucker-Zemlin formulation
 * \param euclidean
 * \details x_ij belongs to the column j* numberOfNodes +i
 * hereby we ommit th x_jj entries and store u_j in that place
 */
std::vector<unsigned int> solve(const Euclidean& euclidean, const ProblemType problemType);

class Index {
public:
  Index(const size_t numberOfNodes, const ProblemType type) : pNumberOfNodes(numberOfNodes), pType(type) {}

  size_t xVariables() const { return pNumberOfNodes * (pNumberOfNodes - 1); }
  size_t uVariables() const { return pNumberOfNodes - 1; }
  size_t cVariables() const { return 1; }
  size_t numVariables() const {return xVariables() + uVariables() + cVariables();}

  size_t variableX(const size_t i, const size_t j) const { return j * pNumberOfNodes + i; }
  size_t variableU(const size_t i) const { return i * pNumberOfNodes + i; }
  size_t variableC() const { return 0; }

  size_t xInConstraints() const { return pNumberOfNodes; }
  size_t xOutConstraints() const { return pNumberOfNodes; }
  size_t xConstraints() const { return xInConstraints() + xOutConstraints(); }
  size_t uConstraints() const { return (pNumberOfNodes - 1) * (pNumberOfNodes - 2); }
  size_t cConstraints() const { return (pType == ProblemType::BTSP_exact ? xVariables() : 0);}
  size_t numConstraints() const {return xConstraints() + uConstraints() + cConstraints();}

  size_t constraintXin(const size_t j) const { return j; }
  size_t constraintXout(const size_t j) const { return pNumberOfNodes + j; }
  size_t constraintU(const size_t i, const size_t j) const {
    return xConstraints() + (j - 1) * (pNumberOfNodes - 2) + (i > j ? i - 2 : i - 1);
  }
  size_t constraintC(const size_t i, const size_t j) const { return xConstraints() + uConstraints() + j * (pNumberOfNodes - 1) + (i>j ? i-1 : i);}

private:
  const size_t pNumberOfNodes;
  const ProblemType pType;
};
