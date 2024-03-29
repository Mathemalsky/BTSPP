/*
 * BTSPP is a tool to solve, approximate and draw instances of BTSVPP,
 * BTSPP, BTSP and TSP. Drawing is limited to euclidean graphs.
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
#include "solve/exactsolver.hpp"

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/SparseCore>

#include <Highs.h>

// graph library
#include "graph.hpp"

#include "exception/exceptions.hpp"

#include "solve/commonfunctions.hpp"

namespace exactsolver {

using Entry = Eigen::Triplet<double>;

static constexpr double M_INFINITY = 1e32;

/***********************************************************************************************************************
 *                                                  output function
 **********************************************************************************************************************/

void printInfo(const exactsolver::Result& res, const ProblemType problemType, const double runtime) {
  std::cout << "-------------------------------------------------------\n";
  std::cout << "Solved an instance of " << problemType << " using HiGHS library." << std::endl;
  std::cout << "OPT                                  : " << res.opt << std::endl;
  if (runtime != -1.0) {
    std::cout << "elapsed time                         : " << runtime << " ms\n";
  }
}

/***********************************************************************************************************************
 *                                            algorithms for BTSP & BTSPP
 **********************************************************************************************************************/

class Index {
public:
  Index(const size_t numberOfNodes, const ProblemType type) : pNumberOfNodes(numberOfNodes), pType(type) {}

  size_t xVariables() const { return pNumberOfNodes * (pNumberOfNodes - 1); }
  size_t uVariables() const { return pNumberOfNodes - 1; }
  size_t cVariables() const { return 1; }
  size_t numVariables() const { return xVariables() + uVariables() + cVariables(); }

  size_t variableX(const size_t i, const size_t j) const { return j * pNumberOfNodes + i; }
  size_t variableU(const size_t i) const { return i * pNumberOfNodes + i; }
  size_t variableC() const { return 0; }

  size_t xInConstraints() const { return pNumberOfNodes; }
  size_t xOutConstraints() const { return pNumberOfNodes; }
  size_t xConstraints() const { return xInConstraints() + xOutConstraints(); }
  size_t uConstraints() const { return (pNumberOfNodes - 1) * (pNumberOfNodes - 2); }
  size_t cConstraints() const { return (pType == ProblemType::BTSP_exact || pType == ProblemType::BTSPP_exact ? xVariables() : 0); }
  size_t numConstraints() const { return xConstraints() + uConstraints() + cConstraints(); }

  size_t constraintXin(const size_t j) const { return j; }
  size_t constraintXout(const size_t j) const { return pNumberOfNodes + j; }
  size_t constraintU(const size_t i, const size_t j) const {
    return xConstraints() + (j - 1) * (pNumberOfNodes - 2) + (i > j ? i - 2 : i - 1);
  }
  size_t constraintC(const size_t i, const size_t j) const {
    return xConstraints() + uConstraints() + j * (pNumberOfNodes - 1) + (i > j ? i - 1 : i);
  }

private:
  const size_t pNumberOfNodes;
  const ProblemType pType;
};

static void setTSPcost(HighsModel& model, const Index& index, const graph::Euclidean& euclidean, const size_t numberOfNodes) {
  model.lp_.col_cost_.resize(model.lp_.num_col_);
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      const double dist                          = euclidean.weight(i, j);
      model.lp_.col_cost_[index.variableX(i, j)] = dist;
      model.lp_.col_cost_[index.variableX(j, i)] = dist;  // exploiting symmetry
    }
    model.lp_.col_cost_[index.variableU(j)] = 0.0;        // here we store u_j instead of x_jj
  }
}

static void setBTSPcost(HighsModel& model, const Index& index) {
  model.lp_.col_cost_                    = std::vector<double>(model.lp_.num_col_, 0.0);  // set vector to 0
  model.lp_.col_cost_[index.variableC()] = 1.0;                                           // objective is c
}

static void setMillerTuckerZemlinBounds(HighsModel& model, const Index& index, const size_t numberOfNodes) {
  const double p = numberOfNodes;

  model.lp_.col_lower_ = std::vector<double>(model.lp_.num_col_, 0.0);  // set lower bound of variables to 0

  model.lp_.col_upper_.resize(model.lp_.num_col_);
  model.lp_.integrality_.resize(model.lp_.num_col_);

  // iterate over all variables
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      model.lp_.col_upper_[index.variableX(i, j)]   = 1.0;                     // set upper bound of varibales to 1
      model.lp_.col_upper_[index.variableX(j, i)]   = 1.0;                     // exploiting symmetry
      model.lp_.integrality_[index.variableX(i, j)] = HighsVarType::kInteger;  // constrain x_ij to be \in {0,1}
      model.lp_.integrality_[index.variableX(j, i)] = HighsVarType::kInteger;  // exploiting symmetry
    }
    model.lp_.col_upper_[index.variableU(j)]   = p - 2;                        // set upper bound of u_ij to p-2
    model.lp_.integrality_[index.variableU(j)] = HighsVarType::kContinuous;    // allow real numbers for u_ij
  }

  model.lp_.row_lower_.resize(model.lp_.num_row_);
  model.lp_.row_upper_.resize(model.lp_.num_row_);

  // iterate over all constraints
  // inequalities for fixing in and out degree to 1
  for (size_t i = 0; i < index.xConstraints(); ++i) {
    model.lp_.row_lower_[i] = 1.0;
    model.lp_.row_upper_[i] = 1.0;
  }
  // inequalities for guaranteeing connectednes
  for (size_t i = 0; i < index.uConstraints(); ++i) {
    model.lp_.row_lower_[index.xConstraints() + i] = -(p - 1);  // lower bound -infinity sharpened to -(p - 1)
    model.lp_.row_upper_[index.xConstraints() + i] = p - 1;     // upper bound p-1
  }
}

static void setCBounds(HighsModel& model, const Index& index) {
  model.lp_.col_upper_[index.variableC()] = M_INFINITY;
  model.lp_.col_lower_[index.variableC()] = 0.0;

  for (size_t i = index.xConstraints() + index.uConstraints(); i < index.numConstraints(); ++i) {
    model.lp_.row_lower_[i] = 0.0;
    model.lp_.row_upper_[i] = M_INFINITY;
  }
  model.lp_.integrality_[index.variableC()] = HighsVarType::kContinuous;
}

static void setPathBounds(HighsModel& model, const Index& index, const size_t s = 0, const size_t t = 1) {
  model.lp_.row_lower_[index.constraintXin(s)]  = 0;
  model.lp_.row_upper_[index.constraintXin(s)]  = 0;
  model.lp_.row_lower_[index.constraintXout(t)] = 0;
  model.lp_.row_upper_[index.constraintXout(t)] = 0;
}

static void setAntiCrossingBounds(HighsModel& model, const size_t numAntiCrossingConstraints) {
  const size_t numberOfPreviousConstraints = model.lp_.row_lower_.size();
  model.lp_.row_lower_.resize(numberOfPreviousConstraints + numAntiCrossingConstraints);
  model.lp_.row_upper_.resize(numberOfPreviousConstraints + numAntiCrossingConstraints);
  for (size_t i = numberOfPreviousConstraints; i < numberOfPreviousConstraints + numAntiCrossingConstraints; ++i) {
    model.lp_.row_lower_[i] = 0.0;
    model.lp_.row_upper_[i] = 1.0;
  }
}

static void setMillerTuckerZemlinMatrix(std::vector<Entry>& entries, const Index& index, const size_t numberOfNodes) {
  // inequalities for fixing in and out degree to 1
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = 0; i < j; ++i) {
      entries.push_back(Entry(index.constraintXin(j), index.variableX(i, j), 1.0));   // sum over in degree
      entries.push_back(Entry(index.constraintXout(j), index.variableX(j, i), 1.0));  // sum over out degree
    }
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      entries.push_back(Entry(index.constraintXin(j), index.variableX(i, j), 1.0));   // sum over in degree
      entries.push_back(Entry(index.constraintXout(j), index.variableX(j, i), 1.0));  // sum over out degree
    }
  }
  // inequalities for guaranteeing connectednes
  const double p = numberOfNodes;
  for (size_t j = 1; j < numberOfNodes; ++j) {
    for (size_t i = 1; i < j; ++i) {
      entries.push_back(Entry(index.constraintU(i, j), index.variableU(i), 1.0));   // +u_i
      entries.push_back(Entry(index.constraintU(i, j), index.variableU(j), -1.0));  // -u_j
      entries.push_back(Entry(index.constraintU(i, j), index.variableX(i, j), p));  // +px_ij
    }
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      entries.push_back(Entry(index.constraintU(i, j), index.variableU(i), 1.0));   // +u_i
      entries.push_back(Entry(index.constraintU(i, j), index.variableU(j), -1.0));  // -u_j
      entries.push_back(Entry(index.constraintU(i, j), index.variableX(i, j), p));  // +px_ij
    }
  }
}

static void setCConstraints(std::vector<Entry>& entries, const Index& index, const graph::Euclidean& euclidean) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      const double dist = euclidean.weight(i, j);
      entries.push_back(Entry(index.constraintC(i, j), index.variableX(i, j), -dist));
      entries.push_back(Entry(index.constraintC(i, j), index.variableC(), 1.0));
      entries.push_back(Entry(index.constraintC(j, i), index.variableX(j, i), -dist));  // exploit symmetry
      entries.push_back(Entry(index.constraintC(j, i), index.variableC(), 1.0));
    }
  }
}

static size_t setAntiCrossingConstraints(std::vector<Entry>& entries, const Index& index, const graph::Euclidean& euclidean) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  size_t row                 = index.numConstraints();
  for (size_t i = 0; i < numberOfNodes; ++i) {
    for (size_t j = i + 1; j < numberOfNodes; ++j) {
      for (size_t k = i + 1; k < numberOfNodes; ++k) {
        if (k == j) {
          continue;
        }
        for (size_t l = k + 1; l < numberOfNodes; ++l) {
          if (l == j) {
            continue;
          }
          if (graph::intersect(graph::LineSegment{euclidean.position(i), euclidean.position(j)},
                               graph::LineSegment{euclidean.position(k), euclidean.position(l)})) {
            entries.push_back(Entry(row, index.variableX(i, j), 1.0));  // Here we are putting 4 constraints into one
            entries.push_back(Entry(row, index.variableX(k, l), 1.0));  // by abusing the fact, that at most one of
            entries.push_back(Entry(row, index.variableX(j, i), 1.0));  // edges {(i,j), (j,i)} and at most on of
            entries.push_back(Entry(row, index.variableX(l, k), 1.0));  // {(k,l), (l,k)} can be part of the solution.
            ++row;
          }
        }
      }
    }
  }
  return row - index.numConstraints();
}

static void forbidCrossing(HighsModel& model, std::vector<Entry>& entries, const graph::Euclidean& euclidean, const Index& index) {
  const size_t numOfAntiCrossingConstraints = setAntiCrossingConstraints(entries, index, euclidean);
  setAntiCrossingBounds(model, numOfAntiCrossingConstraints);
  model.lp_.num_row_ += numOfAntiCrossingConstraints;
}

Result solve(const graph::Euclidean& euclidean, const ProblemType problemType, const bool noCrossing) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  const Index index(numberOfNodes, problemType);

  HighsModel model;
  model.lp_.num_col_ = index.numVariables();
  model.lp_.num_row_ = index.numConstraints();  // may be changed later on by forbidCrossing()
  model.lp_.sense_   = ObjSense::kMinimize;
  model.lp_.offset_  = 0;                       // offset has no effect on optimization

  std::vector<Entry> entries;
  if (problemType == ProblemType::BTSP_exact) {
    setBTSPcost(model, index);
    setMillerTuckerZemlinBounds(model, index, numberOfNodes);
    setMillerTuckerZemlinMatrix(entries, index, numberOfNodes);
    setCBounds(model, index);
    setCConstraints(entries, index, euclidean);
    if (noCrossing) {
      forbidCrossing(model, entries, euclidean, index);
    }
  }
  else if (problemType == ProblemType::BTSPP_exact) {
    setBTSPcost(model, index);
    setMillerTuckerZemlinBounds(model, index, numberOfNodes);
    setPathBounds(model, index);
    setMillerTuckerZemlinMatrix(entries, index, numberOfNodes);
    setCBounds(model, index);
    setCConstraints(entries, index, euclidean);
  }
  else if (problemType == ProblemType::TSP_exact) {
    entries.reserve((numberOfNodes - 1) * index.xConstraints() + 3 * index.uConstraints());
    setTSPcost(model, index, euclidean, numberOfNodes);          // set cost function
    setMillerTuckerZemlinBounds(model, index, numberOfNodes);    // set bounds on variables and constraints
    setMillerTuckerZemlinMatrix(entries, index, numberOfNodes);  // set left hand side of constraints
  }

  Eigen::SparseMatrix<double> A(model.lp_.num_row_, model.lp_.num_col_);
  A.setFromTriplets(entries.begin(), entries.end());

  // copy data from eigen sparse matrix into HiGHs sparse matrix
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;                                    // use column compressed storage order
  model.lp_.a_matrix_.start_.assign(A.outerIndexPtr(), A.outerIndexPtr() + A.cols());      // copy start indices of columns
  model.lp_.a_matrix_.start_.push_back(A.nonZeros());                                      // add number of nonZeros in the end
  model.lp_.a_matrix_.index_.assign(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros());  // copy inner indices
  model.lp_.a_matrix_.value_.assign(A.valuePtr(), A.valuePtr() + A.nonZeros());            // copy values

  Highs highs;
  highs.setOptionValue("output_flag", false);
  [[maybe_unused]] HighsStatus return_status = highs.passModel(model);
  assert(return_status == HighsStatus::kOk);

  return_status = highs.run();  // solve instance
  assert(return_status == HighsStatus::kOk);

  [[maybe_unused]] const HighsModelStatus& model_status = highs.getModelStatus();
  assert(model_status == HighsModelStatus::kOptimal);

  const HighsInfo& info         = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();  // get variables of optimal solution

  std::vector<size_t> tour(numberOfNodes);
  tour[0] = 0;  // circle starts by definition with node 0
  for (size_t i = 1; i < numberOfNodes; ++i) {
    // round because solution might not exactly hit integers
    tour[std::round(solution.col_value[index.variableU(i)]) + 1] = i;
  }

  if (problemType == ProblemType::BTSP_exact) {
    return Result{tour, info.objective_function_value, findBottleneck(euclidean, tour, true)};
  }
  else if (problemType == ProblemType::BTSPP_exact) {
    return Result{tour, info.objective_function_value, findBottleneck(euclidean, tour, false)};
  }
  else if (problemType == ProblemType::TSP_exact) {
    return Result{
        tour,
        info.objective_function_value,
        graph::Edge{0, 0}
    };
  }
  else {
    throw UnknownType("[SOLVE] Unknown problem type.");
  }
}
}  // namespace exactsolver
