#include "exactsolver.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include <Eigen/SparseCore>

#include <Highs.h>

#include "draw/definitions.hpp"

#include "graph/graph.hpp"

using Entry = Eigen::Triplet<double>;

static constexpr double M_INFINITY = 1e32;

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
  size_t cConstraints() const { return (pType == ProblemType::BTSP_exact ? xVariables() : 0); }
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

static void setTSPcost(HighsModel& model, const Index& index, const Euclidean& euclidean, const size_t numberOfNodes) {
  model.lp_.col_cost_.resize(model.lp_.num_col_);
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      const double dist                          = euclidean.distance(i, j);
      model.lp_.col_cost_[index.variableX(i, j)] = dist;
      model.lp_.col_cost_[index.variableX(j, i)] = dist;  // exploiting symmetry
    }
    model.lp_.col_cost_[index.variableU(j)] = 0.0;  // here we store u_j instead of x_jj
  }
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
    model.lp_.col_upper_[index.variableU(j)]   = p - 2;                      // set upper bound of u_ij to p-2
    model.lp_.integrality_[index.variableU(j)] = HighsVarType::kContinuous;  // allow real numbers for u_ij
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

static void setBTSPcost(HighsModel& model, const Index& index) {
  model.lp_.col_cost_                    = std::vector<double>(model.lp_.num_col_, 0.0);  // set vector to 0
  model.lp_.col_cost_[index.variableC()] = 1.0;                                           // objective is c
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

static void setCConstraints(std::vector<Entry>& entries, const Index& index, const Euclidean& euclidean) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      const double dist = euclidean.distance(i, j);
      entries.push_back(Entry(index.constraintC(i, j), index.variableX(i, j), -dist));
      entries.push_back(Entry(index.constraintC(i, j), index.variableC(), 1.0));
      entries.push_back(Entry(index.constraintC(j, i), index.variableX(j, i), -dist));  // exploit symmetry
      entries.push_back(Entry(index.constraintC(j, i), index.variableC(), 1.0));
    }
  }
}

std::vector<unsigned int> solve(const Euclidean& euclidean, const ProblemType problemType) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  const Index index(numberOfNodes, problemType);

  HighsModel model;
  model.lp_.num_col_ = index.numVariables();
  model.lp_.num_row_ = index.numConstraints();
  model.lp_.sense_   = ObjSense::kMinimize;
  model.lp_.offset_  = 0;  // offset has no effect on optimization

  std::vector<Entry> entries;
  if (problemType == ProblemType::BTSP_exact) {
    setBTSPcost(model, index);
    setMillerTuckerZemlinBounds(model, index, numberOfNodes);
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
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;  // use column compressed storage order
  model.lp_.a_matrix_.start_.assign(A.outerIndexPtr(), A.outerIndexPtr() + A.cols());  // copy start indeces of columns
  model.lp_.a_matrix_.start_.push_back(A.nonZeros());  // add number of nonZeros in the end
  model.lp_.a_matrix_.index_.assign(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros());  // copy inner indeces
  model.lp_.a_matrix_.value_.assign(A.valuePtr(), A.valuePtr() + A.nonZeros());            // copy values

  Highs highs;
  HighsStatus return_status = highs.passModel(model);
  assert(return_status == HighsStatus::kOk);

  return_status = highs.run();  // solve instance
  assert(return_status == HighsStatus::kOk);

  const HighsModelStatus& model_status = highs.getModelStatus();
  assert(model_status == HighsModelStatus::kOptimal);
  //  std::cout << "Model status: " << highs.modelStatusToString(model_status) << std::endl;

  const HighsInfo& info = highs.getInfo();
  std::cout << "Simplex iteration count : " << info.simplex_iteration_count << std::endl;
  std::cout << "Objective function value: " << info.objective_function_value << std::endl;
  std::cout << "Primal  solution status : " << highs.solutionStatusToString(info.primal_solution_status) << std::endl;
  std::cout << "Dual    solution status : " << highs.solutionStatusToString(info.dual_solution_status) << std::endl;
  std::cout << "Basis                   : " << highs.basisValidityToString(info.basis_validity) << std::endl;

  // Get the solution values and basis
  const HighsSolution& solution = highs.getSolution();
  // const HighsBasis& basis = highs.getBasis();
  // const HighsLp& lp = highs.getLp();  // get a const reference to the LP data in HiGHS

  // Report the primal solution values
  /*
  for (int col = 0; col < lp.num_col_; col++) {
    std::cout << "Column " << col;
    if (info.primal_solution_status)
      std::cout << "; value = " << solution.col_value[col];
    std::cout << std::endl;
  }

  for (int row = 0; row < lp.num_row_; row++) {
    std::cout << "Row    " << row;
    if (info.primal_solution_status)
      std::cout << "; value = " << solution.row_value[row];
    std::cout << std::endl;
  }
  */

  std::vector<unsigned int> tour(numberOfNodes);
  tour[0] = 0;  // circle starts by definition with node 0
  for (unsigned int i = 1; i < numberOfNodes; ++i) {
    // round because solution might not exactly hit integers
    tour[std::round(solution.col_value[index.variableU(i)]) + 1] = i;
  }

  return tour;
}
