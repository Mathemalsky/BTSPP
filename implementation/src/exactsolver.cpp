#include "exactsolver.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include <Eigen/SparseCore>

#include <Highs.h>

#include "graph/graph.hpp"

using Entry = Eigen::Triplet<double>;

static constexpr double M_INFINITY = 1e32;

class Index {
public:
  Index(size_t numberOfNodes) : pNumberOfNodes(numberOfNodes) {}

  size_t variableX(const size_t i, const size_t j) const { return j * pNumberOfNodes + i; }
  size_t variableU(const size_t i) const { return i * pNumberOfNodes + i; }

  size_t constraintXin(const size_t j) const { return j; }
  size_t constraintXout(const size_t j) const { return pNumberOfNodes + j; }
  size_t constraintU(const size_t i, const size_t j) const {
    return 2 * pNumberOfNodes + (j - 1) * (pNumberOfNodes - 2) + (i > j ? i - 2 : i - 1);
  }

private:
  const size_t pNumberOfNodes;
};

std::vector<size_t> solveTSP(const Euclidean& euclidean) {
  const size_t numberOfNodes = euclidean.numberOfNodes();
  Index index(numberOfNodes);

  HighsModel model;
  model.lp_.num_col_ = numberOfNodes * numberOfNodes;
  model.lp_.num_row_ = 2 * numberOfNodes + (numberOfNodes - 1) * (numberOfNodes - 2);
  model.lp_.sense_   = ObjSense::kMinimize;
  model.lp_.offset_  = 0;  // offset has no effect on optimization

  // set costfunction
  model.lp_.col_cost_.resize(model.lp_.num_col_);
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      const double dist                          = euclidean.distance(i, j);
      model.lp_.col_cost_[index.variableX(i, j)] = dist;
      model.lp_.col_cost_[index.variableX(j, i)] = dist;  // exploiting symmetry
    }
    model.lp_.col_cost_[index.variableU(j)] = 0.0;  // here we store u_j instead of x_jj
  }

  model.lp_.col_lower_ = std::vector(model.lp_.num_col_, 0.0);  // set lower bound to zero

  // iterate over all variables
  const double bigM = numberOfNodes - 2;  //
  model.lp_.col_upper_.resize(model.lp_.num_col_);
  model.lp_.integrality_.resize(model.lp_.num_col_);
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      model.lp_.col_upper_[index.variableX(i, j)]   = 1.0;                     // bound by 1
      model.lp_.col_upper_[index.variableX(j, i)]   = 1.0;                     // exploiting symmetry
      model.lp_.integrality_[index.variableX(i, j)] = HighsVarType::kInteger;  // constrain x_ij to be \in {0,1}
      model.lp_.integrality_[index.variableX(j, i)] = HighsVarType::kInteger;  // exploiting symmetry
    }
    model.lp_.col_upper_[index.variableU(j)]   = bigM;
    model.lp_.integrality_[index.variableU(j)] = HighsVarType::kContinuous;
  }

  // iterate over all constraints
  model.lp_.row_lower_.resize(model.lp_.num_row_);
  model.lp_.row_upper_.resize(model.lp_.num_row_);
  // inequalities for fixing in and out degree to 1
  for (size_t i = 0; i < 2 * numberOfNodes; ++i) {
    model.lp_.row_lower_[i] = 1.0;
    model.lp_.row_upper_[i] = 1.0;
  }
  // inequalities for guaranteeing connectednes
  const double p = numberOfNodes;
  for (size_t i = 2 * numberOfNodes; i < (size_t) model.lp_.num_row_; ++i) {
    model.lp_.row_lower_[i] = -M_INFINITY;  // lower bound -infinity
    model.lp_.row_upper_[i] = p - 1;        // upper bound p-1
  }

  // construct matrix A
  // use handy fill function from eigen to populate the sparse matrix
  std::vector<Entry> entries;
  entries.reserve(2 * numberOfNodes * (numberOfNodes - 1) + 3 * (numberOfNodes - 1) * (numberOfNodes - 2));

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
  // std::cout << "Model status: " << highs.modelStatusToString(model_status) << std::endl;

  /*
  const HighsInfo& info = highs.getInfo();
  std::cout << "Simplex iteration count : " << info.simplex_iteration_count << std::endl;
  std::cout << "Objective function value: " << info.objective_function_value << std::endl;
  std::cout << "Primal  solution status : " << highs.solutionStatusToString(info.primal_solution_status) << std::endl;
  std::cout << "Dual    solution status : " << highs.solutionStatusToString(info.dual_solution_status) << std::endl;
  std::cout << "Basis                   : " << highs.basisValidityToString(info.basis_validity) << std::endl;
  */

  // Get the solution values and basis
  const HighsSolution& solution = highs.getSolution();
  // const HighsBasis& basis = highs.getBasis();
  const HighsLp& lp = highs.getLp();  // get a const reference to the LP data in HiGHS

  /*
  // Report the primal solution values
  for (int col = 0; col < lp.num_col_; col++) {
    std::cout << "Column " << col;
    if (info.primal_solution_status)
      std::cout << "; value = " << solution.col_value[col];
    std::cout << std::endl;
  }
  */
  /*
  for (int row = 0; row < lp.num_row_; row++) {
    std::cout << "Row    " << row;
    if (info.primal_solution_status)
      std::cout << "; value = " << solution.row_value[row];
    std::cout << std::endl;
  }
  */

  std::vector<size_t> tour(numberOfNodes);
  tour[0] = 0;  // circle starts by definition with node 0
  for (size_t i = 1; i < numberOfNodes; ++i) {
    tour[solution.col_value[index.variableU(i)] + 1] = i;
  }

  return tour;
}
