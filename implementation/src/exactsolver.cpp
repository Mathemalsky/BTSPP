#include "exactsolver.hpp"

#include <Highs.h>

#include "graph/graph.hpp"
/*!
 * \brief solveExact solve the tsp to optimality using Miller-Tucker-Zemlin formulation
 * \param euclidean
 * \details x_ij belongs to the column j* numberOfNodes +i
 * hereby we ommit th x_jj entries and store u_j in that place
 */
void solveExact(const Euclidean* euclidean) {
  const size_t numberOfNodes = euclidean->numberOfNodes();

  HighsModel model;
  model.lp_.num_col_ = numberOfNodes * numberOfNodes;
  model.lp_.num_row_ = 2 * numberOfNodes + (numberOfNodes - 1) * (numberOfNodes - 2);
  model.lp_.sense_   = ObjSense::kMinimize;
  model.lp_.offset_  = 0;  // offset has no effect on optimization

  // set costfunction
  model.lp_.col_cost_.resize(model.lp_.num_col_);
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      const double dist                          = euclidean->distance(i, j);
      model.lp_.col_cost_[j * numberOfNodes + i] = dist;
      model.lp_.col_cost_[i * numberOfNodes + j] = dist;  // exploiting symmetry
    }
    model.lp_.col_cost_[j * numberOfNodes + j] = 0.0;  // here we store u_j instead of x_jj
  }

  model.lp_.col_lower_ = std::vector(model.lp_.num_col_, 0.0);  // set lower bound to zero

  // iterate over all variables
  const double bigM = numberOfNodes;  // can probably be shrinked to numberOfNodes -1
  model.lp_.col_upper_.resize(model.lp_.num_col_);
  model.lp_.integrality_.resize(model.lp_.num_col_);
  for (size_t j = 0; j < numberOfNodes; ++j) {
    for (size_t i = j + 1; i < numberOfNodes; ++i) {
      model.lp_.col_upper_[j * numberOfNodes + i]   = 1.0;                     // bound by 1
      model.lp_.col_upper_[i * numberOfNodes + j]   = 1.0;                     // exploiting symmetry
      model.lp_.integrality_[j * numberOfNodes + i] = HighsVarType::kInteger;  // constrain x_ij to be \in {0,1}
      model.lp_.integrality_[i * numberOfNodes + j] = HighsVarType::kInteger;  // exploiting symmetry
    }
    model.lp_.col_upper_[j * numberOfNodes + j]   = bigM;
    model.lp_.integrality_[j * numberOfNodes + j] = HighsVarType::kContinuous;
  }

  // iterate over all constraints
  model.lp_.row_lower_.resize(model.lp_.num_row_);
  model.lp_.row_upper_.resize(model.lp_.num_row_);
  for (size_t i = 0; i < 2 * numberOfNodes; ++i) {
    // fix the sum over the x_i to 1
    model.lp_.row_lower_[i] = 1.0;
    model.lp_.row_upper_[i] = 1.0;
  }
  const double p = numberOfNodes;
  for (size_t j = 0; j < numberOfNodes - 1; ++j) {
    for (size_t i = j + 1; i < numberOfNodes - 1; ++i) {
      model.lp_.row_upper_[2 * numberOfNodes + j * (numberOfNodes - 1) + i] = p - 1;
      model.lp_.row_upper_[2 * numberOfNodes + i * (numberOfNodes - 1) + j] = p - 1;  // exploiting symmetry
    }
  }

  model.lp_.row_lower_ = std::vector(model.lp_.num_row_, -1.0e32);  // set lower bound to -infinty

  // construct matrix A
  // DEFICIENT
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;         // use column compressed storage order
  for (size_t j = 0; j < (size_t) model.lp_.num_col_; ++j) {    // iterate through columns
    for (size_t i = 0; i < (size_t) model.lp_.num_row_; ++i) {  // iterate through rows
    }
  }
  // model.lp_.a_matrix_.start_;
}
