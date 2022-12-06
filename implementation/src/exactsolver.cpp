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
  model.lp_.num_row_ = 2 * numberOfNodes + (numberOfNodes - 1) * (numberOfNodes - 2) / 2;
  model.lp_.sense_   = ObjSense::kMinimize;
  model.lp_.offset_  = 0;  // offset has no effect on optimization

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
  model.lp_.col_upper_ = std::vector(model.lp_.num_col_, 1.0);  // set upper bound to zero // rework

  model.lp_.row_lower_ = std::vector(model.lp_.num_row_, 0.0);  // set lower bound to zero // rework
  model.lp_.row_upper_ = std::vector(model.lp_.num_row_, 0.0);  // set lower bound to zero // rework
}
