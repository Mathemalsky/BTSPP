#pragma once

enum class ProblemType : unsigned int {
  BTSP_approx,
  BTSP_exact,
  BTSPP_exact,
  TSP_exact,
  NUMBER_OF_OPTIONS
};

constexpr unsigned int SEED_LENGTH = 2;
