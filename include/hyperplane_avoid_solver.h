#ifndef HYPERPLANE_AVOID_SOLVER_H_
#define HYPERPLANE_AVOID_SOLVER_H_

#include "global_defs.h"
#include "matrix_tools.h"

namespace GroupIP
{

  class HyperplaneAvoidSolver
  {
  private:
    Vector c_;
    Matrix A_;
    Matrix get_edge_directions(const std::vector<Matrix> &A_cones);

  public:
    HyperplaneAvoidSolver(const Matrix &A) : A_(A) {}
    Vector get_vector(const std::vector<Matrix> &);
  };

  struct LinMultiplier
  {
    std::unordered_map<int, int> coeffs;
    int constant = 0;

    void substitution(int, int);
  };

  Vector choose_c(std::vector<LinMultiplier> &, int);
  int eliminate_variable(std::vector<LinMultiplier> &, int);
}

#endif // HYPERPLANE_AVOID_SOLVER_H_