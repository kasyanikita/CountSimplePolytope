#include "hyperplane_avoid_solver.h"

namespace GroupIP
{

  void LinMultiplier::substitution(int index, int value)
  {
    if (coeffs.find(index) != coeffs.end())
    {
      auto coeff = coeffs[index];
      constant += coeff * value;
      coeffs.erase(index);
    }
  }

  Vector choose_c(Matrix &hyperplanes)
  {
    std::vector<LinMultiplier> lin_mult_vec;
    int n = hyperplanes[0].size();
    for (auto hyperplane : hyperplanes)
    {
      LinMultiplier lin_mult;
      for (int i = 0; i < hyperplane.size(); ++i)
      {
        if (hyperplane[i] != 0)
        {
          lin_mult.coeffs[i] = hyperplane[i].get_si();
        }
      }
      lin_mult_vec.push_back(lin_mult);
    }

    auto res = choose_c(lin_mult_vec, n);
    return res;
  }

  Vector choose_c(std::vector<LinMultiplier> &hyperplanes, int n)
  {
    Vector c;
    for (int i = 0; i < n; ++i)
    {
      auto val = eliminate_variable(hyperplanes, i);
      c.push_back(val);
    }
    return c;
  }

  int eliminate_variable(std::vector<LinMultiplier> &hyperplanes, int index)
  {
    std::unordered_set<int> bans;
    int res_val;

    for (int i = 0; i < hyperplanes.size(); ++i)
    {
      if ((hyperplanes[i].coeffs.size() == 1) &&
          (hyperplanes[i].coeffs.find(index) != hyperplanes[i].coeffs.end()))
      {
        double value = -hyperplanes[i].constant / hyperplanes[i].coeffs[index];
        if (fabs(trunc(value) - value) < 1e-3)
        {
          bans.insert(value);
        }
      }
    }

    for (int val = 0; val < hyperplanes.size() + 1; ++val)
    {
      if (bans.find(val) == bans.end())
      {
        res_val = val;
        break;
      }
      if (bans.find(-val) == bans.end())
      {
        res_val = -val;
        break;
      }
    }

    for (int i = 0; i < hyperplanes.size(); ++i)
    {
      hyperplanes[i].substitution(index, res_val);
    }

    return res_val;
  }

  Vector HyperplaneAvoidSolver::get_vector(const std::vector<Matrix> &Asubs)
  {
    Matrix h = get_edge_directions(Asubs);
    Matrix B = Asubs[0];
    Matrix Bh;

    for (int i = 0; i < h.size(); ++i)
    {
      Vector tmp_vector(h[0].size(), 0);
      for (int m = 0; m < B.size(); ++m)
      {
        for (int n = 0; n < B.size(); ++n)
        {
          tmp_vector[m] = tmp_vector[m] + B[m][n] * h[i][n];
        }
      }
      Bh.push_back(tmp_vector);
    }

    Vector z = choose_c(Bh);
    Vector c(B.size(), 0);
    for (int i = 0; i < c.size(); ++i)
    {
      for (int k = 0; k < c.size(); ++k)
      {
        c[i] = c[i] + B[k][i] * z[k];
      }
    }

    return c;
  }

  Matrix HyperplaneAvoidSolver::get_edge_directions(const std::vector<Matrix> &A_cones)
  {
    Matrix h;
    for (auto A_cone : A_cones)
    {
      auto A_adj = calculate_adjugate_matrix(A_cone);
      auto A_adj_transpose = transpose(A_adj);
      for (auto x : A_adj_transpose)
      {
        h.push_back(x);
      }
    }
    return h;
  }

} // namespace GroupIP