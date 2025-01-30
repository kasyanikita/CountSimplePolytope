#include "dynamic.h"

namespace GroupIP
{

  std::vector<std::vector<ExpPoly>> Dynamic::init_dp(const Vector &snf) const
  {
    mpz_class n_cols = 1;
    for (int i = 0; i < n_; ++i)
    {
      n_cols *= snf[i];
    }
    std::vector<std::vector<ExpPoly>> res(n_, std::vector<ExpPoly>(n_cols.get_si()));
    return res;
  }

  Dynamic::Dynamic(const Vector &c, const std::vector<GroupIP::GroupElement> &g, const Matrix &h)
      : c_(c), g_(g), h_(h), n_(h.size()) {}

  void Dynamic::Init(const Vector &snf)
  {
    dp = init_dp(snf);

    for (size_t i = 0; i < dp.size(); ++i)
    {
      isComputed.push_back(std::vector<bool>(dp[0].size(), false));
    }
  }

  ExpPoly &Dynamic::d(uint_t k, GroupElement g) { return dp[k][g.get_idx()]; }

  ExpPoly Dynamic::operator()(uint_t k, const GroupElement &ge)
  {
    auto g_idx = ge.get_idx();
    if (isComputed[k][g_idx])
    {
      return dp[k][g_idx];
    }
    else
    {
      if (k != 0)
      {
        auto sum = (*this)(k - 1, ge).monomial_multiply(
            -dot_product(c_, scalar_vector_product(0, h_[k])), 1);

        // j = 1 ... rk - 1
        for (int j = 1; j < g_[k].getOrder(); ++j)
        {
          auto sum_part = (*this)(k - 1, ge - j * g_[k])
                              .monomial_multiply(
                                  -dot_product(c_, scalar_vector_product(j, h_[k])), 1);
          sum = sum + sum_part;
        }
        // save d(k, ge)
        dp[k][g_idx] = sum;
        isComputed[k][g_idx] = true;
      }
      else
      {
        // k = 0
        auto s = calc_s(ge, g_[0], g_[0].getOrder());
        mpz_class exp = 0;
        mpz_class coeff = 0;
        if (s != -1)
        {
          exp = -dot_product(c_, scalar_vector_product(s, h_[0]));
          coeff = 1;
        }
        d(k, ge).init({exp}, {coeff});
        isComputed[k][g_idx] = true;
      }
    }
    return dp[k][g_idx];
  }
}