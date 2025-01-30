#include "tools.h"

#include <fstream>
#include <iomanip>
#include <chrono>

namespace GroupIP
{

  mpz_class calc_s(const GroupElement &g, const GroupElement &g0, mpz_class r0)
  {
    mpz_class res = -1;
    for (int i = 0; i < r0; ++i)
    {
      if (i * g0 == g)
      {
        res = i;
        break;
      }
    }
    return res;
  }

  std::vector<mpz_class> calculate_factorials(int_t n)
  {
    std::vector<mpz_class> fact(n + 1, 1);
    for (int i = 2; i <= n; ++i)
    {
      fact[i] = fact[i - 1] * mpz_class(i);
    }
    return fact;
  }

  std::vector<mpq_class> calc_bernoulli(int_t n)
  {
    // Calculate bernoulli numbers from 0 to n

    std::vector<mpq_class> bernoulli;
    bernoulli.push_back(1);
    bernoulli.push_back(-0.5);
    mpq_class sum = 0;
    for (size_t i = 2; i <= n; ++i)
    {
      if (i % 2 == 1)
      {
        bernoulli.push_back(0);
      }
      else
      {
        sum = 0;
        auto pascal = calc_pascal(i + 1);
        for (size_t k = 0; k < i; ++k)
        {
          mpq_class x = pascal[k + 2] * bernoulli[i - k - 1];
          sum += x;
        }
        bernoulli.push_back(-sum / (i + 1));
      }
    }

    return bernoulli;
  }

  std::vector<mpz_class> calc_pascal(size_t n)
  {
    // Get Pascal's triangle n-th row
    std::vector<mpz_class> res(n + 1, 1);
    for (size_t i = 0; i < n; ++i)
    {
      res[i + 1] = res[i] * (n - i) / (i + 1);
    }
    return res;
  }
}
