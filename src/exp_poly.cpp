#include "exp_poly.h"

namespace GroupIP
{

  ExpPoly::ExpPoly(const std::vector<exp_t> &exp, const std::vector<coeff_t> &coeff)
  {
    assert(exp.size() == coeff.size());
    for (size_t i = 0; i < exp.size(); ++i)
    {
      auto it = poly.find(exp[i]);
      if (it != poly.end())
      {
        poly[exp[i]] += coeff[i];
      }
      else
      {
        poly.insert(std::make_pair(exp[i], coeff[i]));
      }
    }
  }

  ExpPoly::ExpPoly(
      const std::unordered_map<exp_t, coeff_t, mpz_hash, mpz_equal> &poly_)
  {
    poly = poly_;
  }

  void ExpPoly::init(const std::vector<exp_t> &exp, const std::vector<coeff_t> &coeff)
  {
    assert(exp.size() == coeff.size());
    for (size_t i = 0; i < exp.size(); ++i)
    {
      auto it = poly.find(exp[i]);
      if (it != poly.end())
      {
        poly[exp[i]] += coeff[i];
      }
      else
      {
        poly.insert(std::make_pair(exp[i], coeff[i]));
      }
    }
  }

  ExpPoly ExpPoly::operator+(const ExpPoly &rhs) const
  {
    ExpPoly res(*this);
    for (auto [k, v] : rhs.poly)
    {
      auto it = res.poly.find(k);
      if (it != res.poly.end())
      {
        res.poly[k] += v;
      }
      else
      {
        res.poly[k] = v;
      }
    }
    return res;
  }

  ExpPoly ExpPoly::operator*(const ExpPoly &rhs)
  {
    ExpPoly res;
    for (auto monoIter = rhs.poly.begin(); monoIter != rhs.poly.end();
         ++monoIter)
    {
      ExpPoly partRes =
          this->monomial_multiply(monoIter->first, monoIter->second);
      res = res + partRes;
    }
    return res;
  }

  ExpPoly ExpPoly::monomial_multiply(exp_t exp, coeff_t coeff)
  {
    std::vector<exp_t> exps;
    std::vector<coeff_t> coeffs;

    std::unordered_map<exp_t, coeff_t, mpz_hash, mpz_equal> new_poly;
    for (const auto &[e, c] : poly)
    {
      exps.push_back(e + exp);
      coeffs.push_back(c);
    }

    for (int i = 0; i < exps.size(); ++i)
    {
      new_poly[exps[i]] = coeffs[i];
    }

    return ExpPoly(new_poly);
  }

  std::vector<ExpPoly::coeff_t> ExpPoly::get_coeffs()
  {
    std::vector<coeff_t> coeffs;
    for (const auto &[exp, coeff] : poly)
    {
      coeffs.push_back(coeff);
    }
    return coeffs;
  }

  std::vector<ExpPoly::exp_t> ExpPoly::get_exps()
  {
    std::vector<exp_t> exps;
    for (const auto &[exp, coeff] : poly)
    {
      exps.push_back(exp);
    }
    return exps;
  }

  std::unordered_map<ExpPoly::exp_t, ExpPoly::coeff_t, mpz_hash, mpz_equal> ExpPoly::get_poly() const
  {
    return poly;
  }

  ExpPoly operator*(ExpPoly::coeff_t c, ExpPoly &exp_poly)
  {
    ExpPoly res(exp_poly);
    for (auto [k, v] : res.poly)
    {
      res.poly[k] *= c;
    }
    return res;
  }

  ExpPoly operator*(ExpPoly &exp_poly, ExpPoly::coeff_t c)
  {
    return c * exp_poly;
  }

  std::ostream &operator<<(std::ostream &out, ExpPoly &exp_poly)
  {
    int i = 0;
    for (auto [exp, coeff] : exp_poly.poly)
    {
      ++i;
      if ((coeff == 1) && (exp == 1))
      {
        out << "e^x";
      }
      else if ((coeff == 1) && (exp != 1))
      {
        out << "e^" << exp << "x";
      }
      else if (((coeff != 1) && (exp == 1)))
      {
        out << coeff << "*e^x";
      }
      else
      {
        out << coeff << "*e^" << exp << "x";
      }
      if (i != exp_poly.poly.size())
      {
        out << "+";
      }
    }
    return out;
  }

}