#include "todd_poly.h"

namespace GroupIP
{
  ToddPoly::ToddPoly(int_t _m, const std::vector<mpz_class> &_xi)
      : m(_m), xi(_xi), todd(_m + 1), todd_part(_m + 1) {}

  void ToddPoly::init()
  {
    bernoulli = calc_bernoulli(m);
    calc_todd_arr();
  }

  void ToddPoly::init_todd()
  {
    todd = todd_part;
  }

  void ToddPoly::calc_todd_arr()
  {
    // Calculate todd polynomials of degree from 0 to m
    update_part(xi[0]);
    init_todd();
    for (size_t i = 1; i < xi.size(); ++i)
    {
      update_part(xi[i]);
      calc_todd();
    }
  }

  void ToddPoly::update_part(mpz_class x)
  {
    todd_part[0] = 1;
    mpz_class fact = 1;
    mpz_class pow_x = 1;
    for (size_t i = 1; i <= m; ++i)
    {
      fact *= i;
      pow_x *= -x;
      todd_part[i] = pow_x * bernoulli[i] / fact;
    }
  }

  void ToddPoly::calc_todd()
  {
    fmpq_poly_t a, b, res;
    size_t size = todd.size();
    fmpq_poly_init2(a, size);
    fmpq_poly_init2(b, size);
    fmpq_poly_init2(res, size);
    for (size_t i = 0; i < size; ++i)
    {
      fmpq_poly_set_coeff_mpq(a, i,
                              mpq_class(todd[i]).get_mpq_t());
      fmpq_poly_set_coeff_mpq(
          b, i, mpq_class(todd_part[i]).get_mpq_t());
    }
    fmpq_poly_mullow(res, a, b, size);
    mpq_t tmp;
    mpq_init(tmp);
    for (size_t i = 0; i < size; ++i)
    {
      fmpq_poly_get_coeff_mpq(tmp, res, i);
      todd[i] = mpq_class(tmp);
    }
  }

  const std::vector<mpq_class> &ToddPoly::get_todd() const
  {
    return todd;
  }

  std::vector<mpq_class> &ToddPoly::get_todd()
  {
    return todd;
  }
}