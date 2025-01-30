#include "counter.h"

#include <chrono>
#include <fstream>

namespace GroupIP
{

  mpq_class cone_evaluation(const Matrix &A, const Vector &b, const Vector &c)
  {

    int_t dim = A.size();
    SNFClass snf(A);
    snf.calculate_snf();

    auto snf_diagonal = snf.get_diagonal();
    Matrix P = snf.get_P();

    auto g = calculate_group_elements(P, snf_diagonal, b);
    auto h = transpose(calculate_adjugate_matrix(A));

    Dynamic d(c, g, h);
    d.Init(snf_diagonal);
    auto generating_function_numerator = d(dim - 1, g[dim]);

    auto numerator_coeffs = generating_function_numerator.get_coeffs();
    auto alphas = generating_function_numerator.get_exps();
    auto betas = calculate_betas(c, g, h);

    alphas_normalize(alphas, A, b, c);
    betas_normalize(betas, A);

    ToddPoly todd_poly(dim, betas);
    todd_poly.init();
    auto todd = todd_poly.get_todd();

    auto cone_value = tailor_series_const_term(alphas, numerator_coeffs, betas, todd, dim);

    return cone_value;
  }

  mpq_class tailor_series_const_term(const std::vector<ExpPoly::exp_t> &alphas,
                                     const std::vector<ExpPoly::coeff_t> &numerator_coeffs,
                                     const Vector &betas,
                                     const std::vector<mpq_class> &todd,
                                     int_t dim)
  {
    mpq_class numerator = 0;
    mpq_class denominator = 1;
    auto factoials = calculate_factorials(dim);

    for (int i = 0; i < betas.size(); ++i)
    {
      denominator *= betas[i];
    }

    for (int i = 0; i < alphas.size(); ++i)
    {
      mpq_class alpha_part_sum = 0;
      mpz_class alpha_power = 1;
      for (int j = 0; j <= dim; ++j)
      {
        alpha_part_sum += (todd[dim - j] / factoials[j]) * alpha_power;
        alpha_power *= alphas[i];
      }
      numerator += (numerator_coeffs[i] * alpha_part_sum);
    }

    return numerator / denominator;
  }

  void alphas_normalize(std::vector<ExpPoly::exp_t> &alphas, const Matrix &A, const Vector &b, const Vector &c)
  {
    auto Aadj = calculate_adjugate_matrix(A);
    mpz_class det = abs(calculate_det(A));

    for (int i = 0; i < alphas.size(); ++i)
    {
      alphas[i] = (dot_product(c, matrix_dot_vector(Aadj, b)) + alphas[i]) / det;
    }
  }

  void betas_normalize(Vector &betas, const Matrix &A)
  {
    mpz_class det = abs(calculate_det(A));

    for (int i = 0; i < betas.size(); ++i)
    {
      betas[i] = betas[i] / det;
    }
  }

  std::vector<GroupIP::GroupElement> calculate_group_elements(const Matrix &P, const Vector &snf_diagonal,
                                                              const Vector &b)
  {

    int n = b.size();
    std::vector<GroupElement> g(n + 1, GroupElement(snf_diagonal));
    g[n].assign(matrix_dot_vector(P, b));

    for (int i = 0; i < n; ++i)
    {
      g[i].assign(get_matrix_column(P, i));
    }

    return g;
  }

  Vector calculate_betas(const Vector &c, std::vector<GroupElement> &g, const Matrix &h)
  {
    Vector betas;
    for (int i = 0; i < h.size(); ++i)
    {
      mpz_class sum = 0;
      for (int j = 0; j < h[0].size(); ++j)
      {
        sum += c[j] * h[i][j];
      }
      betas.push_back(g[i].getOrder() * sum);
    }
    return betas;
  }

}
