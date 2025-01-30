#ifndef COUNTER_H_
#define COUNTER_H_

#include <chrono>
#include <iostream>

#include <gmpxx.h>
#include <flint/fmpz_mat.h>

#include "io_handler.h"
#include "dynamic.h"
#include "exp_poly.h"
#include "todd_poly.h"
#include "tools.h"
#include "matrix_tools.h"
#include "hyperplane_avoid_solver.h"

namespace GroupIP
{
    mpq_class cone_evaluation(const Matrix &A, const Vector &b, const Vector &c);
    Vector calculate_betas(const Vector &c, std::vector<GroupElement> &g, const Matrix &h);

    void alphas_normalize(std::vector<ExpPoly::exp_t> &alphas, const Matrix &A, const Vector &b, const Vector &c);
    void betas_normalize(Vector &b, const Matrix &A);

    std::vector<GroupIP::GroupElement> calculate_group_elements(
        const Matrix &P, const Vector &snf_diagonal, const Vector &b);

    mpq_class tailor_series_const_term(const std::vector<ExpPoly::exp_t> &alphas,
                                       const std::vector<ExpPoly::coeff_t> &numerator_coeffs,
                                       const Vector &betas,
                                       const std::vector<mpq_class> &todd,
                                       int_t dim);
}

#endif // COUNTER_H_