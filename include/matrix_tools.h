#ifndef MATRIX_TOOLS_H_
#define MATRIX_TOOLS_H_

#include <flint/fmpz_mat.h>

#include "global_defs.h"

namespace GroupIP
{
    Vector scalar_vector_product(mpz_class s, const Vector &v);
    mpz_class dot_product(const Vector &a, const Vector &b);
    Vector get_matrix_column(const Matrix &M,
                             int k);
    Vector matrix_dot_vector(const Matrix &M, const Vector &v);
    int_t calculate_det(const GroupIP::Matrix &A);
    Matrix calculate_adjugate_matrix(const Matrix &A);
    Matrix transpose(const Matrix &A);
}

#endif // MATRIX_TOOLS_H_