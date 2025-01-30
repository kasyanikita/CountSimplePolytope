#include "matrix_tools.h"

namespace GroupIP
{
    Vector scalar_vector_product(mpz_class s, const Vector &v)
    {
        Vector res(v.size());
        for (int i = 0; i < v.size(); ++i)
        {
            res[i] = s * v[i];
        }
        return res;
    }

    mpz_class dot_product(const Vector &a, const Vector &b)
    {
        mpz_class res = 0;
        for (int i = 0; i < a.size(); ++i)
        {
            res += a[i] * b[i];
        }
        return res;
    }

    Vector get_matrix_column(const Matrix &M, int k)
    {
        Vector res;
        for (int i = 0; i < M.size(); ++i)
        {
            res.push_back(M[i][k]);
        }
        return res;
    }

    Vector matrix_dot_vector(const Matrix &M, const Vector &v)
    {
        Vector res;
        for (int i = 0; i < M.size(); ++i)
        {
            mpz_class sum = 0;
            for (int j = 0; j < M[i].size(); ++j)
            {
                sum += M[i][j] * v[j];
            }
            res.push_back(sum);
        }
        return res;
    }

    int_t calculate_det(const GroupIP::Matrix &A)
    {
        fmpz_mat_t Af;
        fmpz_t det;
        fmpz_mat_init(Af, A.size(), A[0].size());
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A[i].size(); ++j)
            {
                auto val = fmpz_mat_entry(Af, i, j);
                *val = A[i][j].get_si();
            }
        }
        fmpz_mat_det(det, Af);
        return *det;
    }

    Matrix calculate_adjugate_matrix(const Matrix &A)
    {
        fmpz_mat_t Aadj;
        fmpz_t den;
        fmpz_t det;
        Matrix res(A.size(), Vector(A[0].size()));
        fmpz_mat_init(Aadj, A.size(), A[0].size());
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A[i].size(); ++j)
            {
                auto val = fmpz_mat_entry(Aadj, i, j);
                *val = A[i][j].get_si();
            }
        }

        fmpz_mat_det(det, Aadj);
        int is_not_singular = fmpz_mat_inv(Aadj, den, Aadj);
        if (is_not_singular == 0)
        {
            throw std::domain_error("Matrix A is singular");
        }

        fmpz_divexact(den, det, den);
        fmpz_mat_scalar_divexact_fmpz(Aadj, Aadj, den);

        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A[i].size(); ++j)
            {
                auto val = fmpz_mat_entry(Aadj, i, j);
                res[i][j] = fmpz_get_d(val);
            }
        }
        return res;
    }

    Matrix transpose(const Matrix &matrix)
    {
        Matrix matrix_transpose(matrix[0].size(), Vector(matrix.size(), 0));
        for (int i = 0; i < matrix.size(); ++i)
        {
            for (int j = 0; j < matrix[0].size(); ++j)
            {
                matrix_transpose[j][i] = matrix[i][j];
            }
        }

        return matrix_transpose;
    }
}