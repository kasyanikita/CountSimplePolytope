#include "snf_class.h"

#include <iostream>

SNFClass::SNFClass(const Matrix &A) : A_(A)
{
  initialize_identity_matrices();
}

void SNFClass::calculate_snf()
{
  int size = A_.size();
  eigen_mat A(size, size);
  for (int i = 0; i < size; ++i)
  {
    for (int j = 0; j < size; ++j)
    {
      A(i, j) = A_[i][j];
    }
  }

  eigen_mat P, S, Q;
  my_DF_algo DF_algo;
  std::tie(P, S, Q) = DF_algo.compute(A);

  if (P * A * Q != S)
  {
    std::cerr << "Error: P * A * Q != S";
    std::exit(1);
  }

  P_ = eigen2vector(P);
  Q_ = eigen2vector(Q);
  S_ = eigen2vector(S);

  extract_diagonal();
}

Matrix SNFClass::eigen2vector(eigen_mat matrix)
{
  int size = matrix.rows();
  Matrix res(size, Vector(size));

  for (int i = 0; i < size; ++i)
  {
    for (int j = 0; j < size; ++j)
    {
      res[i][j] = matrix(i, j);
    }
  }

  return res;
}

const Matrix &SNFClass::get_A() const { return A_; }
const Matrix &SNFClass::get_P() const { return P_; }
const Matrix &SNFClass::get_S() const { return S_; }
const Matrix &SNFClass::get_Q() const { return Q_; }
const Vector &SNFClass::get_diagonal() const
{
  return s_diagonal_;
}

void SNFClass::initialize_identity_matrices()
{
  if (!S_.empty())
  {
    size_t size = S_.size();
    P_ = Matrix(size, Vector(size, 0));
    Q_ = Matrix(size, Vector(size, 0));

    for (size_t i = 0; i < size; ++i)
    {
      P_[i][i] = 1;
      Q_[i][i] = 1;
    }
  }
}

void SNFClass::extract_diagonal()
{
  s_diagonal_.clear();
  for (size_t i = 0; i < S_.size(); ++i)
  {
    s_diagonal_.push_back(S_[i][i]);
  }
}