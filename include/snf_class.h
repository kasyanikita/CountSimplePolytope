#ifndef SNF_CLASS_H_
#define SNF_CLASS_H_

#include <vector>

#include "Eigen/Dense"
#include "SNF/DF_module.h"
#include "SNF/gmp_supp.h"
#include "global_defs.h"

using namespace GroupIP;

class SNFClass
{
public:
  explicit SNFClass(const Matrix &);

  void calculate_snf();

  const Matrix &get_A() const;
  const Matrix &get_P() const;
  const Matrix &get_S() const;
  const Matrix &get_Q() const;
  const Vector &get_diagonal() const;

private:
  using snf_int_type = mpz_class;
  using eigen_mat = Eigen::Matrix<snf_int_type, Eigen::Dynamic, Eigen::Dynamic>;
  using my_DF_algo = LatLib::my_diagonal_form_algo<snf_int_type>;

  Matrix A_;
  Matrix P_;
  Matrix S_;
  Matrix Q_;

  Matrix eigen2vector(eigen_mat);

  Vector s_diagonal_;

  void initialize_identity_matrices();
  void extract_diagonal();
};

#endif // SNF_CLASS_H_
