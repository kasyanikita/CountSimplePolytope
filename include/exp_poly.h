#ifndef COUNTINGINTEGERPOINTS_EXPPOLY_H_
#define COUNTINGINTEGERPOINTS_EXPPOLY_H_

#include <gmpxx.h>

#include "global_defs.h"

namespace GroupIP
{
    class ExpPoly
    {
    public:
        using exp_t = mpz_class;
        using coeff_t = mpz_class;

    private:
        std::unordered_map<exp_t, coeff_t, mpz_hash, mpz_equal> poly;

    public:
        friend ExpPoly operator*(coeff_t c, ExpPoly &exp_poly);
        friend ExpPoly operator*(ExpPoly &exp_poly, coeff_t c);
        friend std::ostream &operator<<(std::ostream &out, ExpPoly &exp_poly);

        ExpPoly() {}
        ExpPoly(const std::vector<exp_t> &exp, const std::vector<coeff_t> &coeff);
        ExpPoly(const std::unordered_map<exp_t, coeff_t, mpz_hash, mpz_equal> &);
        void init(const std::vector<exp_t> &exp, const std::vector<coeff_t> &coeff);
        ExpPoly operator+(const ExpPoly &rhs) const;
        ExpPoly operator*(const ExpPoly &rhs);
        ExpPoly monomial_multiply(exp_t exp, coeff_t coeff);
        std::vector<coeff_t> get_coeffs();
        std::vector<exp_t> get_exps();
        std::unordered_map<exp_t, coeff_t, mpz_hash, mpz_equal> get_poly() const;
    };

    ExpPoly operator*(ExpPoly::coeff_t c, ExpPoly &exp_poly);
    ExpPoly operator*(ExpPoly &exp_poly, ExpPoly::coeff_t c);
    std::ostream &operator<<(std::ostream &out, ExpPoly &exp_poly);

} // namespace GroupIP

#endif // COUNTINGINTEGERPOINTS_EXPPOLY_H_