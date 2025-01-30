#ifndef TOOLS_H_
#define TOOLS_H_

#include <stdexcept>
#include <vector>

#include <flint/fmpz_mat.h>
#include <gmpxx.h>

#include "global_defs.h"
#include "group_element.h"
#include "exp_poly.h"

using namespace GroupIP;

namespace GroupIP
{
    mpz_class calc_s(const GroupElement &g, const GroupElement &g0, mpz_class r0);
    std::vector<mpz_class> calculate_factorials(int_t n);
    std::vector<mpz_class> calc_pascal(size_t);
    std::vector<mpq_class> calc_bernoulli(int_t n);
}

#endif // TOOLS_H_