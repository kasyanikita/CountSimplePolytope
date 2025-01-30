#ifndef TODD_POLY_H_
#define TODD_POLY_H_

#include <vector>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <gmpxx.h>

#include "global_defs.h"
#include "tools.h"

namespace GroupIP
{
    class ToddPoly
    {
    protected:
        int_t m;
        std::vector<mpz_class> xi;
        std::vector<mpq_class> bernoulli;
        std::vector<mpq_class> todd;
        std::vector<mpq_class> todd_part;
        void update_part(mpz_class);
        void init_todd();
        void calc_todd();
        void calc_todd_arr();

    public:
        ToddPoly(int_t, const std::vector<mpz_class> &);
        void init();
        const std::vector<mpq_class> &get_todd() const;
        std::vector<mpq_class> &get_todd();
    };
}

#endif // TODD_POLY_H_