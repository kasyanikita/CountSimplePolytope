#pragma once

#include <gmpxx.h>

#include "Eigen/Core"

namespace Eigen
{
    template <>
    struct NumTraits<mpq_class> : GenericNumTraits<mpq_class>
    {
        typedef mpq_class Real;
        typedef mpq_class NonInteger;
        typedef mpq_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 150,
            MulCost = 100
        };
    };

    template <>
    struct NumTraits<mpf_class> : GenericNumTraits<mpf_class>
    {
        typedef mpf_class Real;
        typedef mpf_class NonInteger;
        typedef mpf_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 50,
            MulCost = 75
        };
    };

    template <>
    struct NumTraits<mpz_class> : GenericNumTraits<mpz_class>
    {
        typedef mpz_class Real;
        typedef mpz_class NonInteger;
        typedef mpz_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 1,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 50,
            MulCost = 75
        };
    };
}