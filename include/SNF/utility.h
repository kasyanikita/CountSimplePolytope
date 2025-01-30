#pragma once

#include <vector>
#include <cassert>
#include <random>
#include <iostream>

#include "Eigen/Dense"

namespace LatLib
{
    using namespace Eigen;

    template <typename T>
    T euclid_bezout_nonnegative(
        const T &a,
        const T &b,
        T &u, T &v,
        T &u0, T &v0)
    {
        assert(a >= 0);
        assert(b >= 0);

        u = 1;
        u0 = 0;
        v = 0;
        v0 = 1;

        T r = a;
        T r0 = b;

        while (r0 != 0)
        {
            T q = r / r0;

            T aux = u0;
            u0 = u - q * u0;
            u = aux;

            aux = v0;
            v0 = v - q * v0;
            v = aux;

            aux = r0;
            r0 = r - q * r0;
            r = aux;
        }

        return r;
    }

    template <typename T>
    T euclid_bezout(
        const T &a,
        const T &b,
        T &u,
        T &v,
        T &u0,
        T &v0)
    {
        T res = euclid_bezout_nonnegative<T>(abs(a), abs(b), u, v, u0, v0);

        if (a < 0)
        {
            u = -u;
            u0 = -u0;
        }
        if (b < 0)
        {
            v = -v;
            v0 = -v0;
        }

        return res;
    }

    template <typename vecT>
    size_t find_nonZero_index(const vecT &vec, size_t start)
    {
        size_t n = vec.size();

        for (size_t i = start; i < n; ++i)
        {
            if (vec[i] != 0)
                return i;
        }

        return n;
    }

    template <typename matrixT>
    bool is_diagonal(const matrixT &A)
    {
        size_t n = A.cols();
        assert(n == A.rows());

        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (i != j && A(i, j) != 0)
                    return false;

        return true;
    }

} // namespace LatLib
