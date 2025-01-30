#ifndef DYNAMIC_H_
#define DYNAMIC_H_

#include <stdexcept>
#include <flint/fmpz_mat.h>

#include "global_defs.h"
#include "group_element.h"
#include "snf_class.h"
#include "matrix_tools.h"
#include "tools.h"

namespace GroupIP
{
    class Dynamic
    {
        int_t n_;
        Vector c_;
        Matrix h_;
        std::vector<GroupIP::GroupElement> g_;
        std::vector<std::vector<bool>> isComputed;
        std::vector<std::vector<ExpPoly>> init_dp(const Vector &) const;

    public:
        std::vector<std::vector<ExpPoly>> dp;
        Dynamic(const Vector &, const std::vector<GroupIP::GroupElement> &,
                const Matrix &);
        void Init(const Vector &);
        ExpPoly &d(uint_t k, GroupElement g);
        ExpPoly operator()(uint_t k, const GroupElement &ge);
    };
}

#endif // DYNAMIC_H_