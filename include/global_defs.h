#ifndef _GLOBAL_DEFS_H_
#define _GLOBAL_DEFS_H_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <valarray>
#include <cassert>
#include <functional>
#include <gmpxx.h>

namespace GroupIP
{

    using int_t = int64_t;
    using uint_t = uint64_t;

    using Matrix = std::vector<std::vector<mpz_class>>;
    using Vector = std::vector<mpz_class>;
    using uVector = std::vector<uint_t>;

    struct mpz_hash {
      std::size_t operator()(const mpz_class& z) const {
        std::string str = z.get_str();
        return std::hash<std::string>()(str);
      }
    };

    struct mpz_equal {
      bool operator()(const mpz_class& lhs, const mpz_class& rhs) const {
        return lhs == rhs;
      }
    };

} // namespace GroupIP

#endif
