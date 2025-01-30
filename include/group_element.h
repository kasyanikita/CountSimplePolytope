#ifndef COUNTINGINTEGERPOINTS_GROUPELEMENT_H_
#define COUNTINGINTEGERPOINTS_GROUPELEMENT_H_

#include "global_defs.h"

namespace GroupIP
{
    inline mpz_class modulo(mpz_class x, mpz_class mod);

    class GroupElement
    {
        friend GroupElement operator*(int_t c, const GroupElement &x);
        friend GroupElement operator*(const GroupElement &x, int_t c);
        friend std::ostream &operator<<(std::ostream &, GroupElement &);

        Vector _components;
        Vector _mod;
        mpz_class order = -1;

        void normalize_components();
        bool compare(const GroupElement &ge);

    public:
        GroupElement(Vector mod) : _mod(std::move(mod)) {}
        void assign(Vector comp);
        GroupElement invert() const;
        bool operator==(const GroupElement &ge);
        bool operator!=(const GroupElement &ge);
        GroupElement operator-(const GroupElement &rhs) const;
        GroupElement &operator+=(const GroupElement &rhs);
        GroupElement operator+(const GroupElement &rhs) const;
        bool operator<(const GroupElement &rhs) const;
        mpz_class getOrder();
        size_t get_idx() const;
        const Vector &get_components() const;
        const Vector &get_mod() const;
    };

    GroupElement operator*(int_t c, const GroupElement &x);
    GroupElement operator*(const GroupElement &x, int_t c);
    std::ostream &operator<<(std::ostream &, GroupElement &);

} // namespace GroupIP

#endif // COUNTINGINTEGERPOINTS_GROUPELEMENT_H_
