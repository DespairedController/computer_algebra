#ifndef CANTOR_ZASSENHAUZ_CANTOR_ZASSENHAUS_H
#define CANTOR_ZASSENHAUZ_CANTOR_ZASSENHAUS_H

#include <givaro/gfq.h>
#include <givaro/givpoly1.h>


class CantorZassenhaus {
 public:
    typedef Givaro::GFqDom<int> GaloisField;
    typedef Givaro::Poly1Dom<GaloisField, Givaro::Dense>::Element Polynomial;
    typedef Givaro::Poly1Dom<GaloisField, Givaro::Dense> PolynomialRing;

    static std::vector<std::pair<Polynomial, int>> run(
            const PolynomialRing &Fx,
            const Polynomial &f);

// private:
    CantorZassenhaus() = default;

    static Polynomial rooting(
            const PolynomialRing &Fx,
            const Polynomial &f);

    static Polynomial divide_by_lc(
            const PolynomialRing &Fx,
            const Polynomial &f);

    static Polynomial square_free(
            const PolynomialRing &Fx,
            const Polynomial &f);

    static std::vector<std::pair<Polynomial, int>>
    compute_distinct_degree(
            const PolynomialRing &Fx,
            const Polynomial &f);

    static std::vector<Polynomial> factorize_equal_degree(
            const PolynomialRing &Fx,
            const Polynomial &f, int degree);
};

#endif CANTOR_ZASSENHAUZ_CANTOR_ZASSENHAUS_H