#ifndef CANTOR_ZASSENHAUZ_TEST_H
#define CANTOR_ZASSENHAUZ_TEST_H

#include <vector>
#include "cantor_zassenhaus.h"

class Test {
public:
    static void run_random_tests(const CantorZassenhaus::PolynomialRing &Fx);

    static long measure(const CantorZassenhaus::PolynomialRing &Fx,
                        int degree);

private:
    Test() = default;

    static bool run_random_test(const CantorZassenhaus::PolynomialRing &Fx,
                                int degree);

    static std::vector<std::pair<CantorZassenhaus::Polynomial, int>>
    run_random(const CantorZassenhaus::PolynomialRing &Fx, int degree);

    static const std::vector<int> kDegrees;

    static bool test(const CantorZassenhaus::PolynomialRing &Fx,
                     CantorZassenhaus::Polynomial &f);
};

#endif //CANTOR_ZASSENHAUZ_TEST_H
