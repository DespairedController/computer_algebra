#include <chrono>

#include "test.h"
#include "cantor_zassenhaus.h"

const std::vector<int> Test::kDegrees({1, 2, 4, 8});

void Test::run_random_tests(const CantorZassenhaus::PolynomialRing &Fx) {
    for (auto d : kDegrees) {
        {
            if (!run_random_test(Fx, d)) {
                return;
            };
        }
    }
}

long Test::measure(const CantorZassenhaus::PolynomialRing &Fx, int degree) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    run_random(Fx, degree);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
}


bool Test::run_random_test(const CantorZassenhaus::PolynomialRing &Fx, int degree) {
    CantorZassenhaus::Polynomial P;
    auto gen = Givaro::GivRandom();
    Fx.random(gen, P, Givaro::Degree(degree));
    return test(Fx, P);
}

bool Test::test(const CantorZassenhaus::PolynomialRing &Fx, CantorZassenhaus::Polynomial &f) {
    Fx.write(std::cout << "Testing: ", f) << std::endl;
    auto factorization = CantorZassenhaus::run(Fx, f);
    CantorZassenhaus::Polynomial a, tmp1, tmp2;
    Fx.assign(a, Fx.getDomain().one);
    for (const auto &pair : factorization) {
        Fx.pow(tmp1, pair.first, pair.second);
        Fx.mul(tmp2, a, tmp1);
        a = tmp2;
    }
    if (a == f) {
        Fx.write(std::cout << "\033[1;32mTest passed:\033[0m", a) << std::endl;
        return true;
    } else {
        Fx.write(Fx.write(
                std::cout << "\033[1;31mTest failed:\033[0m\nExpected: ", f) << "\nBut got: ", a) << std::endl;
        return false;
    }
}

std::vector<std::pair<CantorZassenhaus::Polynomial, int>>
Test::run_random(const CantorZassenhaus::PolynomialRing &Fx, int degree) {
    CantorZassenhaus::Polynomial P;
    auto gen = Givaro::GivRandom();
    Fx.random(gen, P, Givaro::Degree(degree));
    return CantorZassenhaus::run(Fx, P);
}
