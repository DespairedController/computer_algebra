#include "cantor_zassenhaus.h"
#include <givaro/givpower.h>
#include <set>
#include <random>

CantorZassenhaus::Polynomial CantorZassenhaus::rooting(const PolynomialRing &Fx, const Polynomial &f) {
    const GaloisField &R = Fx.getDomain();
    long p = R.characteristic();
    long q = R.cardinality();
    uint64_t n = Fx.degree(f).value();
    Polynomial g;
    Fx.init(g);
    uint64_t i = 0;

    for (uint64_t c = 0; c < n + 1; c++) {
        GaloisField::Element coefficient;
        Fx.getEntry(coefficient, Givaro::Degree(c), f);
        if (R.isnzero(coefficient)) {
            Polynomial el;
            GaloisField::Element new_coefficient;
            Givaro::dom_power(new_coefficient, coefficient, q / p, R);
            Fx.init(el, Givaro::Degree(i / p), new_coefficient);
            Fx.addin(g, el);
        }
        i += 1;
    }
    return g;
}

CantorZassenhaus::Polynomial
CantorZassenhaus::square_free(const CantorZassenhaus::PolynomialRing &Fx, const CantorZassenhaus::Polynomial &f) {
    Polynomial ans, d, w, derivative;
    Fx.init(ans, Fx.one);
    Fx.diff(derivative, f);
//    Fx.write(std::cout << "derivative is", derivative) << std::endl;
    Fx.gcd(d, f, derivative);
//    Fx.write(std::cout << "gcd is", d) << std::endl;

    Fx.div(w, f, d);
//    Fx.write(std::cout << "w is", w) << std::endl;

    while (w != Fx.one) {
        Polynomial y, factor;
        Fx.gcd(y, w, d);
        Fx.div(factor, w, y);
        Fx.mulin(ans, factor);
        w = y;
        Fx.divin(d, y);
    }
    if (d != Fx.one) {
        d = rooting(Fx, d);
        auto sf = square_free(Fx, d);
        Fx.mulin(ans, sf);
    }
    return ans;
}

std::vector<std::pair<CantorZassenhaus::Polynomial, int>>
CantorZassenhaus::compute_distinct_degree(const CantorZassenhaus::PolynomialRing &Fx,
                                          const CantorZassenhaus::Polynomial &f) {
    std::vector<std::pair<Polynomial, int>> ans;
    const GaloisField &R = Fx.getDomain();
    uint64_t q = R.cardinality();
    uint64_t i = 1;
    Polynomial v, copy;
    copy = f;
    Fx.init(v, Givaro::Degree(1));
    while (2 * i <= Fx.degree(copy).value()) {
        Polynomial tmp, g, lc, tmp1;
        // tmp = v ^ q ^ i
        Fx.pow(tmp, v, q);
        // tmp = v ^ q ^ i - variable
        Fx.subin(tmp, v);
        Fx.mod(tmp1, tmp, copy);
        // g = gcd(copy, v^q^i - v)
        Fx.gcd(g, copy, tmp1);
        // g = g / g.lc()
//        g = divide_by_lc(Fx, g);

        if (g != Fx.one) {
            ans.emplace_back(g, i);
            //copy = copy / g
            Fx.divin(copy, g);
        }
        i += 1;
        q *= R.cardinality();
    }
    if (copy != Fx.one) {
//        copy = divide_by_lc(Fx, copy);
        ans.emplace_back(copy, Fx.degree(copy).value());
    }
    if (ans.empty()) {
        return {{f, Fx.degree(f).value()}};
    }
    return ans;
}

std::vector<CantorZassenhaus::Polynomial>
CantorZassenhaus::factorize_equal_degree(const PolynomialRing &Fx, const Polynomial &f, int degree) {
    const GaloisField &R = Fx.getDomain();
    uint64_t q = R.cardinality();
    uint64_t n = Fx.degree(f).value();
    if (n == degree) {
        return {f};
    }
    uint64_t r = n / degree;
    std::set<Polynomial> ans, adding, sub, tmp_set;
    ans.insert(f);

    std::random_device device;
    std::mt19937_64 generator(device());
    std::uniform_int_distribution<int> distribution(1, n - 1);
    auto gen = Givaro::GivRandom();

    while (ans.size() < r) {
        Polynomial a, b, c;
        Fx.random(gen, a, Givaro::Degree(distribution(generator)));
        //b = (a^((q^s - 1) // 2) - 1).quo_rem(f)[1]
        Givaro::IntegerDom ID;
        Givaro::IntegerDom::Element el, tmp;
        ID.assign(tmp, Givaro::Integer(q));
        ID.pow(el, tmp, degree);
        ID.sub(tmp, el, ID.one);
        ID.div(el, tmp, Givaro::Integer(2));
        int k;
        Fx.powmod(b, a, ID.convert(k, el), f);
        Fx.subin(b, Fx.one);
        Fx.mod(a, b, f);
        b = a;
        for (const auto &u : ans) {
            if (Fx.degree(u).value() <= degree) {
                continue;
            }
            Polynomial g;

            // g = gcd(u, b)
            Fx.gcd(g, u, b);

            if (!(Fx.isOne(g)) && (g != u)) {
                sub.insert(u);
                adding.insert(divide_by_lc(Fx, g));

                Polynomial dv;
                Fx.div(dv, u, g);
                adding.insert(divide_by_lc(Fx, dv));
            }
        }

        for (const auto &e : ans) {
            if (sub.find(e) == sub.end()) {
                tmp_set.insert(e);
            }
        }
        ans = tmp_set;
        tmp_set.clear();
        ans.insert(adding.begin(), adding.end());
        adding.clear();
        sub.clear();
    }
    std::vector<Polynomial> final_ans(ans.begin(), ans.end());
    return final_ans;
}

std::vector<std::pair<CantorZassenhaus::Polynomial, int>>
CantorZassenhaus::run(const PolynomialRing &Fx, const Polynomial &f) {
    auto copy = divide_by_lc(Fx, f);
    Polynomial lc;
    Fx.div(lc, f, copy);
    auto distinct = compute_distinct_degree(Fx, square_free(Fx, copy));
    std::vector<std::pair<Polynomial, int>> factorization;
    for (const auto &element : distinct) {
        auto equal_degree_factors = factorize_equal_degree(Fx, element.first, element.second);
        for (const auto &factor : equal_degree_factors) {
            if (Fx.degree(factor).value() == 0) {
                continue;
            }
            int power = 1;
            Polynomial g;
            Fx.gcd(g, copy, factor);
            while (Fx.degree(g).value() != 0) {
                power++;
                Polynomial tmp;
                Fx.div(tmp, copy, factor);
                copy = tmp;
                Fx.gcd(g, copy, factor);
            }
            factorization.emplace_back(factor, power - 1);
        }
    }
    Polynomial tmp;
    Fx.mul(tmp, lc, copy);
    factorization.emplace_back(tmp, 1);
    return factorization;
}

CantorZassenhaus::Polynomial
CantorZassenhaus::divide_by_lc(const CantorZassenhaus::PolynomialRing &Fx, const CantorZassenhaus::Polynomial &f) {
    GaloisField::Element lc_int;
    Polynomial res, lc;
    Fx.leadcoef(lc_int, f);
    Fx.assign(lc, lc_int);
    Fx.div(res, f, lc);
    return res;
}
