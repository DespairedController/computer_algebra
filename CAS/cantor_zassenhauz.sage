def __divide(f, g):
   return f.quo_rem(g)[0]


def __rooting(f):
    p = f.base_ring().characteristic()
    q = f.base_ring().order()
    variable = f.variables()[0]
    coeffs = f.list()
    g = 0
    i = 0
    for coef in coeffs:
        if (coef != 0):
            g += coef^(q // p) * variable ^ (i // p)
        i += 1
    return g


def __square_free(f):
    m = 1
    ans = 1
    p = f.base_ring().characteristic()
    
    d = f.gcd(f.derivative())
    w = __divide(f, d)
    
    while w != 1:
        y = w.gcd(d)
        factor = __divide(w, y)
        ans *= factor
        w = y
        d = __divide(d, y)
    
    if d != 1:
        d = __rooting(d)
        ans *= __square_free(d)
    return ans


def __compute_distinct_degree(f):
    R = f.base_ring()
    q = R.order()
    n = f.degree()
    variable = f.variables()[0]
    
    ans = []
    i = 0
    v = variable
    while i + 1 <= f.degree() // 2:
        i += 1
        v = (v ^ q).quo_rem(f)[1]
        g = f.gcd(v - variable)
        g = g / g.lc()
        if g != 1:
            ans.append((g, i))
            f = __divide(f, g)
            v = v.quo_rem(f)[1]
    if f != 1:
        ans.append((f/f.lc(), f.degree()))
    return ans


def __find_divisor(f, s):
    R = f.base_ring()
    q = R.order()
    n = f.degree()
    while (true):
        a = f.parent().random_element((1, n))
        a = a / a.lc()
        if not gcd(f, a).is_constant():
            return a
        b = a^((q^s - 1) // 2).quo_rem(f)[1] - 1
        a2 = gcd(b, f)
        if (not a2.is_constant()) and a2 != f:
            return a2


def __factorize_equal_degree(f, s):
    R = f.base_ring()
    q = R.order()
    n = f.degree()
    r = n // s
    ans = {f}
    while len(ans) < r:
        a = f.parent().random_element((1, n))
        b = (a^((q^s - 1) // 2) - 1).quo_rem(f)[1]
        for u in ans:
            if u.degree() == s:
                continue
            g = u.gcd(b)
            adding = set()
            sub = set()
            if g != 1 and g != u:
                sub.add(u)
                adding.add(g)
                adding.add(__divide(u, g))
        ans = ans.difference(sub)
        ans = ans.union(adding)
    return list(ans)


def cantor_zassenhauz(f):
    lc = f.lc()
    f /= lc
    distinct = __compute_distinct_degree(__square_free(f))
    square_free_factorization = []
    for (factor, power) in distinct:
        square_free_factorization += __factorize_equal_degree(factor, power)
    factorization = []
    for factor in square_free_factorization:
        if factor.is_constant():
            continue
        power = 1
        while not gcd(f, factor).is_constant():
            power += 1
            f = __divide(f, factor)
        factorization += [(factor, power - 1)]
    return [(f * lc, 1)] + factorization


def __test(f):
    g = cantor_zassenhauz(f)
    restored_f = 1
    for (factor, power) in g:
        restored_f *= (factor^power)
    if restored_f != f:
        print("Test failed.\nExpected: " + str(f) + "\nActual: " + str(restored_f))
        return false
    return true


def run_random_tests():
    testing_orders = [3, 5, 7, 9]
    for order in testing_orders:
        for degree in range(3, 5):
            R = PolynomialRing(GF(order), 'x')
            for _ in range(1, 10):
                f = R.random_element((degree, degree))
                print("Testing polynomial " + str(f) + " over GF(" + str(order) +")")
                if not __test(f):
                    return
    print("OK")
    return

