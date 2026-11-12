import math

import pytest

from zeta_explorer.funcfield import (
    Curve,
    LPoly,
    _newton,
    from_power_sums,
    legendre,
    lpoly,
    power_sums,
    ulmer,
)


def ulmer_rank(p, d):
    q = p
    r = sum(_phi(e) // _order(q, e) for e in range(1, d + 1) if d % e == 0 and 6 % e)
    r += 1 if d % 2 == 0 and (q - 1) % 4 == 0 else 0
    if d % 3 == 0:
        r += 2 if (q - 1) % 3 == 0 else 1
    return r


def _phi(e):
    return sum(1 for k in range(1, e + 1) if math.gcd(k, e) == 1)


def _order(q, e):
    k, x = 1, q % e
    while x != 1:
        x, k = x * q % e, k + 1
    return k


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_legendre_curve_has_trivial_l(p):
    E = legendre(p)
    assert E.degree == 0
    assert power_sums(E, 3) == [0, 0, 0]
    assert lpoly(E).c == (1,)


@pytest.mark.parametrize("p, d", [(5, 2), (5, 3), (5, 6), (7, 2), (7, 4), (7, 8), (11, 3),
                                  (11, 4), (11, 6)])
def test_ulmer_ranks(p, d):
    E = ulmer(p, d)
    assert E.degree == (d - 3 if d % 6 == 0 else d - 1)
    L = lpoly(E)
    assert L.order(1) == ulmer_rank(p, d)
    assert L.eps == (-1) ** L.order(1)


def test_l_has_exactly_the_predicted_degree():
    E = legendre(5).twist((2, 1, 1))
    N = E.degree
    assert N == 4
    c = _newton(power_sums(E, N + 2))
    assert c[N + 1:] == [0, 0]
    assert tuple(c[:N + 1]) == lpoly(E).c


def test_functional_equation_check():
    with pytest.raises(ArithmeticError):
        from_power_sums([1, 5, 2], 4, 5)


def test_order_counts_exact_factors():
    q = 7
    L = LPoly(q, (1, -2 * q, q * q), 1)
    assert L.order(1) == 2 and L.order(-1) == 0
    assert L.deflated() == [1]


def test_rejects_bad_curves():
    with pytest.raises(ValueError):
        Curve(5, (), (), ())
    with pytest.raises(ValueError):
        Curve(9, (), (1,), (0, 1))
    with pytest.raises(ValueError):
        _ = Curve(5, (), (1,), (2,)).degree
