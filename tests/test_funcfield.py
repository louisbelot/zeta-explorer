import math

import numpy as np
import pytest

from zeta_explorer.funcfield import (
    Curve,
    Fit,
    LPoly,
    Twist,
    _newton,
    central_k,
    compare,
    family,
    from_power_sums,
    legendre,
    lowest,
    lpoly,
    power_sums,
    twist_set,
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


@pytest.mark.parametrize("n, N", [(1, 1), (2, 4), (3, 5), (4, 8)])
def test_legendre_twist_degree(n, N):
    D = twist_set(5, n, legendre(5).disc)[0]
    assert legendre(5).twist(D).degree == N


def test_family_agrees_with_point_counts_on_each_twist():
    E = legendre(5)
    fam = family(E, 4)
    for t in fam[::40]:
        assert lpoly(E.twist(t.D)) == t.L


@pytest.fixture(scope="module")
def fam5():
    return family(legendre(5), 4)


def test_riemann_hypothesis(fam5):
    for t in fam5:
        z = t.L.roots()
        assert np.all(np.abs(np.abs(z) - 1) < 1e-6)
        assert len(z) + t.L.order(1) + t.L.order(-1) == 8


def test_family_size(fam5):
    q, n = 5, 4
    a = [1, q] + [q**k - q ** (k - 1) for k in range(2, n + 1)]
    count = sum(a[k] * (-1) ** (n - k) * (n - k + 1) for k in range(n + 1))
    assert len(fam5) == count == 345


def test_root_number_is_the_parity_of_the_rank(fam5):
    assert all(t.L.eps == (-1) ** t.rank for t in fam5)
    assert {t.L.eps for t in fam5} == {1, -1}


def test_central_values_are_discrete(fam5):
    for t in fam5:
        if t.L.eps == 1:
            assert t.w >= 0
            assert (t.w == 0) == (t.rank > 0)


def test_so_angles(fam5):
    for t in fam5:
        theta = t.L.so_angles()
        assert len(theta) == (4 if t.L.eps == 1 else 3)
        assert np.all((theta >= 0) & (theta <= np.pi))
        if t.L.eps == 1:
            assert np.prod(2 - 2 * np.cos(theta)) == pytest.approx(float(t.L.central), abs=1e-9)


def test_functional_equation_check():
    with pytest.raises(ArithmeticError):
        from_power_sums([1, 5, 2], 4, 5)


def test_order_counts_exact_factors():
    q = 7
    L = LPoly(q, (1, -2 * q, q * q), 1)
    assert L.order(1) == 2 and L.order(-1) == 0
    assert L.deflated() == [1]


def test_w_needs_integral_value():
    t = Twist((), LPoly(5, (1, 1, 1), 1))
    with pytest.raises(ArithmeticError):
        _ = t.w


def test_rejects_bad_curves():
    with pytest.raises(ValueError):
        Curve(5, (), (), ())
    with pytest.raises(ValueError):
        Curve(9, (), (1,), (0, 1))
    with pytest.raises(ValueError):
        _ = Curve(5, (), (1,), (2,)).degree


def test_central_values_are_four_q_times_a_square(fam5):
    k = central_k(fam5)
    even = [t for t in fam5 if t.L.eps == 1]
    assert len(k) == len(even)
    assert [t.w for t in even] == [4 * 25 * v * v for v in k]
    assert np.sum(k == 0) == sum(t.rank > 0 for t in even)


def test_lowest_angles(fam5):
    even, odd = lowest(fam5, 1), lowest(fam5, -1)
    assert len(even) == sum(t.L.eps == 1 and t.rank == 0 for t in fam5)
    assert len(odd) == sum(t.rank == 1 for t in fam5)
    assert np.all(even > 0) and np.all(odd > 0)


def test_excision_moves_the_model_toward_the_data(fam5):
    fit = compare(fam5, 4 / 25, size=20000)
    assert isinstance(fit, Fit)
    assert (fit.q, fit.N, fit.even + fit.odd) == (5, 8, 311)
    assert fit.ks_excised < fit.ks_so
    assert 0 < fit.rank2 < fit.predicted < 0.5


def test_compare_needs_even_degree():
    with pytest.raises(ValueError):
        compare(family(legendre(5), 1), 0.1)
