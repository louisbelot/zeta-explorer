import cmath
import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.elliptic import Curve, EllipticL, central, completed, exp1, twists
from zeta_explorer.explicit import primes

E11 = Curve(0, -1, 1, 0, 0)
E37 = Curve(0, 0, 1, -1, 0)
E389 = Curve(0, 1, 1, -2, 0)
E5077 = Curve(0, 0, 1, -7, 6)
KAPPA = 6.346046521


def eta_product(M):
    f = np.zeros(M + 1, dtype=object)
    f[1] = 1
    for n in range(1, M + 1):
        for k in (n, n, 11 * n, 11 * n):
            if k <= M:
                f[k:] = f[k:] - f[: M + 1 - k]
    return f


def test_invariants():
    assert E11.discriminant == -11
    assert Curve(0, -1, 1, -10, -20).discriminant == -(11**5)
    assert [E.conductor for E in (E11, E37, E389, E5077)] == [11, 37, 389, 5077]


def test_additive_reduction_is_rejected():
    with pytest.raises(ValueError, match="additive"):
        EllipticL(Curve(0, 0, 0, -1, 0))


def test_coefficients_match_eta_product():
    assert EllipticL(E11).an(300).tolist() == eta_product(300).tolist()


def test_isogenous_curves_share_coefficients():
    a = EllipticL(Curve(0, -1, 1, -10, -20)).an(500)
    assert a.tolist() == EllipticL(E11).an(500).tolist()


def test_37a_q_expansion():
    assert EllipticL(E37).an(11)[1:].tolist() == [1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5]


def test_hasse_and_hecke():
    a = EllipticL(E5077).an(3000)
    for p in primes(50).tolist():
        if p != 5077:
            assert a[p * p] == a[p] ** 2 - p
    ps = primes(3000)
    assert np.all(a[ps] ** 2 <= 4 * ps)
    for m, n in [(4, 9), (7, 25), (12, 35), (16, 81)]:
        assert a[m * n] == a[m] * a[n]


def test_root_numbers_follow_rank_parity():
    assert [EllipticL(E).root_number for E in (E11, E37, E389, E5077)] == [1, -1, 1, -1]


def test_twist_data():
    L = EllipticL(E11, -4)
    assert L.conductor == 11 * 16
    assert L.root_number == 1
    assert EllipticL(E11, 5).root_number == 1
    assert EllipticL(E11, 8).root_number == -1
    chi = [0, 1, 0, -1, 0]
    assert (L.an(40)[1:] == EllipticL(E11).an(40)[1:] * [chi[n % 4] for n in range(1, 41)]).all()
    with pytest.raises(ValueError):
        EllipticL(E11, 33)
    with pytest.raises(ValueError):
        EllipticL(E11, 9)


def test_exp1_against_mpmath():
    x = np.array([1e-3, 0.3, 0.999, 1.0, 2.5, 30.0])
    ref = [float(mpmath.e1(v)) for v in x]
    assert exp1(x) == pytest.approx(ref, rel=1e-14)


def test_central_values():
    assert central(EllipticL(E11)) == pytest.approx(0.2538418608559107, abs=1e-14)
    assert central(EllipticL(E37)) == pytest.approx(0.3059997738340523, abs=1e-14)
    assert abs(central(EllipticL(E389))) < 1e-13
    assert abs(central(EllipticL(E5077))) < 1e-13


@pytest.mark.parametrize("L", [EllipticL(E11), EllipticL(E37), EllipticL(E11, -7),
                               EllipticL(E11, 8)], ids=repr)
def test_split_point_does_not_matter(L):
    s = np.array([0.5, 0.3 + 0.7j, 0.8 + 2.5j])
    ref = completed(s, L)
    for y0 in (1.3, 0.8, cmath.exp(0.6j)):
        assert np.abs(completed(s, L, y0) - ref).max() < 1e-12 * np.abs(ref).max()


def test_wrong_root_number_breaks_the_split():
    L = EllipticL(E11, 8)
    L.__dict__["root_number"] = 1
    s = 0.3 + 0.7j
    assert abs(completed(s, L, 1.3) - completed(s, L)) > 1e-3 * abs(completed(s, L))


def test_waldspurger_squares():
    even = [L for L in twists(E11, 5, 400, sign=1) if L.root_number == 1]
    c = np.array([math.sqrt(max(central(L), 0) * math.sqrt(L.d) / KAPPA) for L in even])
    assert len(even) == 52
    assert np.abs(c - c.round()).max() < 1e-7
    assert c[:7].round().tolist() == [1, 1, 1, 2, 2, 1, 3]
    assert [L.d for L, ci in zip(even, c, strict=True) if ci < 0.5] == [232, 265, 273, 344, 364]


def test_twists_family():
    fam = twists(E11, 1, 40, sign=-1)
    ds = [L.d for L in fam]
    assert ds == [-3, -4, -7, -8, -15, -19, -20, -23, -24, -31, -35, -39, -40]
    assert all(L.root_number in (1, -1) for L in fam)
