
import numpy as np
import pytest

from zeta_explorer.cmtwist import (
    Gauss,
    angles,
    curve,
    frobenius,
    hecke,
    jacobi,
    split,
    units,
)
from zeta_explorer.cubic import index, root, zeta3
from zeta_explorer.funcfield import lpoly
from zeta_explorer.gf import field

CASES = [(7, 6), (13, 6), (5, 4), (13, 4)]


def test_gauss_arithmetic():
    i = Gauss(0, 1)
    assert i * i == Gauss(-1)
    a, b = Gauss(3, -2), Gauss(1, 4)
    assert (a * b).norm == a.norm * b.norm
    assert (a * a.conj()).b == 0 and a.tr == 6
    with pytest.raises(ArithmeticError):
        Gauss(3, 1).div(2)


@pytest.mark.parametrize("p, n", CASES)
def test_root_order(p, n):
    z = root(p, n)
    assert [pow(z, k, p) == 1 for k in range(1, n + 1)] == [False] * (n - 1) + [True]
    if n == 6:
        assert z * z % p == zeta3(p)
    assert units(n)[1] ** n == units(n)[0]


def test_root_needs_congruence():
    with pytest.raises(ValueError):
        root(7, 4)
    with pytest.raises(ValueError):
        root(5, 6)


@pytest.mark.parametrize("p, n", CASES)
def test_index_compatible_with_norm(p, n):
    base = index(field(p, 1), n)
    for m in (2, 3):
        idx = index(field(p, m), n)
        assert all(idx[a] == m * base[a] % n for a in range(1, p))


@pytest.mark.parametrize("p, n", CASES)
def test_hasse_davenport(p, n):
    pi = frobenius(p, n)
    assert pi.norm == p
    for m in (1, 2, 3):
        assert jacobi(p, n, m).norm == p**m
        assert frobenius(p, n, m) == pi**m


@pytest.mark.parametrize("p, n", CASES)
def test_frobenius_trace(p, n):
    x = np.arange(p)
    rhs = (x**3 + 1) if n == 6 else (x**3 + x)
    legendre = np.array([0] + [1 if pow(int(v), (p - 1) // 2, p) == 1 else -1 for v in range(1, p)])
    assert frobenius(p, n).tr == -int(legendre[rhs % p].sum())


@pytest.mark.parametrize("p, n, c", [
    (7, 6, (1, 0, 1)), (7, 6, (3, 2, 0, 1)), (7, 6, (0, 0, 1, 4, 1)), (7, 6, (2, 1, 1, 1, 1, 1, 1)),
    (13, 6, (5, 1, 1)), (5, 4, (1, 0, 1)), (5, 4, (0, 2, 3, 1)), (5, 4, (1, 1, 0, 0, 1)),
    (13, 4, (2, 7, 1)), (17, 4, (3, 0, 1, 1)),
])
def test_hecke_split(p, n, c):
    H = hecke(p, n, [c])[0]
    assert split(H, frobenius(p, n)) == list(lpoly(curve(n, p, c)).c)


def test_hecke_angles_on_circle():
    p, n = 7, 6
    H = hecke(p, n, [(3, 2, 0, 1)])[0]
    th = angles(H, frobenius(p, n), p)
    assert len(th) == len(H) - 1 and np.all(np.abs(th) <= np.pi)
