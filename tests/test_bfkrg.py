import itertools
import math

import numpy as np
import pytest

from zeta_explorer.bfkrg import (
    check_curve,
    chi_M,
    h2,
    local,
    raw2,
    root_number,
    sym2,
)
from zeta_explorer.funcfield import Curve, family, legendre, twisted_legendre
from zeta_explorer.gf import monics, pdivmod
from zeta_explorer.onelevel import GOOD, MULT, multiplicative, places

E5 = twisted_legendre(5)


@pytest.fixture(scope="module")
def fams():
    return {n: family(E5, n) for n in (1, 3, 5)}


@pytest.fixture(scope="module")
def primes():
    out = {}
    for d in (1, 2, 3):
        L = local(E5, d)
        for P, lam, kind in zip(places(E5, d).P, L.lam, L.kind, strict=True):
            out[P] = (lam, kind)
    return out


def _lam(lam, kind, k):
    if kind == GOOD:
        a, b = 1.0, lam
        for _ in range(k):
            a, b = b, lam * b - a
        return a
    return lam**k if kind == MULT else float(k == 0)


def _data(f, primes):
    lam, odd, seen = 1.0, set(), set()
    for P, (lp, kind) in primes.items():
        k = 0
        while len(f) > 1:
            quo, rem = pdivmod(f, P, 5)
            if rem:
                break
            f, k = quo, k + 1
        if k:
            lam *= _lam(lp, kind, k)
            seen.add(P)
            if k % 2:
                odd.add(P)
    return lam, frozenset(odd), seen


def _brute(a, b, twisted, primes):
    bad = {(0, 1), (4, 1)}
    M = frozenset(multiplicative(E5)) if twisted else frozenset()
    F = [_data(f, primes) for f in monics(5, a)]
    H = [_data(h, primes) for h in monics(5, b)]
    total = 0.0
    for (lf, of, pf), (lh, oh, ph) in itertools.product(F, H):
        if of ^ oh ^ M:
            continue
        w = math.prod(1 / (1 + 5.0 ** -(len(P) - 1)) for P in pf | ph if P not in bad)
        total += lf * lh * w
    return total / 5 ** ((a + b) / 2)


@pytest.mark.parametrize("twisted", [False, True])
def test_euler_product_coefficients(twisted, primes):
    n, r = 16, 0.2
    z = r * np.exp(2j * np.pi * np.arange(n) / n)
    V = np.array([[raw2(E5, x, y, 3, twisted) for y in z] for x in z])
    C = (np.fft.fft2(V) / n**2).real / r ** np.add.outer(np.arange(n), np.arange(n))
    for a, b in [(0, 0), (1, 1), (2, 0), (0, 3), (2, 1), (1, 2), (3, 1), (2, 2)]:
        assert C[a, b] == pytest.approx(_brute(a, b, twisted, primes), abs=1e-9)


def test_sym2_factors_reassemble_the_raw_product():
    for twisted in (False, True):
        raw = raw2(E5, 0.3, 0.2, 7, twisted) * (1 - 0.3 * 0.2)
        assert h2(E5, 0.3, 0.2, 7, twisted) == pytest.approx(raw, rel=1e-9)


def test_sym2_of_the_twist_is_trivial():
    assert sym2(E5) == (1,)


def test_curve_checks():
    check_curve(E5)
    check_curve(legendre(5))
    with pytest.raises(ValueError):
        check_curve(Curve(5, (), (), (0, 1)))


def test_chi_M_is_the_residue_of_D_at_one():
    M = multiplicative(E5)
    assert M == [(4, 1)]
    for D in [(1, 1), (2, 0, 1), (1, 1, 0, 1)]:
        v = sum(D) % 5
        assert chi_M(D, M, 5) == (1 if pow(v, 2, 5) == 1 else -1)


def test_twist_degree(fams):
    assert all(t.L.N == 2 * n + 1 for n, fam in fams.items() for t in fam)


def test_root_number_is_minus_chi_at_one(fams):
    M = multiplicative(E5)
    for fam in fams.values():
        assert root_number(fam, M) == -1
        odd = [t for t in fam if t.L.eps == -1]
        assert odd and all(t.L.central == 0 for t in odd)
