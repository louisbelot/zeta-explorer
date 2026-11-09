import math

import numpy as np
import pytest

from zeta_explorer import delta as D
from zeta_explorer.explicit import primes

LMFDB_TAU = [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920, 534612,
             -370944, -577738, 401856, 1217160, 987136, -6905934]



def test_tau_lmfdb():
    assert list(D.tau(17)[1:]) == LMFDB_TAU


def test_tau_against_flint_product():
    flint = pytest.importorskip("flint")
    N = 1500
    eta = flint.fmpz_poly([1])
    for n in range(1, N):
        eta = (eta * flint.fmpz_poly([1] + [0] * (n - 1) + [-1])).truncate(N)
    ref = [int(c) for c in eta.pow_trunc(24, N).coeffs()]
    assert list(D.tau(N)[1:]) == ref


def test_ramanujan_691():
    N = 3000
    t = D.tau(N)
    sigma = [0] * (N + 1)
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            sigma[m] += pow(d, 11, 691)
    assert all((t[n] - sigma[n]) % 691 == 0 for n in range(1, N + 1))


def test_hecke_relations():
    t = D.tau(2000)
    for m in range(2, 40):
        for n in range(2, 2000 // m + 1):
            if math.gcd(m, n) == 1:
                assert t[m * n] == t[m] * t[n]
    for p in primes(44):
        assert t[p * p] == t[p] ** 2 - p**11


def test_deligne_bound_exact():
    t = D.tau(1 << 15)
    for p in primes(1 << 15):
        assert t[p] ** 2 < 4 * int(p) ** 11


def test_moduli_cover_the_bound():
    for N in (1024, 1 << 20):
        assert math.prod(D._moduli(N)) > 4 * N**6


def test_lam_is_normalized_tau():
    n = np.arange(1, 101)
    assert np.allclose(D.lam(100)[1:] * n**5.5, [float(x) for x in D.tau(100)[1:]], rtol=1e-15)
