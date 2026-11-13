import math
from fractions import Fraction

import numpy as np
import pytest
from flint import fmpz_poly

from zeta_explorer import delta
from zeta_explorer.hecke import (
    basis,
    dimension,
    eigenforms,
    hecke_matrix,
    hurwitz,
    kloosterman,
    trace,
)

SMALL = [3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24]
H = [Fraction(1, 3), Fraction(1, 2), 1, 1, 1, Fraction(4, 3), 2, Fraction(3, 2), 1, 2, 3, 2]
ONE_DIM_A2 = {12: -24, 16: 216, 18: -528, 20: 456, 22: -288, 26: -48}


def test_dimension_is_the_trace_of_t1():
    for k in range(4, 120, 2):
        assert dimension(k) == trace(k, 1)
    assert dimension(2) == dimension(13) == 0


def test_basis_is_unitriangular_and_starts_with_delta():
    g = basis(36, 10)
    for j, f in enumerate(g, start=1):
        assert all(f[i] == 0 for i in range(j)) and f[j] == 1
    tau = delta.tau(40)
    D = basis(12, 40)[0]
    assert [int(D[n]) for n in range(41)] == [int(t) for t in tau[:41]]


def test_one_dimensional_spaces():
    for k, a2 in ONE_DIM_A2.items():
        assert hecke_matrix(k).entries() == [a2]
        assert eigenforms(k)[0].a2 == a2


def test_weight_24_matches_lmfdb():
    assert hecke_matrix(24).charpoly() == fmpz_poly([-20468736, -1080, 1])
    beta = 12 * math.sqrt(144169)
    a2 = [f.a2 for f in eigenforms(24)]
    assert a2 == pytest.approx([540 - beta, 540 + beta], rel=1e-15)
    assert trace(24, 2) == 1080 and trace(24, 3) == 339480


def test_hecke_operators_commute():
    for k in (36, 48, 60):
        T2, T3 = hecke_matrix(k, 2), hecke_matrix(k, 3)
        assert T2 * T3 == T3 * T2


def test_t2_is_irreducible_up_to_weight_120():
    for k in range(24, 122, 2):
        if dimension(k) > 1:
            _, factors = hecke_matrix(k).charpoly().factor()
            assert len(factors) == 1 and factors[0][1] == 1


def test_hurwitz_class_numbers():
    assert [hurwitz(n) for n in SMALL] == H
    assert hurwitz(0) == Fraction(-1, 12)
    assert hurwitz(5) == hurwitz(6) == 0


def test_kronecker_hurwitz_relation():
    for n in range(1, 60):
        t = math.isqrt(4 * n)
        left = sum(hurwitz(4 * n - s * s) for s in range(-t, t + 1))
        divisors = [d for d in range(1, n + 1) if n % d == 0]
        assert left == 2 * sum(divisors) - sum(min(d, n // d) for d in divisors)


def test_trace_formula_gives_tau():
    tau = delta.tau(40)
    assert [trace(12, n) for n in range(1, 41)] == [int(t) for t in tau[1:41]]
    assert trace(14, 5) == 0
    assert trace(36, 1) == dimension(36)


@pytest.mark.parametrize("k", [48, 100])
def test_eigenvalues_against_the_trace_formula(k):
    F = eigenforms(k)
    p = F[0].p
    lam = np.array([f.lam for f in F])
    for i in np.flatnonzero(p < 60):
        q = int(p[i])
        scale = q ** ((k - 1) / 2)
        assert lam[:, i].sum() * scale == pytest.approx(trace(k, q), rel=1e-12)
        assert (lam[:, i] ** 2 - 1).sum() * scale**2 == pytest.approx(trace(k, q * q), rel=1e-12)


def test_deligne_bound_and_certified_radius():
    for k in range(12, 102, 2):
        for f in eigenforms(k):
            assert f.rad < 1e-20
            assert np.all(np.abs(f.lam) < 2)


def test_kloosterman_symmetry_and_weil_bound():
    for c in range(1, 60):
        tau = sum(1 for d in range(1, c + 1) if c % d == 0)
        for m, n in [(1, 1), (1, 2), (3, 5)]:
            s = kloosterman(m, n, c)
            assert s == pytest.approx(kloosterman(n, m, c), abs=1e-9)
            if math.gcd(m * n, c) == 1:
                assert abs(s) <= tau * math.sqrt(c) + 1e-9
    assert kloosterman(1, 1, 3) == pytest.approx(-1)
