from collections import Counter
from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.funcfield import family, legendre, traces, twist_set
from zeta_explorer.gf import factor, field, monics, pdivmod, pgcd, squarefree
from zeta_explorer.twistmoments import (
    coefficients,
    diagonal,
    first_moment,
    place_traces,
    recipe_first,
    so_ratios,
    squarefree_counts,
)


@pytest.mark.parametrize("q,n", [(3, 4), (5, 3), (7, 2)])
def test_squarefree_counts(q, n):
    assert squarefree_counts(q, n)[n] == len(twist_set(q, n, legendre(q).disc))


def test_local_coefficients():
    assert coefficients(5, 6, 6) == [3 ** (m + 1) - 2 ** (m + 1) for m in range(7)]


@pytest.mark.parametrize("q,d", [(3, 1), (3, 3), (5, 2), (7, 1)])
def test_place_traces(q, d):
    a = place_traces(q, d)
    assert len(a) == len(_irreducible(q, d)) - 2 * (d == 1)
    assert np.all(a**2 <= 4 * q**d)
    assert np.all((q**d + 1 - a) % 4 == 0)


def _irreducible(q, d):
    return [f for f in monics(q, d) if squarefree(f, q) and len(factor(f, q)) == 1]


def _trace(P, q):
    F = field(q, len(P) - 1)
    x = np.arange(F.q)
    root = x[F.polyval(P, x) == 0][0]
    return int(traces(legendre(q), len(P) - 1)[root])


def _factorize(f, irreducibles, q):
    out = Counter()
    for P in irreducibles:
        while len(f) >= len(P):
            quo, rem = pdivmod(f, P, q)
            if rem:
                break
            f, out[P] = quo, out[P] + 1
    return out


@pytest.mark.parametrize("q,n", [(3, 4), (5, 2)])
def test_diagonal_against_brute_force(q, n):
    Ds = twist_set(q, n, legendre(q).disc)
    irr = [P for d in range(1, n + 1) for P in _irreducible(q, d)]
    M = {(0, 1): 0, (q - 1, 1): 1}
    a0, a1 = (int(x) for x in traces(legendre(q), 1)[[0, 1]])
    want = [0] * (n + 1)
    for j in range(n + 1):
        for f in monics(q, j):
            e = _factorize(f, irr, q)
            odd = {P for P, k in e.items() if k % 2}
            if odd and odd != set(M):
                continue
            a = a0 ** e[(0, 1)] * a1 ** e[(q - 1, 1)]
            for P, k in e.items():
                if P not in M:
                    a *= coefficients(_trace(P, q), q ** (len(P) - 1), k)[k]
            want[j] += a * sum(len(pgcd(D, f, q)) == 1 for D in Ds)
    assert [int(c) for c in diagonal(q, n)] == want


@pytest.mark.parametrize("q", [5, 7])
def test_root_number_is_the_symbol_of_t_and_t_minus_1(q):
    chi = field(q, 1).chi
    for t in family(legendre(q), 4):
        assert t.L.eps == chi[t.D[0]] * chi[sum(t.D) % q]


def test_recipe_converges_in_the_degree():
    gaps = [abs(recipe_first(3, n) - first_moment(family(legendre(3), n))) for n in (2, 4, 6)]
    assert gaps[0] > 2 * gaps[1] > 4 * gaps[2]
    assert gaps[2] < Fraction(1, 15)
    with pytest.raises(ValueError):
        recipe_first(3, 3)


def test_recipe_small_case_by_hand():
    fam = family(legendre(3), 2)
    assert [t.D for t in fam] == [(1, 0, 1), (2, 1, 1), (2, 2, 1)]
    assert first_moment(fam) == Fraction(4, 9)
    assert recipe_first(3, 2) == Fraction(8, 9)


def test_large_q_moments_approach_so4():
    r = so_ratios(family(legendre(13), 2), [1, 2])
    assert r[0] == pytest.approx(1, abs=0.04)
    assert r[1] == pytest.approx(1, abs=0.15)
