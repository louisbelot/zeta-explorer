import math

import pytest

from zeta_explorer.gf import factor, monics, squarefree
from zeta_explorer.hyperelliptic import irreducible_counts, lpolys
from zeta_explorer.hypermoments import (
    _symbol,
    dual_coefficients,
    dual_mean,
    irreducibles,
    lpoly,
    squarefrees,
    squares,
)


def test_irreducibles():
    for q in (3, 5):
        for d in range(1, 5):
            Ps = irreducibles(q, d)
            assert len(Ps) == irreducible_counts(q, d)[d]
            assert all(P[-1] == 1 and len(P) == d + 1 for P in Ps)
            assert all(len(factor(P, q)) == 1 for P in Ps[:30])


def test_squarefrees_with_factor_degrees():
    q, n = 3, 5
    got = {f: degs for f, degs in squarefrees(q, n)}
    want = [f for d in range(n + 1) for f in monics(q, d) if squarefree(f, q)]
    assert set(got) == set(want) and len(got) == len(want)
    for f in want[1:200]:
        assert sorted(got[f]) == sorted(len(P) - 1 for P in factor(f, q))


def test_every_monic_is_squarefree_times_square():
    q = 3
    counts = [0] * 7
    hs = squares(q, 6)
    for f, _ in squarefrees(q, 6):
        for h, _ in hs:
            j = len(f) - 1 + 2 * (len(h) - 1)
            if j <= 6:
                counts[j] += 1
    assert counts == [q**j for j in range(7)]


def _brute(q, f, k):
    Ps = factor(f, q)
    return sum(math.prod(_symbol(D, P, q) for P in Ps) for D in monics(q, k))


@pytest.mark.parametrize("q", [3, 5])
def test_lpoly_against_dirichlet_series(q):
    for d in range(1, 6):
        fs = [f for f in monics(q, d) if squarefree(f, q)][:6]
        for f, L in zip(fs, lpoly(q, fs), strict=True):
            assert list(L) == [_brute(q, f, k) for k in range(d)]


@pytest.mark.parametrize("q, g", [(3, 2), (3, 3), (5, 2), (7, 2)])
def test_dual_coefficients_match_brute_force(q, g):
    E = lpolys(q, g)
    assert dual_coefficients(q, g) == [int(x) for x in E.c[:, : g + 1].sum(axis=0)]


def test_dual_mean_matches_every_curve():
    for q, g in [(3, 4), (5, 3)]:
        assert dual_mean(q, g) == pytest.approx(lpolys(q, g).central.mean(), rel=1e-13)
