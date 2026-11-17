import math

import numpy as np
import pytest

from zeta_explorer.gf import factor, monics, squarefree
from zeta_explorer.hyperelliptic import andrade_keating, irreducible_counts, lpolys, usp_moment
from zeta_explorer.hypermoments import (
    _symbol,
    arithmetic_factor,
    dual_coefficients,
    dual_mean,
    florea_D,
    florea_R,
    irreducibles,
    lpoly,
    recipe,
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


def test_florea_product_closed_form():
    q = 5
    z0 = q ** (-4 / 3)
    pi = irreducible_counts(q, 80)
    prod = sum(n * math.log1p(-(q ** (4 * d / 3) + q ** (2 * d / 3) + q ** (d / 3) + 1)
                              / (q ** (4 * d / 3) + q**d) ** 2) for d, n in enumerate(pi) if d)
    assert florea_D(q, z0, q * z0).real == pytest.approx(prod, abs=1e-12)


def test_florea_R_is_linear():
    r = [florea_R(13, x) for x in (5, 7, 9, 11)]
    assert np.diff(r, 2) == pytest.approx([0, 0], abs=1e-14)
    with pytest.raises(ValueError):
        florea_R(3, 5)


def test_recipe_first_moment_is_andrade_keating():
    for q in (3, 5):
        gs = [2, 5, 10]
        for got, g in zip(recipe(q, gs, 1), gs, strict=True):
            assert got == pytest.approx(andrade_keating(q, g), rel=1e-12)


@pytest.mark.parametrize("k", [2, 3])
def test_recipe_is_a_polynomial_in_g(k):
    m = k * (k + 1) // 2
    vals = recipe(5, list(range(2, m + 4)), k)
    assert abs(np.diff(vals, m + 1)).max() < 1e-9 * max(vals)
    assert abs(np.diff(vals, m)).min() > 1e-3


def test_recipe_leading_term_is_keating_snaith():
    a = arithmetic_factor(3, 2)
    gaps = [abs(r / (a * usp_moment(g, 2)) - 1) for r, g in zip(recipe(3, [200, 800], 2),
                                                             [200, 800], strict=True)]
    assert gaps[1] < 0.01 and gaps[0] / gaps[1] == pytest.approx(4, rel=0.05)


def test_second_moment_against_recipe():
    gaps = []
    for g in (3, 4):
        c = lpolys(3, g).central
        gaps.append(abs(np.mean(c**2) / recipe(3, [g], 2)[0] - 1))
    assert gaps[1] < gaps[0] < 2e-3
