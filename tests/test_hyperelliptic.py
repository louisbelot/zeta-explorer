import numpy as np
import pytest

from zeta_explorer.gf import factor, field, monics, peval, squarefree
from zeta_explorer.hyperelliptic import (
    andrade_keating,
    ensemble,
    irreducible_counts,
    lpolys,
    recipe_first,
    usp_moment,
)
from zeta_explorer.rmt import char_poly_at_one, o_minus_angles


@pytest.fixture(scope="module")
def E32():
    return lpolys(3, 2)


def test_ensemble_size():
    for q, g in [(3, 1), (3, 2), (5, 1)]:
        assert len(ensemble(q, g)) == (q - 1) * q ** (2 * g)


def test_irreducible_counts():
    for q in (3, 5):
        want = [0] + [sum(squarefree(f, q) and len(factor(f, q)) == 1 for f in monics(q, d))
                      for d in range(1, 5)]
        assert irreducible_counts(q, 4) == want


def _symbol(D, f, q):
    roots = [a for a in range(q) if peval(f, a, q) == 0]
    chi = field(q, 1).chi
    if not roots:
        F = field(q, 2)
        x = np.arange(F.q)
        theta = x[F.polyval(f, x) == 0][0]
        return int(F.chi[F.polyval(D, theta)])
    if len(f) == 2:
        return int(chi[peval(D, roots[0], q)])
    b = roots[1] if len(roots) == 2 else roots[0]
    return int(chi[peval(D, roots[0], q)]) * int(chi[peval(D, b, q)])


def test_lpolys_against_the_dirichlet_series(E32):
    for D, c in list(zip(E32.D, E32.c, strict=True))[::17]:
        assert c[1] == sum(_symbol(D, (-a % 3, 1), 3) for a in range(3))
        assert c[2] == sum(_symbol(D, f, 3) for f in monics(3, 2))


def test_riemann_hypothesis_and_symmetry(E32):
    q, g = E32.q, E32.g
    for c in E32.c[::23]:
        full = [*c, *(c[j] * q ** (g - j) for j in range(g - 1, -1, -1))]
        z = np.roots(full[::-1])
        assert np.allclose(np.abs(z), q**-0.5, rtol=1e-9)
    assert E32.central.min() > -1e-12
    assert E32.angles().shape == (len(E32.D), g)


def test_first_coefficient_averages_to_zero():
    for q, g in [(3, 2), (5, 1)]:
        assert lpolys(q, g).c[:, 1].sum() == 0
    assert np.mean(lpolys(5, 1).central) == pytest.approx(2, abs=1e-12)


def test_recipe_by_hand(E32):
    H = len(E32.D)
    coprime = sum(peval(D, a, 3) != 0 for D in E32.D for a in range(3))
    assert recipe_first(3, 2) == pytest.approx(2 + coprime / 3 / H, rel=1e-13)


def test_andrade_keating_first_moment():
    E = lpolys(7, 2)
    assert np.mean(E.central) == pytest.approx(andrade_keating(7, 2), abs=1e-4)
    assert abs(np.mean(lpolys(3, 3).central) - andrade_keating(3, 3)) < 0.005


def test_usp_moments():
    assert usp_moment(1, 1) == pytest.approx(2)
    assert usp_moment(1, 2) == pytest.approx(5)
    v = char_poly_at_one(o_minus_angles(3, 40_000, seed=2))
    assert np.mean(v) == pytest.approx(usp_moment(2, 1), rel=0.02)
    assert np.mean(v**2) == pytest.approx(usp_moment(2, 2), rel=0.05)
