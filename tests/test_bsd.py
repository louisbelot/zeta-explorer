import math

import pytest

from zeta_explorer.bsd import (
    _torsion_images,
    fibers,
    klass,
    legendre_roots,
    local,
    local_image,
    place_degree,
)
from zeta_explorer.funcfield import family, legendre
from zeta_explorer.gf import factor, monics, squarefree


@pytest.fixture(scope="module")
def fam5():
    return family(legendre(5), 4)


def _is_square(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


def test_local_valuation_and_residue():
    p = 5
    assert local((0, 0, 1, 1), (0, 1), p) == (2, 1)
    assert local((0, 2), (0, 1), p) == (1, -1)
    assert local((3, 0, 0, 2), None, p) == (-3, -1)
    assert klass((0, 2), (0, 1), p) == (1, 1)
    with pytest.raises(ZeroDivisionError):
        local((), (0, 1), p)


@pytest.mark.parametrize("p,d", [(5, 2), (5, 3), (7, 2), (3, 4)])
def test_residue_character_is_the_norm_character(p, d):
    for P in monics(p, d):
        if not squarefree(P, p) or len(factor(P, p)) != 1:
            continue
        norm = (-1) ** d * P[0] % p
        assert local((0, 1), P, p)[1] == (1 if pow(norm, (p - 1) // 2, p) == 1 else -1)


def test_legendre_fibers():
    f = fibers(legendre_roots((1,), 7), 7)
    assert {P: x.kind for P, x in f.items()} == {(0, 1): "I2", (6, 1): "I2", None: "I2*"}
    assert [x.c for x in f.values()] == [2, 2, 4]
    assert sum(x.disc for x in f.values()) == 12
    assert sum(x.cond for x in f.values()) == 4


@pytest.mark.parametrize("p,n", [(5, 3), (7, 2), (3, 5)])
def test_conductor_matches_the_curve(p, n):
    E = legendre(p)
    for t in family(E, n)[::7]:
        f = fibers(legendre_roots(t.D, p), p)
        cond = sum(x.cond * place_degree(P) for P, x in f.items())
        twisted = E.twist(t.D)
        assert cond == twisted.f_finite + twisted.f_inf == t.L.N + 4
        assert sum(x.disc * place_degree(P) for P, x in f.items()) % 12 == 0


def test_local_images_have_order_four(fam5):
    t = fam5[17]
    e = legendre_roots(t.D, 5)
    for P in fibers(e, 5):
        W = local_image(e, P, 5)
        assert len(W) == 2
        span = {tuple((x * a + y * b) % 2 for a, b in zip(*W, strict=True))
                for x in (0, 1) for y in (0, 1)}
        assert all(v in span for v in _torsion_images(e, P, 5))
