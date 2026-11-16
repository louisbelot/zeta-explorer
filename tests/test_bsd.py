import math

import numpy as np
import pytest

from zeta_explorer.bsd import (
    _torsion_images,
    _v2,
    arith,
    fibers,
    klass,
    legendre_roots,
    local,
    local_image,
    place_degree,
    point_counts,
    selmer_rank,
    v2k,
)
from zeta_explorer.funcfield import family, legendre
from zeta_explorer.gf import factor, monics, peval, squarefree


@pytest.fixture(scope="module")
def fam5():
    return family(legendre(5), 4)


@pytest.fixture(scope="module")
def info5(fam5):
    return [(t, arith(t)) for t in fam5 if t.L.eps == 1]


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


def test_tamagawa_product(info5):
    for _, A in info5:
        assert A.tamagawa == 2 * 2 * 4 * 4**A.omega
        assert A.height == 3


def test_sha_is_a_square(info5):
    for t, A in info5:
        assert A.torsion == 4
        if t.rank == 0:
            assert A.sha.denominator == 1 and A.sha > 0
            assert _is_square(A.sha.numerator)
        else:
            assert A.sha is None


def test_sha_squares_in_characteristic_three():
    for t in family(legendre(3), 6):
        if t.L.eps == 1 and t.rank == 0:
            sha = arith(t).sha
            assert sha.denominator == 1 and _is_square(sha.numerator)


def test_nonsplit_i2_has_two_components(info5):
    split = set()
    for t, A in info5:
        fiber = A.fibers[(0, 1)]
        assert (fiber.kind, fiber.c) == ("I2", 2)
        if t.rank == 0:
            split.add(pow(-t.D[0] % 5, 2, 5) == 1)
    assert split == {True, False}


def test_selmer_against_the_l_function(info5):
    for t, A in info5:
        assert A.selmer >= A.rank + 2
        assert A.sha2 % 2 == 0
        if t.rank == 0:
            n = A.sha.numerator
            assert n % 2**A.sha2 == 0
            assert (A.sha2 == 0) == (n % 2 == 1)
    assert {A.selmer for t, A in info5 if t.rank == 2} == {4}


def test_floor_bounds_the_two_adic_part_of_k(info5):
    gaps = []
    for t, A in info5:
        if t.rank == 0:
            gap = v2k(t) - A.floor2
            assert gap == (_v2(A.sha) - A.sha2) // 2
            gaps.append(gap)
    assert min(gaps) == 0 and np.mean(np.array(gaps) == 0) > 0.85


def test_local_images_have_order_four(fam5):
    t = fam5[17]
    e = legendre_roots(t.D, 5)
    for P in fibers(e, 5):
        W = local_image(e, P, 5)
        assert len(W) == 2
        span = {tuple((x * a + y * b) % 2 for a, b in zip(*W, strict=True))
                for x in (0, 1) for y in (0, 1)}
        assert all(v in span for v in _torsion_images(e, P, 5))
    assert selmer_rank(e, 5) == arith(t).selmer


def test_point_counts():
    D = (1, 0, 1)
    for m in (1, 2, 3):
        c = point_counts(D, 5, m)
        assert np.all(np.abs(c - 5**m - 1) <= 2 * 5 ** (m / 2))
        assert np.all(c % 4 == 0)
    good = [t for t in range(2, 5) if peval(D, t, 5)]
    assert len(point_counts(D, 5, 1)) == len(good)
