import itertools
from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.bsd import _rank2, selmer_rank
from zeta_explorer.funcfield import family, legendre
from zeta_explorer.gf import norm, pgcd, pmul
from zeta_explorer.heights import (
    MAX_CONTR,
    Point,
    _kernel,
    delta,
    lstar,
    model,
    points,
    psqrt,
    rank_one,
    sqclass,
    torsion_classes,
)


@pytest.fixture(scope="module")
def solved():
    out = []
    for q, n in ((5, 2), (3, 3)):
        for t in family(legendre(q), n):
            if t.rank == 1:
                out.append((t, rank_one(t)))
    return out


def test_psqrt():
    rng = np.random.default_rng(3)
    for _ in range(200):
        s = norm([int(c) for c in rng.integers(0, 5, rng.integers(1, 8))], 5)
        if not s:
            continue
        r = psqrt(pmul(s, s, 5), 5)
        assert pmul(r, r, 5) == pmul(s, s, 5)
        assert psqrt(pmul((0, 1), pmul(s, s, 5), 5), 5) is None
    assert psqrt((2,), 5) is None
    assert psqrt((4,), 5) in ((2,), (3,))


def test_kernel_is_the_nullspace():
    rng = np.random.default_rng(0)
    for _ in range(100):
        n = int(rng.integers(1, 9))
        rows = rng.integers(0, 2, (int(rng.integers(1, 7)), n)).tolist()
        K = _kernel(rows, n)
        assert len(K) == n - _rank2(rows, n)
        for v in K:
            assert all(sum(a * b for a, b in zip(r, v, strict=True)) % 2 == 0 for r in rows)


def test_sqclass():
    assert sqclass((4,), 5) == (1,)
    assert sqclass((2,), 5) == (2,)
    assert sqclass(pmul((0, 1), (0, 1), 5), 5) == (1,)
    assert sqclass(pmul((0, 2), (1, 1), 5), 5) == pmul((2,), pmul((0, 1), (1, 1), 5), 5)


@pytest.mark.parametrize("q,n", [(5, 2), (3, 3), (7, 2)])
def test_selmer_classes(q, n):
    for t in family(legendre(q), n):
        M = model(t.D, q)
        S = set(M.selmer)
        assert len(S) == len(M.selmer) == 2 ** selmer_rank(M.e, q)
        assert torsion_classes(M) <= S
        assert len(torsion_classes(M)) == 4


@pytest.mark.parametrize("q", [5, 7, 11, 13])
def test_rational_surface_generator(q):
    seen = 0
    for t in family(legendre(q), 1):
        if t.rank != 1:
            continue
        r = rank_one(t)
        assert r.chi == 1
        assert r.generator.height == Fraction(1, 2)
        assert r.sha == 1
        seen += 1
    assert seen


def test_every_rank_one_twist_is_solved(solved):
    assert len(solved) == 8 + 5
    assert all(r is not None for _, r in solved)


def test_sha_is_one(solved):
    for _, r in solved:
        assert r.sha == 1
        assert r.sha2 == 0
        assert r.lstar > 0


def test_height_from_doubling_matches_shioda(solved):
    for t, r in solved:
        M = model(t.D, t.L.q)
        P = r.generator.point
        for v, f in M.fibers.items():
            c = M.contr(P, v)
            if f.kind == "I2*":
                assert c in (0, 1, Fraction(3, 2))
            else:
                assert c in (0, MAX_CONTR[f.kind])
            assert (c > 0) == M.singular(P, v)


def test_height_is_quadratic(solved):
    for t, r in solved:
        M = model(t.D, t.L.q)
        P = r.generator.point
        h = M.height(P)
        Q = M.double(P)
        assert psqrt(M.cubic(Q.A, Q.B), M.p) is not None
        assert M.height(Q) == 4 * h
        assert M.height(M.double(Q)) == 16 * h
        for i in range(3):
            R = M.translate(P, i)
            assert psqrt(M.cubic(R.A, R.B), M.p) is not None
            assert M.height(R) == h


def test_delta_is_a_homomorphism(solved):
    for t, r in solved:
        M = model(t.D, t.L.q)
        P, p = r.generator.point, M.p
        d = delta(M, P)
        assert d in M.selmer
        images = {delta(M, M.translate(P, i)) for i in range(3)} | {d}
        expect = {(sqclass(pmul(d[0], a, p), p), sqclass(pmul(d[1], b, p), p))
                  for a, b in torsion_classes(M)}
        assert images == expect
        assert delta(M, M.double(P)) == ((1,), (1,))


def test_generator_is_minimal_in_every_class(solved):
    for t, r in solved[:8]:
        M = model(t.D, t.L.q)
        B = r.generator.bound
        for cls in set(M.selmer) - torsion_classes(M):
            hs = [M.height(P) for P in points(M, B, cls)]
            assert all(h >= r.generator.height for h in hs)


@pytest.mark.parametrize("D", [(2, 1, 0, 1), (1, 1, 2, 1)])
def test_points_match_brute_force(D):
    M = model(D, 3)
    brute = set()
    for w in ((1,), (0, 1), (1, 1), (2, 1)):
        for A in itertools.product(range(3), repeat=2 * M.chi + 3):
            P = Point(norm(A, 3), pmul(w, w, 3))
            c = M.cubic(P.A, P.B)
            if c and len(pgcd(P.A, w, 3)) == 1 and psqrt(c, 3) is not None:
                brute.add((delta(M, P), P))
    fast = {(cls, P) for cls in M.selmer for P in points(M, 1, cls)}
    assert brute
    assert fast == brute


def test_lstar():
    for t in family(legendre(5), 2):
        if t.rank == 0:
            assert lstar(t.L) == t.L.central
        else:
            d = sum(Fraction(k * ck, 5 ** (k - 1)) for k, ck in enumerate(t.L.c) if k)
            assert lstar(t.L) == -d / 5


@pytest.mark.parametrize("q,D", [(7, (5, 1, 2, 1)), (5, (1, 0, 2, 4, 1)), (5, (3, 0, 1, 4, 1))])
def test_sha_four_when_sha2_has_rank_two(q, D):
    t = next(t for t in family(legendre(q), len(D) - 1) if t.D == D)
    r = rank_one(t, 1)
    assert r.selmer == 5 and r.sha2 == 2
    assert r.sha == 4
