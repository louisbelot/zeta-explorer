from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.brauer import ratio
from zeta_explorer.funcfield import family, legendre
from zeta_explorer.gf import peval
from zeta_explorer.heights import Point, model, rank_one, torsion_classes
from zeta_explorer.mordellweil import (
    add,
    basis,
    cosets,
    det,
    gram,
    kmul,
    lift,
    pairing,
    reduce,
    solve,
)


def twists(q, n, rank, Ds=None):
    return [t for t in family(legendre(q), n)
            if t.rank == rank and (Ds is None or t.D in Ds)]


@pytest.fixture(scope="module")
def rank_two():
    ts = twists(3, 5, 2) + twists(7, 2, 2) + twists(11, 2, 2, {(1, 7, 1), (6, 6, 1)})
    return [(t, solve(t)) for t in ts]


@pytest.fixture(scope="module")
def rank_three():
    return [(t, solve(t, 1)) for t in twists(13, 2, 3)]


def test_det():
    rng = np.random.default_rng(5)
    for n in range(1, 5):
        for _ in range(20):
            G = rng.integers(-4, 5, (n, n))
            got = det([[Fraction(int(x), 3) for x in row] for row in G])
            assert abs(float(got) - np.linalg.det(G / 3)) < 1e-9


def test_cosets():
    for t in family(legendre(7), 3):
        M = model(t.D, 7)
        tors = torsion_classes(M)
        reps = cosets(M)
        assert len(reps) == len(M.selmer) // 4 - 1
        cover = {kmul(c, s, 7) for c in reps for s in tors} | tors
        assert cover == set(M.selmer)


def test_every_rank_two_twist_is_solved(rank_two):
    assert len(rank_two) == 8 + 2 + 2
    assert all(r is not None for _, r in rank_two)


def test_sha_is_one_in_rank_two(rank_two):
    for t, r in rank_two:
        assert r.sha == 1
        assert r.sha2 == 0
        assert ratio(t).sha_reg == r.basis.regulator


def test_basis_is_reduced(rank_two):
    for _, r in rank_two:
        (a, b), (c, d) = r.basis.gram
        assert b == c
        assert 0 < a <= d
        assert 2 * abs(b) <= a
        assert r.basis.regulator == a * d - b * b


def test_reduce_undoes_a_skewed_basis(rank_two):
    for t, r in rank_two[:6]:
        M = model(t.D, t.L.q)
        P, Q = r.basis.points
        for _ in range(2):
            Q = max((add(M, Q, P, s) for s in (1, -1)), key=M.height)
        assert M.height(Q) > r.basis.gram[1][1]
        G = gram(M, reduce(M, [Q, P]))
        assert det(G) == r.basis.regulator
        assert [G[0][0], G[1][1], abs(G[0][1])] == [r.basis.gram[0][0], r.basis.gram[1][1],
                                                    abs(r.basis.gram[0][1])]


def test_known_lattices(rank_two):
    reg = {(r.q, r.D): r.basis.regulator for _, r in rank_two}
    assert reg[7, (5, 0, 1)] == Fraction(3, 2)
    assert reg[7, (4, 6, 1)] == 7
    assert reg[11, (1, 7, 1)] == Fraction(7, 4)
    assert sorted(set(v for (q, _), v in reg.items() if q == 3)) == [4, 14, 20]


def test_rank_three(rank_three):
    assert len(rank_three) == 2
    for _, r in rank_three:
        assert r.basis.regulator == Fraction(13, 8)
        assert r.sha == 1
        G = r.basis.gram
        assert all(G[i][j] == G[j][i] for i in range(3) for j in range(3))


def test_parallelogram(rank_two):
    for t, r in rank_two:
        M = model(t.D, t.L.q)
        P, Q = r.basis.points
        hs = M.height(add(M, P, Q)) + M.height(add(M, P, Q, -1))
        assert hs == 2 * M.height(P) + 2 * M.height(Q)


def test_pairing_ignores_torsion(rank_two):
    for t, r in rank_two[:6]:
        M = model(t.D, t.L.q)
        P, Q = r.basis.points
        b = abs(pairing(M, P, Q))
        for i in range(3):
            assert abs(pairing(M, P, M.translate(Q, i))) == b


def test_add_translates_rank_one_generators():
    for t in twists(5, 2, 1):
        M = model(t.D, 5)
        P = rank_one(t).generator.point
        for i in range(3):
            T = Point(M.e[i], (1,))
            assert add(M, P, T) == add(M, P, T, -1) == M.translate(P, i)
        with pytest.raises(ArithmeticError):
            add(M, P, P)


def test_add_matches_translation(rank_two):
    for t, r in rank_two:
        M = model(t.D, t.L.q)
        P = r.basis.points[0]
        for i in range(3):
            assert add(M, P, Point(M.e[i], (1,))) == M.translate(P, i)


def chord(x1, y1, x2, y2, a2, p):
    lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    return (lam * lam - a2 - x1 - x2) % p


def test_add_reduces_to_the_fiber_group_law(rank_two):
    seen = 0
    for t, r in rank_two:
        q = t.L.q
        M = model(t.D, q)
        P, Q = r.basis.points
        S = add(M, P, Q)
        (A1, w1, C1), (A2, w2, C2) = lift(M, P), lift(M, Q)
        for t0 in range(q):
            e = [peval(ei, t0, q) for ei in M.e]
            dens = [peval(w, t0, q) for w in (w1, w2)] + [peval(S.B, t0, q)]
            if len(set(e)) < 3 or 0 in dens:
                continue
            u1, u2 = (pow(peval(w, t0, q), -1, q) for w in (w1, w2))
            x1, y1 = peval(A1, t0, q) * u1**2 % q, peval(C1, t0, q) * u1**3 % q
            x2, y2 = peval(A2, t0, q) * u2**2 % q, peval(C2, t0, q) * u2**3 % q
            if x1 == x2:
                continue
            x3 = peval(S.A, t0, q) * pow(dens[2], -1, q) % q
            assert x3 == chord(x1, y1, x2, y2, -sum(e) % q, q)
            seen += 1
    assert seen > 20


def test_rank_one_agrees_with_generator():
    for t in twists(5, 2, 1):
        a, b = solve(t), rank_one(t)
        assert a.basis.regulator == b.generator.height
        assert a.sha == b.sha == 1


def test_rank_four_is_refused():
    t = twists(7, 2, 2)[0]
    with pytest.raises(ValueError):
        basis(model(t.D, 7), 4)
