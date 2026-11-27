from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.funcfield import family, legendre
from zeta_explorer.heights import Point, model, rank_one
from zeta_explorer.mordellweil import add, det


def twists(q, n, rank):
    return [t for t in family(legendre(q), n) if t.rank == rank]


def test_det():
    rng = np.random.default_rng(5)
    for n in range(1, 5):
        for _ in range(20):
            G = rng.integers(-4, 5, (n, n))
            got = det([[Fraction(int(x), 3) for x in row] for row in G])
            assert abs(float(got) - np.linalg.det(G / 3)) < 1e-9


def test_add_translates_rank_one_generators():
    for t in twists(5, 2, 1):
        M = model(t.D, 5)
        P = rank_one(t).generator.point
        for i in range(3):
            T = Point(M.e[i], (1,))
            assert add(M, P, T) == add(M, P, T, -1) == M.translate(P, i)
        with pytest.raises(ArithmeticError):
            add(M, P, P)
