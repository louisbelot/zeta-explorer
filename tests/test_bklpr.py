import itertools
from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.bklpr import (
    Space,
    alternating,
    alternating_count,
    central_integer,
    det_valuation,
    divides,
    first_step,
    kernel_law,
    nonsingular,
    odd_sha,
    orthogonal_valuations,
    special_orthogonal,
    square_class,
    valuation,
    vp,
)
from zeta_explorer.bsd import arith
from zeta_explorer.funcfield import family, legendre


def test_vp():
    assert vp(np.array([0, 1, 3, 9, 18, 81]), 3, 3).tolist() == [3, 0, 1, 2, 2, 3]


def test_det_valuation_against_integer_determinants():
    rng = np.random.default_rng(0)
    M = rng.integers(-4, 5, (400, 4, 4))
    exact = [round(np.linalg.det(m)) for m in M]
    want = [min(valuation(d, 3), 6) if d else 6 for d in exact]
    assert det_valuation(M, 3, 6).tolist() == want


@pytest.mark.parametrize("n,q", [(4, 2), (4, 3), (5, 3), (6, 2)])
def test_macwilliams_counts(n, q):
    total = sum(alternating_count(n, r, q) for r in range(n // 2 + 1))
    assert total == q ** (n * (n - 1) // 2)


def test_macwilliams_by_brute_force():
    pairs = list(itertools.combinations(range(4), 2))
    counts = [0, 0, 0]
    for entries in itertools.product(range(3), repeat=6):
        A = np.zeros((4, 4), dtype=np.int64)
        for (i, j), a in zip(pairs, entries, strict=True):
            A[i, j], A[j, i] = a, -a
        r = 0 if not A.any() else (2 if round(np.linalg.det(A)) % 3 == 0 else 4)
        counts[r // 2] += 1
    assert counts == [alternating_count(4, r, 3) for r in range(3)]


def test_kernel_law_is_the_limit_of_counts():
    n, q = 20, 3
    for k in range(3):
        finite = alternating_count(n, n // 2 - k, q) / q ** (n * (n - 1) // 2)
        assert finite == pytest.approx(kernel_law(k, q), rel=1e-6)
    assert sum(kernel_law(k, 5) for k in range(10)) == pytest.approx(1)


def test_alternating_nonsingular_fraction():
    v = det_valuation(alternating(6, 3, 1, 20_000, seed=1), 3, 1)
    assert np.mean(v == 0) == pytest.approx(nonsingular(6, 3), abs=0.012)


def test_delaunay_first_step_against_the_matrix_model():
    v = det_valuation(alternating(24, 3, 6, 6000, seed=2), 3, 6)
    assert np.all(v % 2 == 0)
    assert np.mean(v == 0) == pytest.approx(1 - divides(3), abs=0.02)
    assert np.mean(v == 2) == pytest.approx(first_step(3), abs=0.02)
    assert first_step(3) == pytest.approx((1 - divides(3)) * 3 / 8)


def test_divides_values():
    assert 1 - divides(3) == pytest.approx(0.639005, abs=1e-6)
    assert 1 - divides(5) == pytest.approx(0.793335, abs=1e-6)
    assert divides(3, rank=1) < divides(3)


def _so_exact(V):
    m, G = V.p**V.K, V.gram
    out = []
    for x in itertools.product(range(m), repeat=V.n * V.n):
        g = np.array(x).reshape(V.n, V.n)
        if ((g.T * G) @ g % m == np.diag(G)).all() and round(np.linalg.det(g)) % m == 1:
            out.append(g)
    return np.array(out)


@pytest.mark.parametrize("V", [Space(2, 3, 2, 1), Space(2, 3, 2, -1), Space(3, 3, 1, 1)])
def test_orthogonal_sampler_is_haar(V):
    exact = _so_exact(V)
    eye = np.eye(V.n, dtype=np.int64)
    want = np.bincount(det_valuation((eye - exact) % V.p**V.K, V.p, V.K), minlength=V.K + 1)
    got = np.bincount(orthogonal_valuations(V, 30_000, seed=4), minlength=V.K + 1)
    assert got / got.sum() == pytest.approx(want / want.sum(), abs=0.012)


def test_sampler_lands_in_the_group_and_coset():
    V = Space(4, 5, 2, -1)
    m = V.p**V.K
    eye = np.eye(4, dtype=np.int64)
    for theta in (1, -1):
        g, th = special_orthogonal(V, 300, spinor=theta, seed=5)
        assert np.all(th == theta)
        gram = np.einsum("sji,j,sjk->sik", g, V.gram, g) % m
        assert np.all(gram == np.diag(V.gram))
        for h in g:
            d = round(np.linalg.det((eye + h).astype(float))) % m
            if d % 5:
                assert square_class(Fraction(d), 5)[1] == theta


def test_square_class():
    assert square_class(Fraction(18), 3) == (0, -1)
    assert square_class(Fraction(4, 27), 3) == (1, 1)


@pytest.fixture(scope="module")
def fam5():
    return family(legendre(5), 4)


def test_odd_part_of_sha_is_the_odd_part_of_q_n_l(fam5):
    for t in fam5:
        if t.L.eps == 1 and t.rank == 0:
            sha = arith(t).sha
            assert sha.denominator == 1
            for p in (3, 7):
                assert valuation(int(sha), p) == valuation(central_integer(t), p)


def test_odd_sha_summary(fam5):
    S = odd_sha(fam5, 3, size=4000)
    assert S.twists == 135
    assert S.observed.sum() == pytest.approx(1)
    assert S.bklpr[0] == pytest.approx(1 - divides(3))
    assert S.orthogonal[0] == pytest.approx(S.bklpr[0], abs=0.04)
    with pytest.raises(ValueError):
        odd_sha(fam5, 5)
