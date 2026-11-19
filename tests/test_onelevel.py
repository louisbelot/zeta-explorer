import itertools
from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.funcfield import Curve, family, legendre
from zeta_explorer.gf import deg, monics, pgcd, squarefree
from zeta_explorer.hypermoments import _symbol, irreducibles, lpoly
from zeta_explorer.onelevel import (
    char_sum,
    density_gap,
    exact_traces,
    family_traces,
    limit_traces,
    ms_traces,
    multiplicative,
    power_sums,
    rmt_trace,
    sym2_degree,
    sym2_lpoly,
    sym2_sums,
)

E2 = Curve(3, (-1, 0, -1), (0, 0, 1), ())
E3 = Curve(5, (-1, 0, 0, -1), (0, 0, 0, 1), ())


def test_power_sums():
    assert power_sums((1, -5, 6), 3) == [5, 13, 35]


def _haar(dim, size, rng):
    Z = rng.standard_normal((size, dim, dim))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]


@pytest.mark.parametrize("dim", [6, 7])
def test_rmt_trace_against_haar(dim):
    Q = _haar(dim, 40_000, np.random.default_rng(dim))
    lam = np.linalg.eigvals(Q)
    det = np.sign(np.linalg.det(Q))
    for n in range(1, dim + 3):
        tr = (lam**n).sum(axis=1).real
        for s in (1, -1):
            x = tr[det == s]
            assert abs(x.mean() - rmt_trace(dim, s, n)) < 5 * x.std() / np.sqrt(len(x)) + 1e-9


def test_rmt_trace_needs_positive_n():
    with pytest.raises(ValueError):
        rmt_trace(4, 1, 0)


def _brute(q, N, f, S):
    total = 0
    for D in monics(q, N):
        if squarefree(D, q) and all(deg(pgcd(D, Q, q)) == 0 for Q in S):
            total += np.prod([_symbol(D, P, q) for P in f], dtype=np.int64)
    return int(total)


@pytest.mark.parametrize("q, N", [(3, 4), (5, 3)])
def test_char_sum_against_enumeration(q, N):
    S = [(0, 1), (q - 1, 1)]
    assert char_sum(q, N, None, [], [(1, 1), (1, 1)]) == _brute(q, N, [], S)
    for P in [(1, 1), irreducibles(q, 2)[0], irreducibles(q, 3)[-1]]:
        avoid = [(1, _symbol(Q, P, q)) for Q in S]
        assert char_sum(q, N, lpoly(q, [P])[0], [deg(P)], avoid) == _brute(q, N, [P], S)


@pytest.mark.parametrize("E, N", [(legendre(3), 3), (legendre(3), 4), (legendre(5), 2),
                                  (legendre(5), 3), (E2, 3)])
def test_exact_traces_match_the_family(E, N):
    fam = family(E, N)
    n = 2 * N + 1
    assert exact_traces(E, N, n) == family_traces(fam, n)


@pytest.mark.parametrize("E, N, shift", [(legendre(3), 3, 1), (legendre(3), 4, 1),
                                         (legendre(5), 2, 1), (legendre(5), 3, -1), (E2, 3, 1)])
def test_sign_classes_match_root_numbers(E, N, shift):
    fam = family(E, N)
    for s in (1, -1):
        assert exact_traces(E, N, 2 * N, s) == family_traces(fam, 2 * N, s * shift)


@pytest.mark.parametrize("N", [2, 3, 4])
def test_character_at_m_is_the_root_number_over_f3(N):
    E = legendre(3)
    fam = family(E, N)
    M = multiplicative(E)
    for t in fam:
        chi = np.prod([_symbol(t.D, Q, 3) for Q in M])
        assert chi == t.L.eps


def test_sign_needs_multiplicative_places():
    E = Curve(5, (), (), (0, 1))
    assert multiplicative(E) == []
    with pytest.raises(ValueError):
        exact_traces(E, 2, 2, 1)


def test_odd_traces_cancel_when_q_is_3_mod_4():
    assert all(x == 0 for x in exact_traces(legendre(3), 4, 9)[::2])
    assert any(x != 0 for x in exact_traces(legendre(5), 4, 5)[::2])


@pytest.mark.parametrize("q", [3, 5, 7])
def test_legendre_sym2_is_trivial(q):
    assert sym2_degree(legendre(q)) == 0
    assert sym2_lpoly(legendre(q)) == (1,)
    assert sym2_sums(legendre(q), 5) == [0] * 5


@pytest.mark.parametrize("E, d", [(E2, 2), (E3, 4)])
def test_sym2_riemann_hypothesis(E, d):
    c = sym2_lpoly(E)
    q = E.p
    assert len(c) == d + 1
    eps = 1 if c[d] == q ** (3 * d // 2) * c[0] else -1
    for k in range(d // 2 + 1):
        assert c[d - k] == eps * q ** (3 * (d - 2 * k) // 2) * c[k]
    gamma = 1 / np.roots(np.array(c[::-1], dtype=float))
    assert np.allclose(np.abs(gamma), q**1.5)


@pytest.mark.parametrize("E, n", [(E2, 7), (E3, 6)])
def test_sym2_newton_continues_past_the_degree(E, n):
    assert power_sums(sym2_lpoly(E), n) == [-s for s in sym2_sums(E, n)]


def test_sym2_rejects_potentially_good_fibers():
    with pytest.raises(ValueError):
        sym2_degree(Curve(5, (), (), (0, 1)))


@pytest.mark.parametrize("E, N", list(itertools.product([legendre(3), E2], [3, 4])))
def test_ms_formula_is_the_limit(E, N):
    assert limit_traces(E, N, 8) == ms_traces(E, N, 8)


def test_exact_converges_to_the_limit():
    E = legendre(3)
    ex, lim = exact_traces(E, 20, 6), limit_traces(E, 20, 6)
    assert max(abs(float(a - b)) for a, b in zip(ex, lim, strict=True)) < 1e-6
    far = exact_traces(E, 4, 6)
    assert max(abs(float(a - b)) for a, b in zip(far, lim, strict=True)) > 1e-3


def test_lower_order_term_at_n_2():
    assert limit_traces(legendre(3), 4, 2)[1] == Fraction(5, 18)


def test_density_gap():
    tr = [Fraction(0), Fraction(1, 2), Fraction(0), Fraction(1)]
    assert density_gap(tr, 8, lambda u: 1.0) == pytest.approx(2 / 8 * -0.5)
