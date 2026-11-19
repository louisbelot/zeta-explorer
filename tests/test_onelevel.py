
import numpy as np
import pytest

from zeta_explorer.funcfield import Curve, family, legendre
from zeta_explorer.gf import deg, monics, pgcd, squarefree
from zeta_explorer.hypermoments import _symbol, irreducibles, lpoly
from zeta_explorer.onelevel import (
    char_sum,
    exact_traces,
    family_traces,
    multiplicative,
    power_sums,
    rmt_trace,
)

E2 = Curve(3, (-1, 0, -1), (0, 0, 1), ())


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
