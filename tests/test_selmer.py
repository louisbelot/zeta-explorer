import numpy as np
import pytest

from zeta_explorer.bsd import legendre_roots, selmer_rank
from zeta_explorer.funcfield import family, legendre, twist_set
from zeta_explorer.selmer import (
    chain,
    equilibrium,
    histogram,
    law,
    markov,
    poonen_rains,
    random_symbols,
    symbolic_rank,
    symbols,
)


@pytest.mark.parametrize("p", [2, 3, 5])
def test_poonen_rains_mean(p):
    P = np.array([poonen_rains(n, p) for n in range(40)])
    assert P.sum() == pytest.approx(1, abs=1e-12)
    assert (P * float(p) ** np.arange(40)).sum() == pytest.approx(p + 1, abs=1e-9)


@pytest.mark.parametrize("parity", [0, 1])
def test_equilibrium_by_parity(parity):
    E = np.array([equilibrium(n, parity) for n in range(40)])
    assert E.sum() == pytest.approx(1, abs=1e-12)
    assert (E * 2.0 ** np.arange(40)).sum() == pytest.approx(3, abs=1e-9)
    assert E[1 - parity::2].sum() == 0


def test_heath_brown_values():
    assert equilibrium(0, 0) == pytest.approx(0.419422, abs=1e-6)
    assert equilibrium(2, 0) / equilibrium(0, 0) == pytest.approx(4 / 3)
    assert equilibrium(3, 1) / equilibrium(1, 1) == pytest.approx(4 / 21)


def test_markov_is_stochastic_and_stationary():
    M = markov(30)
    assert np.allclose(M[:-1].sum(axis=1), 1)
    for parity in (0, 1):
        x = chain(parity, 60, 30)
        assert np.allclose(x, [equilibrium(n, parity) for n in range(30)], atol=1e-10)
    assert chain(0, 1)[1] == 1
    assert chain(1, 1)[[0, 2]] == pytest.approx([0.5, 0.5])


@pytest.mark.parametrize("q,n", [(3, 4), (5, 3), (7, 2), (11, 2), (13, 2)])
def test_symbolic_descent_matches_full_descent(q, n):
    for D in twist_set(q, n, legendre(q).disc):
        assert symbolic_rank(symbols(D, q)) == selmer_rank(legendre_roots(D, q), q)


@pytest.mark.parametrize("q", [3, 5, 7, 11])
def test_symbols_obey_reciprocity(q):
    rng = np.random.default_rng(q)
    real = [symbols(D, q) for D in twist_set(q, 4, legendre(q).disc)[:60]]
    fake = [random_symbols(q, (1, 2, 3, 2), rng) for _ in range(20)]
    for S in real + fake:
        for i in range(S.omega):
            for j in range(i + 1, S.omega):
                assert (S.mutual[i][j] + S.mutual[j][i]) % 2 == S.sign(S.degs[i] * S.degs[j])


def test_random_symbols_reach_equilibrium():
    rng = np.random.default_rng(3)
    s = np.array([symbolic_rank(random_symbols(5, (1,) * 20, rng)) - 2 for _ in range(1500)])
    even = s[s % 2 == 0]
    assert np.mean(s % 2) == pytest.approx(0.5, abs=0.05)
    assert np.mean(even == 0) == pytest.approx(equilibrium(0, 0), abs=0.04)
    assert np.mean(even == 2) == pytest.approx(equilibrium(2, 0), abs=0.04)


def test_law_of_a_family():
    L = law(family(legendre(7), 3), draws=8)
    assert L.twists == 227
    assert L.observed.sum() == pytest.approx(1)
    assert L.given(0)[1::2].sum() == 0
    assert L.given(1, "predicted")[0::2].sum() == 0
    assert 0 <= L.tv <= 1
    assert L.mean == pytest.approx(np.dot(L.observed, 2.0 ** np.arange(8)))


def test_histogram_caps_the_tail():
    assert histogram([0, 1, 9, 9], 4).tolist() == [0.25, 0.25, 0, 0.5]
