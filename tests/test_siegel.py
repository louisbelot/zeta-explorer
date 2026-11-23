import math

import numpy as np
import pytest

from zeta_explorer.characters import is_fundamental, kronecker
from zeta_explorer.siegel import (
    bound,
    constants,
    full_sum,
    fundamental,
    kronecker_column,
    nth_primes,
    scan,
    weights,
)

TABLE = {1.6: (0.22675, 1.1614), 1.3: (0.23083, 1.1805), 1.2: (0.23222, 1.1870),
         1.1: (0.23362, 1.1936)}


def r_of(lam):
    return lam / (10 * math.log(10))


@pytest.mark.parametrize("lam", sorted(TABLE))
def test_constants_match_table(lam):
    k = constants(r_of(lam))
    phi, E = TABLE[lam]
    assert math.ceil(k.phi * 1e5) / 1e5 == pytest.approx(phi, abs=1e-12)
    assert math.ceil(k.E * 1e4) / 1e4 == pytest.approx(E, abs=1e-12)


def test_phi_tends_to_quarter():
    assert constants(1e-6).phi == pytest.approx(0.25, abs=1e-5)
    rs = [0.02, 0.05, 0.1, 0.2]
    phis = [constants(r).phi for r in rs]
    assert phis == sorted(phis, reverse=True)


def test_fundamental_matches_scan():
    X = 3000
    want = sorted(d for d in range(-X, X + 1) if is_fundamental(d))
    assert sorted(fundamental(X).tolist()) == want


def test_fundamental_density():
    X = 10**6
    d = fundamental(X)
    assert len(d) == 607_925
    assert len(d) / X == pytest.approx(6 / math.pi**2, rel=2e-3)


def test_kronecker_column():
    rng = np.random.default_rng(1)
    d = fundamental(5000)
    for p in nth_primes(40):
        col = kronecker_column(d, int(p))
        sample = rng.choice(d.size, 200, replace=False)
        assert [int(col[i]) for i in sample] == [kronecker(int(d[i]), int(p)) for i in sample]


def test_kronecker_column_euler_path():
    d = fundamental(400)
    p = 1_000_003
    assert [int(x) for x in kronecker_column(d, p)] == [kronecker(int(x), p) for x in d]


def test_weights_sum_to_log_derivatives():
    ps = nth_primes(20000)
    sigma = 2.0
    for d in (-3, 5, -163, 12):
        chi = np.array([kronecker(d, int(p)) for p in ps])
        S = weights(ps, sigma)[np.arange(ps.size), chi + 1].sum()
        assert S == pytest.approx(float(full_sum(d, sigma)), abs=1e-5)


@pytest.mark.parametrize("d", [-3, -4, 5, -163, -420, 1009])
def test_bound_holds_without_real_zero(d):
    lam = 1.6
    k = constants(r_of(lam))
    S = float(full_sum(d, 1 + k.r))
    assert S <= 1 / k.r + k.phi * math.log(abs(d)) + k.E


def test_bound_is_theorem_at_the_edge():
    k = constants(r_of(1.3))
    R = k.r + 7 / 8
    for q in (10**5, 10**8, 10**10):
        L = math.log(q)
        beta = 1 - 0.2 / L
        gap = 1 + k.r - beta
        want = 1 / k.r - 1 / gap + gap / R**2 + k.phi * L + k.E
        assert float(bound(q, k)) == pytest.approx(want, rel=1e-13)


def test_small_discriminants_fail_the_criterion():
    k = constants(r_of(1.6))
    s = scan(np.array([-3, -4, 5, -163]), k, 60000)
    assert (s.needed == 0).all()
