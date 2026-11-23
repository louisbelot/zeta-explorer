import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.characters import is_fundamental, kronecker, totient
from zeta_explorer.siegel import (
    TOL,
    bound,
    character_table,
    constants,
    flint_character,
    full_sum,
    fundamental,
    kronecker_column,
    nth_primes,
    positive_on,
    resum,
    scan,
    scan_one,
    verify,
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


def test_character_table():
    ps = nth_primes(400)
    for d in fundamental(2000):
        want = np.array([kronecker(int(d), int(p)) for p in ps])
        assert (character_table(int(d))[ps % abs(d)] == want).all()


def test_character_table_is_multiplicative():
    for d in (-4, 8, -8, 12, -24, 40, -163, 1005):
        assert is_fundamental(d)
        t = character_table(d).astype(int)
        q = abs(d)
        a, b = np.meshgrid(np.arange(q), np.arange(q))
        assert (t[a * b % q] == t[a] * t[b]).all()
        assert t.sum() == 0


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


def test_scan_matches_scan_one():
    k = constants(r_of(1.6))
    d = fundamental(3000)
    d = d[np.abs(d) > 300]
    ps = nth_primes(20000)
    s = scan(d, k, 20000, ps=ps)
    w = weights(ps, 1 + k.r)
    for i in range(0, d.size, 37):
        n, m = scan_one(int(d[i]), k, ps, w, chunk=1000)
        assert n == s.needed[i]
        assert m == pytest.approx(s.margin[i], abs=1e-9)


def test_clearing_is_certified_in_arb():
    k = constants(r_of(1.6))
    d = fundamental(20000)
    d = d[np.abs(d) > 5000]
    s = scan(d, k, 30000)
    cleared = np.flatnonzero(s.needed)
    for i in cleared[np.argsort(s.margin[cleared])[:3]]:
        S = resum(int(d[i]), k.r, int(s.needed[i]))
        target = float(bound(abs(int(d[i])), k)) + TOL
        assert S > target
        assert float(S.mid()) - target + TOL == pytest.approx(s.margin[i], abs=1e-9)


def test_small_discriminants_fail_the_criterion():
    k = constants(r_of(1.6))
    s = scan(np.array([-3, -4, 5, -163]), k, 60000)
    assert (s.needed == 0).all()


def test_flint_character():
    for d in (-3, 5, -8, 8, -163, 12, -24):
        chi = flint_character(d)
        assert chi.modulus() == abs(d) and chi.is_real() and chi.conductor() == abs(d)


def test_l1_class_number_formula():
    for d, h, w in ((-163, 1, 2), (-3, 1, 6), (-4, 1, 4), (-23, 3, 2), (-232, 2, 2)):
        v = flint_character(d).l(1)
        assert float(v.real.mid()) == pytest.approx(2 * math.pi * h / (w * math.sqrt(-d)))


@pytest.mark.parametrize("d", [-3, 5, -163, 232])
def test_derivative_bound(d):
    vals = [kronecker(d, a) for a in range(abs(d))]
    for s in np.linspace(0.3, 0.98, 12):
        assert abs(mpmath.dirichlet(float(s), vals, 1)) <= totient(abs(d)) / (math.e * s)


def test_positive_on_small():
    for d in (-3, -4, 5, -7, 8, -163, -232):
        assert positive_on(d, 1 - 0.2 / math.log(abs(d)), 1.0)
    assert positive_on(-3, 0.05, 1.0)


def test_verify_small():
    v = verify(1000, passes=((1.6, 60_000), (1.3, 2_000_000)))
    assert v.ok
    direct = set(v.d[v.stage == 3].tolist())
    assert {-3, -4, 5, -163} <= direct
    assert v.counts()[0] > 0.8 * v.d.size


def test_larger_r_suits_small_q():
    ps = nth_primes(200_000)
    for lam, clears in ((3.5, True), (1.3, False)):
        k = constants(r_of(lam))
        n, _ = scan_one(-163, k, ps, weights(ps, 1 + k.r))
        assert (n > 0) == clears
