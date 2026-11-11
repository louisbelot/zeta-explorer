import math

import mpmath
import numpy as np
import pytest

from zeta_explorer import barrier as P

N0, N1 = 69_098, 1_500_000


def test_erfi_tail_is_the_integral():
    sigma, t = 1.3, 0.2
    exact = mpmath.quad(lambda a: mpmath.exp(t / 4 * mpmath.log(a) ** 2) / a**sigma, [100, 5000])
    assert P.erfi_tail(100, 5000, sigma, t) == pytest.approx(float(exact), rel=1e-12)


def test_tails_bound_the_sum():
    sigma, t, M, N = 1.3, 0.2, 100, 5000
    total = P.bsum(N, sigma, t)
    head = P.bsum(M, sigma, t)
    assert total <= head + P.erfi_tail(M, N, sigma, t) <= head + P.lemma_tail(M, N, sigma, t)
    assert total - head > P.erfi_tail(M + 1, N, sigma, t)


def test_lemma_needs_room():
    with pytest.raises(ValueError):
        P.lemma_tail(10, 1e6, 0.5, 0.2)


def test_far_region_from_the_paper():
    A, B = P.far_bound(N1, 0.2, 0.2, N0)
    assert B == pytest.approx(0.0740, abs=1e-4)
    assert A == pytest.approx(1.906, abs=1e-3)
    assert P.far_sup(N1, 1e30, 0.2, 0.2, N0) == pytest.approx(A + B, rel=1e-9)
    assert A + B < 2 - 1.25e-3


def test_top_edge_from_the_paper():
    assert 0.6 < P.top_bound(N0, N1, 0.2) < 0.8


def test_mollifier_range_from_the_paper():
    sigma = P.sigma_floor(N0, 0.2, 0.2)
    assert sigma == pytest.approx(1.7143, abs=1e-4)
    lo, hi, arg = P.mollifier_range(sigma, 0.2)
    assert (lo, hi, arg) == pytest.approx((0.534, 1.635, 0.553), abs=2e-3)


def test_mollifier_kills_small_primes():
    t = 0.2
    n = np.arange(1, 200)
    beta, lam = P.coefficients(n, np.inf, t, 0.2)
    assert beta[0] == 1
    assert np.allclose(beta[[1, 2, 4]], 0, atol=1e-15)
    coprime = np.gcd(n, 30) == 1
    assert np.allclose(beta[coprime], np.exp(P.log_b(n[coprime], t)))
    assert np.allclose(lam[coprime], beta[coprime] * n[coprime] ** 0.2)


def test_mobius_at_time_zero():
    n = np.arange(1, 400)
    full, _ = P.coefficients(n, np.inf, 0.0, 0.2)
    assert np.allclose(full, np.gcd(n, 30) == 1)
    cut, _ = P.coefficients(np.array([210]), 50, 0.0, 0.2)
    assert cut[0] == pytest.approx(1)


def test_lower_edge_from_the_paper():
    assert P.lower_edge(N0, 80_000, every_n=False, c=1.03) == pytest.approx(0.0263, abs=1e-3)
    assert P.lower_edge(N0, 80_000) > 2.15e-3
    assert P.lower_edge(N0, 80_000, improved=False) < 0


def test_polar_distance():
    assert P.polar_distance(2.0, 1.0) == 2.0
    assert P.polar_distance(2.0, 3 * math.pi / 4) == pytest.approx(math.sqrt(2))
    assert P.polar_distance(2.0, 4.0) == 0.0
