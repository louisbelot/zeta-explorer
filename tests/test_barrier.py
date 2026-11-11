import math

import mpmath
import numpy as np
import pytest

from zeta_explorer import barrier as P
from zeta_explorer.effective import _log_mt_mp

N0, N1 = 69_098, 1_500_000


def reference_f(X, u, y, t, N):
    with mpmath.workdps(40):
        s = (1 + mpmath.mpf(y) - 1j * (mpmath.mpf(X) + mpmath.mpf(u))) / 2
        tm = mpmath.mpf(t)

        def a(z):
            return 1 / (2 * z) + 1 / (z - 1) + mpmath.log(z / (2 * mpmath.pi)) / 2

        s1, s2 = s + tm / 2 * a(s), 1 - s + tm / 2 * a(1 - s)
        g = mpmath.exp(_log_mt_mp(1 - s, tm) - _log_mt_mp(s, tm))
        acc = mpmath.mpc(0)
        for n in range(1, N + 1):
            ln = mpmath.log(n)
            acc += mpmath.exp(tm / 4 * ln * ln) * (mpmath.exp(-s1 * ln) + g * mpmath.exp(-s2 * ln))
        return complex(acc)


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


def test_series_matches_direct_sum():
    X = 1e6 + 0.3
    bar = P.Barrier(X=X)
    for u, y, t in [(-0.5, 0.2, 0.0), (0.31, 0.77, 0.013), (0.5, 1.0, 0.2), (-0.2, 0.25, 0.157)]:
        assert abs(bar.f(u, y, t)[0] - reference_f(X, u, y, t, bar.N)) < 1e-13


def test_log_gamma_matches_mpmath():
    bar = P.Barrier(X=P.X0)
    for u, y, t in [(-0.5, 0.2, 0.0), (0.4, 0.9, 0.2)]:
        with mpmath.workdps(40):
            s = (1 + mpmath.mpf(y) - 1j * (mpmath.mpf(P.X0) + mpmath.mpf(u))) / 2
            d = _log_mt_mp(1 - s, mpmath.mpf(t)) - _log_mt_mp(s, mpmath.mpf(t))
            ref = complex(float(d.real), float(d.imag % (2 * mpmath.pi)))
        got = complex(bar.log_gamma(u, y, t))
        assert abs(got.real - ref.real) < 1e-13
        assert abs((got.imag - ref.imag + math.pi) % (2 * math.pi) - math.pi) < 1e-12


def test_bounds_dominate_differences():
    bar = P.Barrier(X=1e6 + 0.3)
    rng = np.random.default_rng(3)
    for t1, t2 in [(0.0, 0.01), (0.15, 0.16)]:
        dz, dt = bar.bounds(t1, t2)
        u, y = rng.uniform(-0.5, 0.5, 300), rng.uniform(0.2, 1.0, 300)
        t = (t1 + t2) / 2
        fz = np.abs(bar.f(u + 1e-6, y, t) - bar.f(u - 1e-6, y, t)) / 2e-6
        ft = np.abs(bar.f(u, y, t + 1e-6) - bar.f(u, y, t - 1e-6)) / 2e-6
        assert fz.max() < dz
        assert ft.max() < dt


def test_error_at_the_barrier():
    bar = P.Barrier()
    assert bar.err(0.0, 0.0) == pytest.approx(1.249e-3, rel=1e-3)


def test_barrier_location_is_large():
    bar = P.Barrier()
    assert bar.N == N0
    assert np.abs(bar.f([-0.5, 0.0, 0.5], 1.0, 0.0)).min() == pytest.approx(4.32, abs=0.01)
