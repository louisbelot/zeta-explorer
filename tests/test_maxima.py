import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.maxima import (
    TWO_PI,
    bessel_k,
    cue_maxima,
    fhk_cdf,
    fhk_density,
    haar_angles,
    ks_distance,
    log_char_poly,
    prime_shift_sample,
    prime_shift_variance,
    window_prime_sum,
    zeta_maxima,
)
from zeta_explorer.moments import EULER_GAMMA
from zeta_explorer.multieval import ZWindow


def test_bessel_k_matches_mpmath():
    x = np.array([0.005, 0.1, 1.0, 3.7, 20.0, 60.0])
    for nu in (0, 1):
        ref = np.array([float(mpmath.besselk(nu, v)) for v in x])
        assert np.allclose(bessel_k(nu, x), ref, rtol=1e-12, atol=0)


def test_fhk_law_moments():
    y = np.linspace(-4, 14, 36001)
    h = y[1] - y[0]
    p = fhk_density(y)
    assert p.sum() * h == pytest.approx(1, abs=1e-8)
    assert (y * p).sum() * h == pytest.approx(EULER_GAMMA, abs=1e-8)
    assert ((y - EULER_GAMMA) ** 2 * p).sum() * h == pytest.approx(math.pi**2 / 12, abs=1e-7)


def test_fhk_cdf_is_integral_of_density():
    y = np.linspace(-2, 5, 701)
    dF = np.gradient(fhk_cdf(y), y)
    assert np.max(np.abs(dF - fhk_density(y))[5:-5]) < 1e-4
    assert fhk_cdf(np.array([-6.0]))[0] < 1e-30
    assert fhk_cdf(np.array([12.0]))[0] == pytest.approx(1, abs=1e-9)


def test_fhk_law_is_two_gumbels():
    rng = np.random.default_rng(0)
    y = np.sort(rng.gumbel(size=(2, 200_000)).sum(axis=0) / 2)
    ks = np.max(np.abs(fhk_cdf(y) - np.arange(1, len(y) + 1) / len(y)))
    assert ks < 1.63 / math.sqrt(len(y))


def test_fhk_right_tail():
    y = np.array([6.0, 9.0, 12.0])
    tail = 1 - fhk_cdf(y)
    assert np.allclose(tail / (y * np.exp(-2 * y)), 2, rtol=0.12)


def test_zeta_maxima_against_dense_scan():
    T = 123456.0
    M = zeta_maxima(T, 6)
    win = ZWindow(T - 1, T + 6 * TWO_PI + 1)
    for k in range(6):
        t = np.linspace(T + k * TWO_PI, T + (k + 1) * TWO_PI, 100_001)
        assert M[k] == pytest.approx(np.log(np.abs(win.z(t - win.t0))).max(), abs=1e-8)


def test_zeta_maxima_chunks_agree():
    a = zeta_maxima(5e5, 30)
    b = zeta_maxima(5e5, 30, chunk=7)
    assert np.allclose(a, b, atol=1e-9)


def test_haar_traces():
    phi = haar_angles(6, 20_000, np.random.default_rng(3))
    for j, ref in ((1, 1), (2, 2), (7, 6)):
        tr = np.abs(np.exp(1j * j * phi).sum(axis=1)) ** 2
        assert tr.mean() == pytest.approx(ref, rel=0.05)


def test_log_char_poly_matches_polynomial():
    phi = haar_angles(5, 1, np.random.default_rng(5))[0]
    theta = np.array([0.3, 2.0, -1.1])
    ref = np.log(np.abs(np.polyval(np.poly(np.exp(1j * phi)), np.exp(1j * theta))))
    assert np.allclose(log_char_poly(phi, theta), ref, atol=1e-12)


def test_cue_maxima_against_dense_scan():
    phi = haar_angles(9, 4, np.random.default_rng(8))
    theta = np.linspace(0, TWO_PI, 400_001)
    ref = [log_char_poly(p, theta).max() for p in phi]
    assert np.allclose(cue_maxima(9, 4, np.random.default_rng(8)), ref, atol=1e-9)


def test_cue_maxima_large_n():
    phi = haar_angles(120, 2, np.random.default_rng(9))
    theta = np.linspace(0, TWO_PI, 2_000_001)
    ref = [log_char_poly(p, theta).max() for p in phi]
    assert np.allclose(cue_maxima(120, 2, np.random.default_rng(9)), ref, atol=1e-6)


def test_cue_one_by_one():
    assert np.allclose(cue_maxima(1, 5, np.random.default_rng(0)), math.log(2), atol=1e-12)


def test_window_prime_sum_is_the_window_average():
    lo = np.array([1000.0, 98765.4])
    t = lo[:, None] + np.linspace(0, TWO_PI, 200_001)
    p = np.array([2, 3, 5, 7, 11, 13])
    f = (np.cos(t[..., None] * np.log(p)) / np.sqrt(p)).sum(axis=-1)
    ref = np.trapezoid(f, t, axis=1) / TWO_PI
    assert np.allclose(window_prime_sum(lo, P=13), ref, atol=1e-9)


def test_prime_shift_variance_over_heights():
    lo = np.random.default_rng(4).uniform(1e6, 1e9, 100_000)
    v = window_prime_sum(lo, P=1000).var()
    assert v == pytest.approx(prime_shift_variance(P=1000), rel=0.03)
    assert prime_shift_variance(P=1000) < prime_shift_variance() < 0.0431


def test_prime_shift_sample_variance():
    x = prime_shift_sample(200_000, np.random.default_rng(6), P=200)
    assert x.var() == pytest.approx(prime_shift_variance(P=200), rel=0.02)
    assert abs(x.mean()) < 0.002


def test_ks_distance():
    x = np.arange(10.0)
    assert ks_distance(x, x) == 0
    assert ks_distance(x, x + 0.5) == pytest.approx(0.1)
    assert ks_distance(x, x + 100) == 1
