import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.maxima import (
    TWO_PI,
    bessel_k,
    fhk_cdf,
    fhk_density,
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
