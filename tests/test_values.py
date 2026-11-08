import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.explicit import primes
from zeta_explorer.moments import prime_square_tail
from zeta_explorer.spacing import first_zeros, zeros_between
from zeta_explorer.values import (
    arithmetic_cumulants,
    cue_cumulants,
    cue_density,
    hybrid_density,
    log_z_sample,
    records,
    sample_cumulants,
    z_extrema,
)


def test_extrema_are_critical_points_of_z():
    t, z = z_extrema(first_zeros(30))
    assert len(t) == 29
    slopes = [float(mpmath.siegelz(x, derivative=1)) for x in t]
    assert np.max(np.abs(slopes)) < 1e-6
    assert np.all(np.sign(z[:-1]) != np.sign(z[1:]))


def test_extrema_above_the_mpmath_seam():
    gamma = zeros_between(5000, 5040)
    t, z = z_extrema(gamma)
    assert np.all((gamma[:-1] < t) & (t < gamma[1:]))
    for i in (0, 17, len(t) - 1):
        assert z[i] == pytest.approx(float(mpmath.siegelz(t[i])), abs=1e-9)
        assert abs(float(mpmath.siegelz(t[i], derivative=1))) < 1e-5


def test_records():
    assert records([1, -3, 2, 3, -5, 0.5, 7]).tolist() == [0, 1, 4, 6]


def test_cue_cumulants_match_polygamma():
    N = 12
    k = cue_cumulants(N)
    for order, scale in ((1, 0.5), (2, 0.75), (3, 0.875)):
        ref = scale * float(mpmath.fsum(mpmath.psi(order, j) for j in range(1, N + 1)))
        assert k[order] == pytest.approx(ref, rel=1e-13)
    assert k[0] == 0


def test_cue_variance_large_n():
    N = 2000
    assert cue_cumulants(N)[1] == pytest.approx(0.5 * (math.log(N) + 0.5772156649 + 1), abs=1e-3)


def test_cue_density_moments():
    x, dx = np.linspace(-30, 8, 7601, retstep=True)
    p = cue_density(x, 12)
    k = cue_cumulants(12)
    assert p.sum() * dx == pytest.approx(1, abs=1e-9)
    assert (x * p).sum() * dx == pytest.approx(0, abs=1e-8)
    assert (x**2 * p).sum() * dx == pytest.approx(k[1], rel=1e-8)
    assert (x**3 * p).sum() * dx == pytest.approx(k[2], rel=1e-7)


def _haar_unitary(n, rng):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    return Q * (d / np.abs(d))


def test_random_unitaries_match_cue_density():
    rng = np.random.default_rng(5)
    N = 12
    eye = np.eye(N)
    x = np.array([np.linalg.slogdet(eye - _haar_unitary(N, rng))[1] for _ in range(20_000)])
    k = cue_cumulants(N)
    assert sample_cumulants(x)[:3] == pytest.approx(k[:3], abs=0.1)
    hist, edges = np.histogram(x, bins=24, range=(-6, 3), density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    assert np.max(np.abs(hist - cue_density(centers, N))) < 0.015


def test_log_z_beats_selberg_with_cue():
    T = 2 * math.pi * math.exp(10)
    x = log_z_sample(T, 50_000, np.random.default_rng(0))
    k = sample_cumulants(x)
    assert abs(k[0]) < 0.03
    cue = cue_cumulants(10)[1]
    selberg = 0.5 * math.log(math.log(T))
    assert abs(k[1] - cue) < 0.15
    assert abs(k[1] - cue) < abs(k[1] - selberg) / 4


def _x_p_cumulants(p, nodes=512):
    theta = 2 * math.pi * np.arange(nodes) / nodes
    x = -np.log(np.abs(1 - np.exp(1j * theta) / math.sqrt(p)))
    c = x - x.mean()
    k2 = np.mean(c**2)
    return x.mean(), k2, np.mean(c**3), np.mean(c**4) - 3 * k2**2


def test_arithmetic_variance_is_dilog_sum():
    ps = primes(1000)
    ref = sum(0.5 * float(mpmath.polylog(2, 1 / int(p))) + 0.5 * math.log1p(-1 / p) for p in ps)
    got = arithmetic_cumulants(1000)
    assert got[1] == pytest.approx(ref - prime_square_tail(1000) / 8, rel=1e-12)
    assert got[0] == 0


def test_arithmetic_higher_cumulants_by_quadrature():
    ks = np.array([_x_p_cumulants(p) for p in primes(500)])
    assert np.all(np.abs(ks[:, 0]) < 1e-15)
    got = arithmetic_cumulants(500)
    tail = prime_square_tail(500)
    assert got[2] == pytest.approx(ks[:, 2].sum() + 3 * tail / 8, rel=1e-10)
    assert got[3] == pytest.approx(ks[:, 3].sum() - 3 * tail / 8, rel=1e-10)


def test_hybrid_density_moments():
    x, dx = np.linspace(-30, 8, 3801, retstep=True)
    p = hybrid_density(x, 12)
    k = cue_cumulants(12) + arithmetic_cumulants()
    assert p.sum() * dx == pytest.approx(1, abs=1e-8)
    assert (x * p).sum() * dx == pytest.approx(0, abs=1e-7)
    assert (x**2 * p).sum() * dx == pytest.approx(k[1], rel=1e-6)
    assert (x**3 * p).sum() * dx == pytest.approx(k[2], rel=1e-5)
    assert p.min() > -1e-6 * p.max()
