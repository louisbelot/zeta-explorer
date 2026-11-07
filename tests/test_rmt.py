import math

import numpy as np
import pytest

from zeta_explorer.rmt import (
    even_orthogonal,
    even_orthogonal_n,
    first_density,
    fit_n,
    fredholm_det,
    gap,
    sine,
    symplectic,
    symplectic_n,
)

S = np.array([0.1, 0.4, 0.9])


def weyl_gap(alpha, density, nodes=80):
    u, w = np.polynomial.legendre.leggauss(nodes)

    def integral(lo):
        a = (np.pi - lo) / 2 * u + (np.pi + lo) / 2
        wa = w * (np.pi - lo) / 2
        t1, t2 = np.meshgrid(a, a)
        return np.sum(np.outer(wa, wa) * density(t1, t2))

    return integral(alpha) / integral(0.0)


def test_usp2_closed_form():
    k = symplectic_n(1)
    assert gap(k, S) == pytest.approx(1 - S + np.sin(2 * np.pi * S) / (2 * np.pi), abs=1e-14)
    assert first_density(k, S) == pytest.approx(1 - np.cos(2 * np.pi * S), abs=1e-14)


def test_so2_is_uniform():
    assert gap(even_orthogonal_n(1), S) == pytest.approx(1 - S, abs=1e-14)


@pytest.mark.parametrize("s", [0.3, 0.8, 1.5])
def test_rank_two_against_weyl_integral(s):
    alpha = np.pi * s / 2
    usp = weyl_gap(alpha, lambda a, b: (np.cos(a) - np.cos(b)) ** 2 * (np.sin(a) * np.sin(b)) ** 2)
    so = weyl_gap(alpha, lambda a, b: (np.cos(a) - np.cos(b)) ** 2)
    assert gap(symplectic_n(2), s) == pytest.approx(usp, abs=1e-13)
    assert gap(even_orthogonal_n(2), s) == pytest.approx(so, abs=1e-13)


def test_sine_gap_series():
    s = 0.1
    series = 1 - s + np.pi**2 * s**4 / 36 - np.pi**4 * s**6 / 675 + np.pi**6 * s**8 / 17640
    assert fredholm_det(sine, -s / 2, s / 2) == pytest.approx(series, abs=1e-11)


@pytest.mark.parametrize("s", [0.2, 0.7, 1.6])
def test_parity_split_of_the_sine_kernel(s):
    whole = fredholm_det(sine, -s, s)
    assert gap(symplectic, s) * gap(even_orthogonal, s) == pytest.approx(whole, abs=1e-14)


def test_first_zero_near_the_origin():
    assert first_density(symplectic, 1e-3) / 1e-6 == pytest.approx(2 * np.pi**2 / 3, rel=1e-5)
    assert first_density(even_orthogonal, 0.0) == 2


def test_first_density_is_a_density():
    s = np.linspace(0, 4, 2001)
    p = first_density(symplectic, s)
    assert np.all(p >= 0)
    assert np.trapezoid(p, s) == pytest.approx(1, abs=1e-6)
    assert gap(symplectic, 4.0) < 1e-12


def test_density_is_minus_gap_derivative():
    h = 1e-5
    for k in (symplectic, symplectic_n(3), even_orthogonal):
        fd = -(gap(k, 0.8 + h) - gap(k, 0.8 - h)) / (2 * h)
        assert first_density(k, 0.8) == pytest.approx(fd, abs=1e-8)


def test_finite_n_converges_like_one_over_n():
    limit = gap(symplectic, 1.0)
    err = [gap(symplectic_n(n), 1.0) - limit for n in (100, 1000)]
    assert err[0] / err[1] == pytest.approx(10, rel=0.02)


def test_quadrature_has_converged():
    assert gap(symplectic, 2.0, m=20) == pytest.approx(gap(symplectic, 2.0, m=60), abs=1e-14)


def test_fit_recovers_matrix_size():
    rng = np.random.default_rng(7)
    n = 3
    s = np.linspace(0, 1.6, 4001)
    cdf = 1 - gap(symplectic_n(n), s)
    sample = np.interp(rng.uniform(0, cdf[-1], 4000), cdf, s)
    best, lo, hi, _ = fit_n(sample, np.linspace(1.5, 8, 27))
    assert lo <= n <= hi
    assert hi - lo < 1.5
    assert math.isclose(best, n, abs_tol=0.6)


def test_kernels_reject_tiny_n():
    with pytest.raises(ValueError):
        symplectic_n(0.5)
