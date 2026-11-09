import math

import numpy as np
import pytest

from zeta_explorer import effective as E
from zeta_explorer.heat import ht
from zeta_explorer.riemann_siegel import theta, z_rs
from zeta_explorer.spacing import first_zeros

CASES = [(0.2, 0.0), (0.2, 0.2), (0.4, 0.4), (0.05, 0.1), (0.5, 1.0)]


@pytest.mark.parametrize("t, y", CASES)
def test_kernel_matches_the_integral(t, y):
    h = complex(ht(complex(300.0, y), t)) / np.exp(complex(E.log_bt(300.0, y, t)))
    assert abs(E.kernel_ratio(300.0, y, t) - h) < 3e-14


def test_kernel_converges_in_nodes():
    a, b = E.kernel_ratio(1e4 + 0.3, 0.2, 0.2, 32), E.kernel_ratio(1e4 + 0.3, 0.2, 0.2, 64)
    assert abs(a - b) < 1e-14


def test_kernel_at_t0_is_xi():
    x = 1000.0
    h = complex(ht(x, 0.0)) / np.exp(complex(E.log_bt(x, 0.0, 0.0)))
    assert abs(E.kernel_ratio(x, 0.0, 0.0) - h) < 1e-12


@pytest.mark.parametrize("t, y", CASES)
@pytest.mark.parametrize("x", [200.0, 250.0, 1234.5, 1e4 + 0.3])
def test_theorem_bounds_hold(x, y, t):
    k = E.kernel_ratio(x, y, t)
    assert abs(k - E.ratio(x, y, t, corrected=False)) <= E.error(x, y, t, corrected=False)
    assert abs(k - E.ratio(x, y, t)) <= E.error(x, y, t)


@pytest.mark.parametrize("x", [1e4 + 0.3, 1e6 + 0.1])
def test_c_term_is_the_rs_correction(x):
    k = E.kernel_ratio(x, 0.2, 0.2)
    assert abs(k - E.ratio(x, 0.2, 0.2)) < abs(k - E.ratio(x, 0.2, 0.2, corrected=False)) / 50


def test_bounds_decay():
    x = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
    e = E.error(x, 0.2, 0.2)
    assert np.all(np.diff(e) < 0)
    assert e[-1] < 1e-5
    assert np.all(E.error_c(x, 0.2, 0.2) < E.error_c0(x, 0.2, 0.2) / 3)


def test_vectorized_matches_scalar():
    x = np.array([[250.0, 999.0], [2e4, 3e5 + 0.5]])
    r = E.ratio(x, 0.3, 0.1)
    assert r.shape == (2, 2)
    assert r[1, 0] == E.ratio(2e4, 0.3, 0.1)
    assert E.error(x, 0.3, 0.1).shape == (2, 2)


@pytest.mark.parametrize("bad", [(150.0, 0.0, 0.2), (300.0, 1.5, 0.2), (300.0, 0.0, 0.6)])
def test_outside_the_theorem(bad):
    with pytest.raises(ValueError):
        E.ratio(*bad)


def test_c0_is_bounded_and_smooth():
    p = np.linspace(-1, 1, 20001)
    v = E.c0(p)
    assert np.abs(v).max() <= 0.5
    assert np.abs(np.diff(v)).max() < 1e-3
    assert E.c0(0.5) == pytest.approx((1 - 1j) / 4)
    assert abs(E.c0(0.5 + 1e-7) - (1 - 1j) / 4) < 1e-7


def test_theorem_bounds_hold_high():
    x, y, t = 1e6 + 0.1, 0.2, 0.2
    k = E.kernel_ratio(x, y, t)
    assert abs(k - E.ratio(x, y, t)) <= E.error(x, y, t)


def test_theta_t0_is_riemann_siegel_theta():
    x = np.array([300.0, 1e4, 1e6])
    assert np.all(np.abs(E.theta_t(x, 0.0) - theta(x / 2) - math.pi) < 1 / x)


def test_theta_shift_is_t_over_16_log():
    x = np.array([1e4, 1e6, 1e8])
    for t in (0.1, 0.4):
        shift = (E.theta_t(x, t) - E.theta_t(x, 0.0)) / math.pi
        assert np.all(np.abs(shift - t / 16 * np.log(x / (4 * math.pi))) < 10 / x)


def test_theta_prime():
    x, h = 5e4, 1e-3
    for t in (0.0, 0.3):
        fd = (E.theta_t(x + h, t) - E.theta_t(x - h, t)) / (2 * h)
        assert E.theta_t_prime(x, t) == pytest.approx(fd, rel=1e-8)


def test_ideal_points():
    n = np.array([100, 10_000, 10**7])
    for t in (0.0, 0.2):
        x = E.ideal_points(n, t)
        assert np.allclose(E.theta_t(x, t), (n - 0.5) * math.pi, rtol=1e-14, atol=1e-6)


def test_zt_at_t0_is_hardy_z():
    x = np.array([2e4 + 0.37, 3e5 + 0.11, 4e6 + 0.9])
    assert np.allclose(E.zt(x, 0.0), -z_rs(x / 2), atol=2e-3)


def test_real_axis_is_real():
    x = np.linspace(1e4, 1e4 + 30, 101)
    for t in (0.0, 0.2):
        im = np.imag(np.exp(-1j * E.theta_t(x, t)) * E.ratio(x, 0.0, t))
        assert np.abs(im).max() < 1e-12


def test_zeros_at_t0_are_twice_gamma():
    g = first_zeros(200)[100:]
    x = E.real_zeros(0.0, 2 * g[0] - 0.3, 2 * g[-1] + 0.3)
    assert len(x) == len(g)
    assert np.abs(x - 2 * g).max() < 0.02
    n, d = E.displacement(2 * g, 0.0)
    off = n - np.arange(101, 201)
    assert np.all(np.abs(off) <= 1)
    assert np.mean(off == 0) > 0.8
    assert np.abs(d).max() <= 0.5


def test_displacement_of_ideal_points():
    n = np.arange(5000, 5010)
    m, d = E.displacement(E.ideal_points(n, 0.3), 0.3)
    assert np.array_equal(m, n)
    assert np.abs(d).max() < 1e-9
