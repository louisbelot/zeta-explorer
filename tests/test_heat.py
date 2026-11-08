import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.heat import (
    collision,
    ht,
    ht_prime,
    ht_second,
    newton,
    phi,
    real_zeros,
    zero_paths,
)
from zeta_explorer.zeta import xi

GAMMA = [14.134725141734693, 21.022039638771555, 25.010857580145688, 30.424876125859513,
         32.935061587739189]


@pytest.mark.parametrize("z", [0.0, 1.0, 13.7, 100.0, 400.0, 30 + 0.7j, 50 - 2j, 7.5 + 3j])
def test_h0_is_xi(z):
    assert ht(z, 0.0) == pytest.approx(complex(xi(0.5 + 0.5j * z)) / 8, rel=2e-12)


def test_h0_at_the_origin():
    assert ht(0.0, 0.0).real == pytest.approx(0.497120778188314 / 8, rel=1e-14)


def defining_integral(z, t):
    mpmath.mp.dps = 30

    def f(u):
        e = mpmath.exp(4 * u)
        s = sum((2 * mpmath.pi**2 * n**4 * mpmath.exp(9 * u) - 3 * mpmath.pi * n**2
                 * mpmath.exp(5 * u)) * mpmath.exp(-mpmath.pi * n**2 * e) for n in range(1, 25))
        return mpmath.exp(t * u * u) * s * mpmath.cos(z * u)

    return float(mpmath.quad(f, [0, 0.25, 0.5, 1, 2]))


@pytest.mark.parametrize("z, t", [(3.0, -0.3), (25.0, 0.25), (12.0, 0.0), (40.0, -1.0)])
def test_against_the_defining_integral(z, t):
    assert ht(z, t).real == pytest.approx(defining_integral(z, t), rel=1e-10)


def test_phi_is_even():
    u = np.linspace(0.01, 0.6, 7)
    assert phi(u) == pytest.approx(phi(-u), rel=1e-12)


def test_backward_heat_equation():
    z, t, e, h = 20.5, 0.1, 1e-4, 1e-2
    dt = (ht(z, t + e) - ht(z, t - e)) / (2 * e)
    dzz = (ht(z + h, t) - 2 * ht(z, t) + ht(z - h, t)) / h**2
    assert dt.real == pytest.approx(-dzz.real, rel=1e-5)


def test_symmetries():
    z, t = 17.3 + 1.1j, -0.15
    assert ht(-z, t) == pytest.approx(ht(z, t), rel=1e-13)
    assert ht(z.conjugate(), t) == pytest.approx(ht(z, t).conjugate(), rel=1e-13)
    assert ht_prime(-z, t) == pytest.approx(-ht_prime(z, t), rel=1e-13)
    x = np.linspace(-60, 60, 41)
    v = np.asarray(ht(x, 0.2))
    assert np.all(np.abs(v.imag) <= 1e-12 * np.abs(v).max())


def test_derivatives():
    z, t, h = 33.0 + 0.4j, 0.05, 1e-5
    fd = (ht(z + h, t) - ht(z - h, t)) / (2 * h)
    assert ht_prime(z, t) == pytest.approx(fd, rel=1e-7)
    fd2 = (ht_prime(z + h, t) - ht_prime(z - h, t)) / (2 * h)
    assert ht_second(z, t) == pytest.approx(fd2, rel=1e-7)


def test_zeros_at_time_zero_are_twice_the_zeta_zeros():
    assert real_zeros(0.0, 20.0, 70.0) == pytest.approx(2 * np.array(GAMMA), abs=1e-9)


def test_collision_of_a_close_pair_is_a_double_zero():
    tc, xc = collision(1128.2, 1129.2, -0.1, -0.05)
    assert -0.1 < tc < -0.05
    scale = abs(ht(1128.0, tc))
    assert abs(ht(xc, tc)) < 1e-9 * scale
    assert abs(ht_prime(xc, tc)) < 1e-6 * abs(ht_prime(1128.0, tc))


def test_square_root_law_at_the_collision():
    tc, xc = collision(1128.2, 1129.2, -0.1, -0.05)
    dt = 1e-4
    a, b = real_zeros(tc + dt, 1128.2, 1129.2, 0.002)
    assert (b - a) ** 2 / (8 * dt) == pytest.approx(1, rel=0.02)
    z = newton(xc + 1j * math.sqrt(2 * dt), tc - dt)
    assert z.imag**2 / (2 * dt) == pytest.approx(1, rel=0.02)
    assert abs(ht(z, tc - dt)) < 1e-9 * abs(ht(1128.0, tc - dt))


def test_zero_paths_follow_the_pair_off_the_line():
    before, after = zero_paths([-0.05, -0.07], 1127.5, 1130.0)
    assert np.all(before.imag == 0) and len(before) == 2
    assert len(after) == 1 and after[0].imag > 0
    assert abs(ht(after[0], -0.07)) < 1e-9 * abs(ht(1127.5, -0.07))


def test_zero_paths_wants_decreasing_time():
    with pytest.raises(ValueError):
        zero_paths([0.0, 0.1], 20.0, 30.0)
