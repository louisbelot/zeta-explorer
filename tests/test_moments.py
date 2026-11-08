import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.moments import (
    EULER_GAMMA,
    arithmetic_factor,
    gregory_weights,
    leading_coefficient,
    moment_polynomial,
    polynomial_integral,
    smoothed_polynomial_integral,
    z_power_integrals,
)


def _mp_block(C, D, ks, nodes=40):
    x, w = np.polynomial.legendre.leggauss(nodes)
    total = np.zeros(len(ks))
    for a in np.arange(C, D, 1.0):
        z = np.array([float(mpmath.siegelz(t)) for t in a + (x + 1) / 2])
        total += [0.5 * np.dot(w, z ** (2 * k)) for k in ks]
    return total


def test_leading_coefficients():
    assert leading_coefficient(1) == pytest.approx(1, rel=1e-12)
    assert leading_coefficient(2) == pytest.approx(1 / (2 * math.pi**2), rel=1e-10)
    assert leading_coefficient(3) == pytest.approx(5.708527034652788e-6, rel=1e-9)
    assert moment_polynomial(2)[0] == pytest.approx(leading_coefficient(2), rel=1e-10)


def test_arithmetic_factor_closed_form():
    assert arithmetic_factor(2) == pytest.approx(6 / math.pi**2, rel=1e-10)


def test_p1_is_ingham():
    T = 12345.0
    exact = T * math.log(T / (2 * math.pi)) + (2 * EULER_GAMMA - 1) * T
    assert polynomial_integral(1, 0, T) == pytest.approx(exact, rel=1e-14)


def test_polynomial_blocks_match_cfkrs():
    assert polynomial_integral(3, 0, 50000) == pytest.approx(7236872972.7, rel=1e-11)
    assert polynomial_integral(3, 50000, 100000) == pytest.approx(15696470555.3, rel=1e-11)


def test_smoothed_polynomial_matches_cfkrs():
    assert smoothed_polynomial_integral(1, 1e4) == pytest.approx(79496.7897047, rel=1e-11)
    assert smoothed_polynomial_integral(2, 1e4) == pytest.approx(5088332.55512, rel=1e-11)
    assert smoothed_polynomial_integral(3, 1e4) == pytest.approx(708965694.5, rel=1e-9)


def test_smoothed_polynomial_against_quadrature():
    T = 300.0
    p = list(moment_polynomial(3)[::-1])

    def f(t):
        return mpmath.polyval(p, mpmath.log(t / (2 * mpmath.pi)), asc=True) * mpmath.exp(-t / T)

    ref = float(mpmath.quad(f, [0, 2 * mpmath.pi, T, 10 * T, mpmath.inf]))
    assert smoothed_polynomial_integral(3, T) == pytest.approx(ref, rel=1e-12)


def test_gregory_weights():
    assert gregory_weights(3) == pytest.approx((3 / 8, 7 / 6, 23 / 24), abs=1e-15)
    m, n = 8, 40
    w = np.ones(n + 1)
    w[:m] = gregory_weights(m)
    w[-m:] = gregory_weights(m)[::-1]
    x = np.arange(n + 1) / n
    for d in range(m):
        assert np.dot(w, x**d) / n == pytest.approx(1 / (d + 1), rel=1e-12)


def test_block_against_mpmath():
    ref = _mp_block(77402.0, 77405.0, (3,))
    assert z_power_integrals(77402.0, 77405.0, ks=(3,), h=0.005) == pytest.approx(ref, rel=1e-10)


def test_block_across_the_seam():
    got = z_power_integrals(998.0, 1002.0, ks=(1, 2), h=0.005)
    assert got == pytest.approx(_mp_block(998.0, 1002.0, (1, 2)), rel=1e-10)


def test_windows_do_not_matter():
    a = z_power_integrals(30000.0, 31000.0, ks=(2,), window=150.0)
    b = z_power_integrals(30000.0, 31000.0, ks=(2,), window=5000.0)
    assert a[0] == pytest.approx(b[0], rel=1e-12)


def test_second_moment_follows_p1():
    C, D = 1e5, 1.2e5
    assert z_power_integrals(C, D, ks=(1,))[0] / polynomial_integral(1, C, D) == pytest.approx(
        1, abs=2e-3
    )
