import math
from fractions import Fraction

import mpmath
import numpy as np
import pytest

from zeta_explorer.delta import lam, tau
from zeta_explorer.selfdual import smoothed
from zeta_explorer.sympower import (
    ZAGIER,
    coefficients,
    critical_ratio,
    delta_power,
    euler,
    gamma_shifts,
    petersson,
    rational,
)

YOSHIDA_SYM3 = [
    4.1558656464, 5.5491219562, 8.1117756122, 10.8952834492, 12.0523651120, 13.4542992617,
    14.9275108496, 16.3036898019, 17.7350625418, 18.837088412,
]
YOSHIDA_SYM4 = [2.3864500, 4.3752457, 6.0435487, 7.571907, 8.841633, 10.605890, 11.437474]


def test_gamma_shifts():
    assert gamma_shifts(1) == (5.5, 6.5)
    assert gamma_shifts(2) == (1.0, 11.0, 12.0)
    assert gamma_shifts(3) == (5.5, 6.5, 16.5, 17.5)
    assert gamma_shifts(4) == (0.0, 11.0, 12.0, 22.0, 23.0)


def test_prime_coefficients_are_chebyshev_in_lambda():
    x = lam(64)
    for p in (2, 3, 5, 7, 61):
        assert coefficients(2, 64)[p] == pytest.approx(x[p] ** 2 - 1, abs=1e-13)
        assert coefficients(3, 64)[p] == pytest.approx(x[p] ** 3 - 2 * x[p], abs=1e-13)
        assert coefficients(4, 64)[p] == pytest.approx(
            x[p] ** 4 - 3 * x[p] ** 2 + 1, abs=1e-13
        )
    assert coefficients(1, 64) == pytest.approx(x, abs=1e-13)


def test_sym2_coefficients_from_tau_at_squares():
    N = 30
    t = tau(N * N)
    c = coefficients(2, N)
    for n in range(1, N + 1):
        want = sum(
            float(t[(n // (d * d)) ** 2]) / (n // (d * d)) ** 11
            for d in range(1, math.isqrt(n) + 1) if n % (d * d) == 0
        )
        assert c[n] == pytest.approx(want, abs=1e-12)


@pytest.mark.parametrize("m", [2, 3, 4])
@pytest.mark.parametrize("s", [0.5 + 3j, 0.5 + 20j, 0.3 + 7j, 1.0, 0.5 - 12j])
def test_weight_does_not_matter(m, s):
    L = delta_power(m)
    a = L._a(complex(s))
    ref = smoothed(complex(s), a, 1.0, L.mu, L.eps)
    for b, h in ((14, 0.1), (60, 0.15)):
        assert abs(smoothed(complex(s), a, 1.0, L.mu, L.eps, b=b, h=h) - ref) < 1e-11


@pytest.mark.parametrize("m", [2, 3, 4])
def test_wrong_root_number_or_shift_breaks_it(m):
    L = delta_power(m)
    s = 0.5 + 6j
    a = L._a(s)
    shifted = (L.mu[0] + 1,) + L.mu[1:]
    for mu, eps in ((L.mu, -L.eps), (shifted, L.eps)):
        gap = abs(smoothed(s, a, 1.0, mu, eps) - smoothed(s, a, 1.0, mu, eps, b=14, h=0.1))
        assert gap > 1e-4


@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_agrees_with_euler_product(m):
    L = delta_power(m)
    assert L(2.5 + 1j) == pytest.approx(euler(m, 2.5 + 1j, 1 << 13), abs=2e-9)


def test_petersson_norm():
    known = mpmath.mpf("1.03536205680432092234781681222516459322e-6")
    assert abs(petersson() / known - 1) < 1e-29


def test_rankin_value_at_the_edge():
    L1 = delta_power(2)(1.0).real
    want = 2**23 * math.pi**13 * float(petersson()) / math.factorial(11)
    assert L1 == pytest.approx(want, rel=1e-13)


@pytest.mark.parametrize("s", [1, 3])
def test_low_critical_values_from_the_functional_equation(s):
    assert float(critical_ratio(s)) == pytest.approx(float(ZAGIER[s]), rel=1e-12)


@pytest.mark.parametrize("s, rtol", [(5, 1e-17), (7, 1e-25), (9, 1e-27), (11, 1e-28)])
def test_critical_values_are_recognized_as_zagiers_rationals(s, rtol):
    r = critical_ratio(s)
    assert rational(r, rtol) == ZAGIER[s]
    assert abs(r * ZAGIER[s].denominator / ZAGIER[s].numerator - 1) < rtol


def test_rational_takes_the_first_convergent_that_fits():
    assert rational(0.3333333333, 1e-9) == Fraction(1, 3)
    assert rational(-22 / 7, 1e-12) == Fraction(-22, 7)
    assert rational(math.pi, 1e-6) == Fraction(355, 113)


def test_critical_ratio_rejects_other_points():
    with pytest.raises(ValueError):
        critical_ratio(2)
    with pytest.raises(ValueError):
        critical_ratio(13)


def test_sym3_zeros_against_yoshida():
    L = delta_power(3)
    assert L.eps == -1
    assert abs(L(0.5)) < 1e-13
    z = L.zeros(19.5)
    assert len(z) == len(YOSHIDA_SYM3)
    assert z[:9] == pytest.approx(YOSHIDA_SYM3[:9], abs=2e-10)
    assert z[9] == pytest.approx(YOSHIDA_SYM3[9], abs=2e-9)


def test_sym4_zeros_against_yoshida():
    L = delta_power(4)
    z = L.zeros(12.0)
    assert len(z) == len(YOSHIDA_SYM4)
    assert z[:3] == pytest.approx(YOSHIDA_SYM4[:3], abs=2e-7)
    assert z[3:] == pytest.approx(YOSHIDA_SYM4[3:], abs=2e-6)


def test_counts_follow_the_main_term():
    for m in (2, 3):
        L = delta_power(m)
        z = L.zeros(19.5)
        s = np.arange(1, len(z) + 1) - L.phi(z)
        assert np.all((s > -0.5) & (s < 1.5))


def test_only_automorphic_powers():
    with pytest.raises(ValueError):
        delta_power(5)
