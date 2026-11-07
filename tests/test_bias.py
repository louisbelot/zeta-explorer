import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.bias import (
    bessel_j0,
    l_log_derivative_one,
    square_roots,
    zero_sum,
)
from zeta_explorer.characters import conrey, quadratic
from zeta_explorer.lfunction import l_zeros


def test_j0_matches_mpmath():
    z = np.random.default_rng(3).uniform(0, 80, 400)
    ref = np.array([float(mpmath.besselj(0, v)) for v in z])
    assert np.max(np.abs(bessel_j0(z) - ref)) < 1e-14


def test_square_roots():
    assert [square_roots(8, a) for a in (1, 3, 5, 7)] == [4, 0, 0, 0]
    assert [square_roots(7, a) for a in range(1, 7)] == [2, 2, 0, 2, 0, 0]


def test_log_derivative_closed_form():
    gamma = float(mpmath.euler)
    ref = gamma + math.log(2 * math.pi) + 2 * math.log(math.gamma(0.75) / math.gamma(0.25))
    assert abs(l_log_derivative_one(quadratic(-4)) - ref) < 1e-13


@pytest.mark.parametrize(
    "d, ref",
    [(-3, 0.056615), (-4, 0.077784), (-8, 0.158037), (12, 0.165083)],
)
def test_zero_sum_matches_feuerverger_martin(d, ref):
    assert abs(zero_sum(quadratic(d)) / 2 - ref) < 1e-6


def test_zero_sum_tail_matches_density_of_zeros():
    chi, T = quadratic(-4), 300.0
    g = l_zeros(chi, 0, T)
    tail = zero_sum(chi) / 2 - np.sum(1 / (0.25 + g**2))
    expected = (math.log(4 * T / (2 * math.pi)) + 1) / (2 * math.pi * T)
    assert abs(tail / expected - 1) < 0.02


def test_complex_zero_sum_is_real_and_conjugate_invariant():
    chi = conrey(5, 2)
    assert not chi.is_real
    assert abs(zero_sum(chi) - zero_sum(chi.conj())) < 1e-12