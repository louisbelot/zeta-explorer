import math

import numpy as np
import pytest

from zeta_explorer.li import li_coefficients, li_exact, li_tail, li_trend
from zeta_explorer.moments import EULER_GAMMA
from zeta_explorer.spacing import first_zeros


@pytest.fixture(scope="module")
def gamma():
    return first_zeros(1000)


def test_lambda_1_closed_form():
    ref = 1 + EULER_GAMMA / 2 - math.log(4 * math.pi) / 2
    assert li_exact(1)[0] == pytest.approx(ref, abs=1e-15)


def test_exact_is_stable_in_the_contour():
    a = li_exact(25)
    b = li_exact(25, nodes=220, r=0.4)
    assert np.max(np.abs(a - b) / a) < 1e-12


def test_zero_sum_matches_exact(gamma):
    n = np.arange(1, 21)
    assert np.max(np.abs(li_coefficients(n, gamma) - li_exact(20))) < 1e-6


def test_tail_replaces_zeros(gamma):
    n = np.array([1, 10, 100, 400])
    full = li_coefficients(n, gamma)
    half = li_coefficients(n, gamma[:600])
    assert np.max(np.abs(full - half) / full) < 1e-5


def test_tail_is_the_smooth_count():
    T = 5000.0
    n = np.array([3.0])
    t = np.linspace(T, 4e6, 4_000_001)
    f = 4 * np.sin(n[0] * np.arctan(0.5 / t)) ** 2 * np.log(t / (2 * math.pi)) / (2 * math.pi)
    approx = np.trapezoid(f, t) + 9 * (np.log(4e6 / (2 * math.pi)) + 1) / (2 * math.pi * 4e6)
    assert li_tail(n, T)[0] == pytest.approx(approx, rel=1e-6)


def test_positive_increasing_near_trend(gamma):
    n = np.arange(1, 401)
    lam = li_coefficients(n, gamma)
    assert np.all(lam > 0)
    assert np.all(np.diff(lam) > 0)
    assert lam[-1] == pytest.approx(li_trend(400), rel=0.05)
