from fractions import Fraction

import mpmath
import numpy as np
import pytest

from zeta_explorer.jensen import (
    is_hyperbolic,
    jensen,
    power_sums,
    taylor_gamma,
    theta_moments,
)


@pytest.fixture(scope="module")
def gamma():
    return taylor_gamma(range(232), dps=80)


def test_theta_moments_match_quadrature():
    ms = [0, 3, 40]
    F = theta_moments(ms, dps=40)
    with mpmath.workdps(50):
        theta = lambda t: mpmath.jtheta(3, 0, mpmath.exp(-mpmath.pi * t)) / 2 - 0.5  # noqa: E731
        pts = [mpmath.mpf(k) / 8 for k in range(81)]
        for m, f in zip(ms, F, strict=True):
            ref = mpmath.quad(lambda u, m=m: u**m * mpmath.exp(u / 4) * theta(mpmath.exp(u)), pts)
            assert abs(f / ref - 1) < 1e-35


def test_gamma_against_gorz_table():
    g = taylor_gamma([10, 100, 1000], dps=30)
    table = ["1.6323380490e-17", "6.5777263785e-205", "3.8760340890e-2567"]
    for x, ref in zip(g, table, strict=True):
        assert abs(x / mpmath.mpf(ref) - 1) < 1e-10


def test_gamma_against_taylor_series():
    g = taylor_gamma(range(4), dps=30)
    with mpmath.workdps(30):
        half = mpmath.mpf(1) / 2
        f = lambda z: (4 * z * z - 1) * mpmath.pi ** (-(half + z) / 2) * mpmath.gamma(  # noqa: E731
            (half + z) / 2) * mpmath.zeta(half + z)
        c = mpmath.taylor(f, 0, 6)
        for n in range(4):
            assert abs(g[n] / (c[2 * n] * mpmath.factorial(n)) - 1) < 1e-20


def test_power_sums():
    a = [Fraction(-6), Fraction(11), Fraction(-6), Fraction(1)]
    assert power_sums(a, 5) == [3, 6, 14, 36, 98]


def test_hyperbolic_simple_cases():
    assert is_hyperbolic([Fraction(-6), Fraction(11), Fraction(-6), Fraction(1)])
    assert not is_hyperbolic([Fraction(1), Fraction(0), Fraction(1)])
    assert not is_hyperbolic([Fraction(1), Fraction(-2), Fraction(1)])
    eps = Fraction(1, 10**30)
    assert is_hyperbolic([1 + eps, -2 - eps, Fraction(1)])
    assert not is_hyperbolic([1 + eps, Fraction(-2), Fraction(1)])


def test_hyperbolic_random_polynomials():
    rng = np.random.default_rng(7)
    for _ in range(20):
        roots = rng.normal(size=rng.integers(2, 9))
        a = [Fraction(1)]
        for r in roots:
            a = [Fraction(0)] + a
            for i in range(len(a) - 1):
                a[i] -= Fraction(float(r)) * a[i + 1]
        assert is_hyperbolic(a)
        q = [c * Fraction(1) for c in a] + [Fraction(0), Fraction(0)]
        for i in range(len(a)):
            q[i + 2] += a[i]
            q[i] += Fraction(1, 4) * a[i]
        assert not is_hyperbolic(q)


def test_degree_two_is_turan(gamma):
    for n in range(60):
        turan = gamma[n + 1] ** 2 > gamma[n] * gamma[n + 2]
        assert turan and is_hyperbolic(jensen(gamma, 2, n))


def test_small_degrees_hyperbolic(gamma):
    for d in range(1, 9):
        for n in range(0, 40, 3):
            assert is_hyperbolic(jensen(gamma, d, n))


def test_degree_twelve_far_out(gamma):
    assert is_hyperbolic(jensen(gamma, 12, 200))
