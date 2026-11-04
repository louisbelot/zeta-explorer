import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.characters import is_fundamental, kronecker
from zeta_explorer.explicit import primes
from zeta_explorer.lowlying import w_symplectic
from zeta_explorer.ratios import arithmetic, quadratic_density, zeta_logderiv
from zeta_explorer.zeta import digamma


def a_d(alpha, gamma, p):
    p = p.astype(complex)
    local = (1 - p ** (-1 - alpha - gamma)) ** -1 * (
        1 - 1 / ((p + 1) * p ** (1 + 2 * alpha)) - 1 / ((p + 1) * p ** (alpha + gamma))
    )
    return np.prod(local)


@pytest.mark.parametrize("s", [1 + 0.01j, 1 + 0.7j, 1 + 3j, 1 - 12j, 0.5 + 14.2j, 2 + 0.1j])
def test_zeta_logderiv_matches_mpmath(s):
    expected = complex(mpmath.zeta(s, derivative=1) / mpmath.zeta(s))
    assert abs(zeta_logderiv(s) - expected) < 1e-13 * abs(expected)


def test_zeta_logderiv_laurent():
    u = 1e-6j
    assert zeta_logderiv(1 + u) + 1 / u == pytest.approx(np.euler_gamma, abs=1e-6)
    with pytest.raises(ValueError):
        zeta_logderiv(1.0)


def test_arithmetic_at_zero():
    a, a_prime = arithmetic(0.0, 10**6)
    with mpmath.workdps(25):
        expected = float(-mpmath.zeta(2, derivative=1) / mpmath.zeta(2))
    assert a[0] == 1
    assert a_prime[0].real == pytest.approx(expected, abs=2e-5)


@pytest.mark.parametrize("t", [0.4, 1.3, 7.0])
def test_arithmetic_matches_conrey_snaith_product(t):
    p = primes(3000)
    a, a_prime = arithmetic(t, 3000)
    r = 1j * t
    assert a[0] == pytest.approx(a_d(-r, r, p), rel=1e-12)
    h = 1e-5
    deriv = (a_d(r + h, r, p) - a_d(r - h, r, p)) / (2 * h)
    assert a_d(r, r, p) == pytest.approx(1, abs=1e-14)
    assert a_prime[0] == pytest.approx(deriv, abs=1e-8)


def test_family_matches_local_densities():
    ds = [d for m in range(1, 30001) for d in (-m, m) if is_fundamental(d)]
    for p in (2, 3, 5, 7):
        signs = np.array([kronecker(d, p) for d in ds])
        assert np.mean(signs == 0) == pytest.approx(1 / (p + 1), abs=0.01)
        assert abs(np.mean(signs[signs != 0])) < 0.02


def test_density_vanishes_at_the_central_point():
    x = np.array([1e-4, 1e-3, 0.5])
    rho = quadratic_density([-4003, 5001], x, p_max=10**4, n=800)
    assert rho[0] < 1e-6
    assert rho[1] < 1e-4
    assert rho[2] == pytest.approx(1, abs=0.6)


def test_leading_correction_is_closed_form():
    with mpmath.workdps(25):
        c = np.euler_gamma - float(mpmath.zeta(2, derivative=1) / mpmath.zeta(2))
    x = np.linspace(0.01, 3, 300)
    lead = 2 * c * (1 - np.cos(2 * np.pi * x))
    worst = []
    for d in (1e100, 1e200, -1e200):
        rho = quadratic_density([d], x, p_max=10**4, n=800)
        L = math.log(abs(d) / math.pi) + digamma(0.25 + (d < 0) / 2).real
        worst.append(np.max(np.abs((rho - w_symplectic(x)) * L - lead)))
    assert worst[1] < 0.05
    assert worst[2] < 0.05
    assert worst[1] < worst[0] / 2


def test_density_rejects_small_conductors():
    with pytest.raises(ValueError):
        quadratic_density([105], [0.5])
