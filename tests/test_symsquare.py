import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.elliptic import Curve, EllipticL
from zeta_explorer.symsquare import MU, SymSquare, coefficients, smoothed
from zeta_explorer.zeta import zeta

E11 = Curve(0, -1, 1, 0, 0)
E37 = Curve(0, 0, 1, -1, 0)


def test_coefficients_are_multiplicative_from_ap():
    c = coefficients(E11, 128)
    a = EllipticL(E11).an(128)
    for p in (2, 3, 5, 7, 13):
        assert c[p] == pytest.approx(a[p] ** 2 / p - 1, abs=1e-14)
    assert c[11] == pytest.approx(1 / 11, abs=1e-15)
    assert c[121] == pytest.approx(1 / 121, abs=1e-15)
    assert c[6] == pytest.approx(c[2] * c[3], abs=1e-14)
    assert c[44] == pytest.approx(c[4] * c[11], abs=1e-14)


def test_agrees_with_euler_product_right_of_one():
    L = SymSquare(E11)
    assert L(2.0) == pytest.approx(L.euler(2.0), abs=1e-8)
    assert L(2.5 - 4j) == pytest.approx(L.euler(2.5 - 4j), abs=1e-7)


@pytest.mark.parametrize("s", [0.5 + 10j, 1 + 2j, 1.0, 0.3 + 5j, 1 + 30j])
def test_weight_does_not_matter(s):
    L = SymSquare(E11)
    a = L._a(s)
    ref = smoothed(s, a, 121, MU)
    for b, h in ((14, 0.1), (60, 0.15)):
        assert abs(smoothed(s, a, 121, MU, b=b, h=h) - ref) < 1e-9


@pytest.mark.parametrize("mu, q", [((1.0, 1.0, 2.0), 11.0), ((0.0, 1.0, 1.0), 121.0),
                                   ((1.0, 2.0, 3.0), 121.0)])
def test_wrong_gamma_factor_breaks_it(mu, q):
    s = 1 + 2j
    a = SymSquare(E11)._a(s)
    assert abs(smoothed(s, a, q, mu) - smoothed(s, a, q, mu, b=14, h=0.1)) > 1e-3


def test_real_on_the_real_axis_and_self_dual():
    L = SymSquare(E11)
    assert abs(L(1.0).imag) < 1e-14
    assert L(1 - 3j) == pytest.approx(np.conj(L(1 + 3j)), abs=1e-12)
    z = L(np.array([0.5 + 7j, 0.5 + 7.5j]))
    assert z.shape == (2,)


def test_rankin_selberg_mean_of_lambda_squared():
    X = 500.0
    n = np.arange(1, 40 * int(X))
    lam2 = EllipticL(E11).an(len(n))[1:].astype(float) ** 2 / n
    L = SymSquare(E11)
    main = X * L(1.0).real * 6 / math.pi**2 * 11 / 12 + L(0.0).real / 2
    w = 1j * math.pi * (2 * np.arange(-8, 8) + 1) / math.log(11)
    poles = np.exp(np.array([complex(mpmath.loggamma(z)) for z in w]) + w * math.log(X))
    poles *= zeta(w) * L(w) / zeta(2 * w) / math.log(11)
    total = main + poles.sum().real
    assert abs(poles.sum()) > 0.1
    assert np.sum(lam2 * np.exp(-n / X)) == pytest.approx(total, rel=2e-6)


@pytest.mark.parametrize("s", [1 + 4j, 0.5 + 10j, 1.0])
def test_derivative_against_cauchy_integral(s):
    L = SymSquare(E37)
    value, dvalue = L.with_derivative(s)
    r = 0.05
    ring = np.exp(2j * math.pi * (np.arange(32) + 0.5) / 32)
    assert dvalue[0] == pytest.approx(np.mean(L(s + r * ring) / ring) / r, abs=1e-11)
    assert value[0] == pytest.approx(L(s), abs=1e-14)


def test_needs_semistable_curve():
    with pytest.raises(ValueError):
        SymSquare(Curve(0, 0, 0, -1, 0))
