import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.characters import is_fundamental, kronecker
from zeta_explorer.elliptic import Curve, EllipticL
from zeta_explorer.explicit import primes
from zeta_explorer.lowlying import w_symplectic
from zeta_explorer.ratios import (
    arithmetic,
    quadratic_density,
    twist_arithmetic,
    twist_density,
    twist_effective_n,
    twist_profile,
    twist_r1,
    twist_rho,
    zeta_logderiv,
)
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


E11 = Curve(0, -1, 1, 0, 0)


def raw_local(lam, x, y, p, terms=200):
    hecke = [1.0, lam]
    for _ in range(terms):
        hecke.append(lam * hecke[-1] - hecke[-2])
    mobius = [1.0, -lam, 1.0]
    total = sum(hecke[m] * mobius[h] * x**m * y**h
                for m in range(terms) for h in range(3) if (m + h) % 2 == 0 and m + h)
    return 1 + p / (p + 1) * total


def raw_log_g(alpha, gamma, eta, P):
    a = EllipticL(E11).an(P)
    out = 0.0
    for p in primes(P).tolist():
        x, y = p ** (-0.5 - alpha), p ** (-0.5 - gamma)
        lam = a[p] / math.sqrt(p)
        if p == 11:
            out += math.log((1 - eta * lam * y) / (1 - eta * lam * x))
        else:
            out += math.log(raw_local(lam, x, y, p))
    return out


def euler_logderivs(r, P):
    a = EllipticL(E11).an(P)
    out = 0.0
    for p in primes(P).tolist():
        X = p ** (-1 - 2 * r)
        lp = math.log(p)
        out += lp * X / (1 - X)
        lam2 = a[p] ** 2 / p
        if p == 11:
            out -= lp * lam2 * X / (1 - lam2 * X)
        else:
            e = lam2 - 1
            out -= lp * X * (e - 2 * e * X + 3 * X * X) / (1 - e * X + e * X * X - X**3)
    return out


@pytest.mark.parametrize("eta", [1, -1])
def test_twist_derivative_against_raw_sum(eta):
    r, h, P = 0.3, 1e-5, 2000
    fd = (raw_log_g(r + h, r, eta, P) - raw_log_g(r - h, r, eta, P)) / (2 * h)
    a1, a = twist_arithmetic(E11, -1j * r, eta, P)
    assert a1[0].real + euler_logderivs(r, P) == pytest.approx(fd, abs=1e-8)
    assert raw_log_g(r, r, eta, P) == pytest.approx(0, abs=1e-12)


@pytest.mark.parametrize("eta", [1, -1])
def test_twist_ratio_factor_against_raw_sum(eta):
    r, P = 0.2, 2000
    a = EllipticL(E11).an(P)
    log_y = 0.0
    for p in primes(P).tolist():
        lam2 = a[p] ** 2 / p
        q = [1 - lam2 * z if p == 11 else 1 - (lam2 - 1) * z + (lam2 - 1) * z * z - z**3
             for z in (p ** (-1 + 2 * r), 1 / p)]
        log_y += math.log((1 - p ** (-1 - 2 * r)) * q[0] / ((1 - 1 / p) * q[1]))
    _, A = twist_arithmetic(E11, -1j * r, eta, P)
    assert A[0].real == pytest.approx(math.exp(raw_log_g(-r, r, eta, P) + log_y), rel=1e-12)


@pytest.mark.parametrize("eta", [1, -1])
def test_twist_constants_cancel_at_the_center(eta):
    a1, A = twist_arithmetic(E11, np.array([0.0, 1e-4, -1e-4]), eta, 2**14)
    assert A[0] == pytest.approx(1, abs=1e-12)
    slope = (np.log(A[1]) - np.log(A[2])) / 2e-4
    assert (-1j * slope).real == pytest.approx(-2 * a1[0].real, abs=1e-7)


def test_twist_r1_matches_dhkms():
    r1 = twist_r1(E11, 1, 2**14) + np.euler_gamma + 1
    assert r1 == pytest.approx(2.8600, abs=1e-3)
    assert math.log(math.sqrt(11) * 400_000 / (2 * math.pi)) / (2 * r1) == pytest.approx(
        2.14, abs=5e-3)


def test_twist_effective_n():
    for d, eta in ((89, 1), (1201, -1)):
        L = EllipticL(E11, d)
        assert kronecker(d, 11) == eta
        expect = (math.log(L.A) - np.euler_gamma) / (2 * twist_r1(E11, eta, 2**14))
        assert twist_effective_n(L, 2**14) == pytest.approx(expect, rel=1e-14)


@pytest.mark.parametrize("w", [1, -1])
def test_twist_rho_expansion(w):
    r1 = {1: 1.283598747733984, -1: 0.8439846143876162}
    for la in (60.0, 400.0):
        L1 = la - np.euler_gamma
        t = np.linspace(1e-4, 3.5 * math.pi / L1, 60)
        prof = twist_profile(E11, t, P=2**14)
        for eta in (1, -1):
            x, rho = twist_rho(t, *prof[eta], la, w)
            first = 1 + w * np.sinc(2 * x)
            pred = first - r1[eta] * (1 + w * np.cos(2 * math.pi * x)) / L1
            assert np.max(np.abs(rho - pred)) < 100 / L1**2
            assert np.max(np.abs(rho - first)) > r1[eta] / L1


def test_twist_density_ends():
    t = np.linspace(1e-5, 0.5, 30)
    prof = twist_profile(E11, t, P=2**14)
    la = math.log(math.sqrt(11) * 2000 / (2 * math.pi))
    x, rho = twist_rho(t, *prof[-1], la, -1)
    assert rho[0] == pytest.approx(0, abs=1e-3)
    x, rho = twist_rho(t, *prof[1], la, 1)
    r1 = twist_r1(E11, 1, 2**14)
    assert rho[0] == pytest.approx(2 - 2 * r1 / (la - np.euler_gamma), abs=1e-3)


def test_twist_density_is_the_family_mean():
    x = np.linspace(0, 2, 41)
    one = [twist_density(E11, [d], x, P=2**12, n=150) for d in (89, 1201)]
    both = twist_density(E11, [89, 1201], x, P=2**12, n=150)
    assert np.allclose(both, (one[0] + one[1]) / 2, atol=2e-3)
