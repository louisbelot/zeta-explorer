import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.bias import (
    bessel_j0,
    l_log_derivative_one,
    observed_density,
    prime_race,
    race,
    square_roots,
    zero_sum,
)
from zeta_explorer.characters import conrey, quadratic, totient
from zeta_explorer.lfunction import l_zeros


@pytest.fixture(scope="module")
def mod4():
    return race(4, 3, 1, 200)


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


@pytest.mark.parametrize(
    "q, a, ref",
    [
        (3, 2, 0.999063),
        (4, 3, 0.995928),
        (5, 2, 0.952175),
        (8, 3, 0.999569),
        (8, 5, 0.997395),
        (8, 7, 0.998939),
        (12, 5, 0.999206),
        (12, 11, 0.999977),
    ],
)
def test_density_matches_fiorilli_martin(q, a, ref):
    assert abs(race(q, a, 1, 200).density() - ref) < 1.5e-6


def test_density_flips_with_the_race(mod4):
    assert abs(race(4, 1, 3, 200).density() + mod4.density() - 1) < 1e-12


def test_symmetric_races_mod_5():
    ref = race(5, 2, 1, 100).density()
    for a, b in ((3, 1), (2, 4), (3, 4)):
        assert abs(race(5, a, b, 100).density() - ref) < 1e-10
    assert race(5, 2, 3, 100).density() == 0.5


def test_density_agrees_with_monte_carlo():
    R = race(5, 2, 1, 40)
    x = R.sample(400_000, np.random.default_rng(11))
    p = R.density()
    sigma = math.sqrt(p * (1 - p) / len(x))
    assert abs(np.mean(x > 0) - p) < 4 * sigma
    assert abs(np.var(x) / R.variance - 1) < 0.01


def test_pdf_integrates_to_density(mod4):
    u = np.linspace(-8, 12, 4001)
    f = mod4.pdf(u)
    assert abs(np.trapezoid(f, u) - 1) < 1e-9
    assert abs(np.trapezoid(u * f, u) - mod4.bias) < 1e-8
    assert abs(np.trapezoid(f[u >= 0], u[u >= 0]) - mod4.density()) < 1e-5


def test_leech_first_lead_change():
    p, d = prime_race(4, 3, 1, 30_000)
    assert p[np.argmax(d < 0)] == 26861
    assert np.all(d[p < 26861] >= 0)


def test_mod_3_race_never_flips_early():
    p, d = prime_race(3, 2, 1, 10**6)
    assert np.all(d[p >= 5] > 0)
    assert observed_density(p, d) > 0.999


@pytest.mark.parametrize("q, a", [(4, 3), (5, 2), (12, 11)])
def test_explicit_formula_tracks_the_race(q, a):
    R = race(q, a, 1, 200)
    p, d = prime_race(q, a, 1, 10**5)
    x = np.exp(np.linspace(math.log(200), math.log(1e5), 600))
    actual = totient(q) * d[np.searchsorted(p, x, side="right") - 1] * np.log(x) / np.sqrt(x)
    err = R.explicit(x) - actual
    assert abs(np.mean(err)) < 0.03
    assert np.mean(np.abs(err)) < 0.25


def test_bad_residues():
    with pytest.raises(ValueError):
        race(4, 2, 1)
    with pytest.raises(ValueError):
        race(5, 2, 7)
