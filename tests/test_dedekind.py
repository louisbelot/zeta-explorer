import math
from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer import dedekind as K
from zeta_explorer.characters import is_fundamental, unit_count

REAL = [d for d in range(5, 1000) if is_fundamental(d)]


@pytest.mark.parametrize("D", [-3, -4, -7, -23, -84, -163])
def test_dirichlet_form_count(D):
    assert np.array_equal(K.form_counts(D, 400), unit_count(D) * K.ideal_counts(D, 400))


def test_gaussian_ideals_are_sums_of_two_squares():
    N = 500
    r2 = np.zeros(N + 1, dtype=np.int64)
    for x in range(-23, 24):
        for y in range(-23, 24):
            if 0 < x * x + y * y <= N:
                r2[x * x + y * y] += 1
    assert np.array_equal(r2, 4 * K.ideal_counts(-4, N))


@pytest.mark.parametrize("D", [-4, 5, 12])
def test_zeta_k_against_dirichlet_series(D):
    N = 20000
    a = K.ideal_counts(D, N)[1:]
    n = np.arange(1, N + 1, dtype=float)
    assert K.dedekind_zeta(3.0, D).real == pytest.approx(np.sum(a / n**3), abs=1e-7)


@pytest.mark.parametrize("D", [-3, -4, -15, -23, -163])
def test_zeta_k_at_zero_imaginary(D):
    assert K.dedekind_zeta(0.0, D).real == pytest.approx(-K.field_class_number(D) / unit_count(D),
                                                        abs=1e-12)


@pytest.mark.parametrize("D", [5, 12, 40, 229])
def test_zeta_k_derivative_at_zero_real(D):
    h = 1e-5
    deriv = (K.dedekind_zeta(h, D) - K.dedekind_zeta(-h, D)).real / (2 * h)
    assert deriv == pytest.approx(-K.real_class_number(D) * K.regulator(D) / 2, rel=1e-8)


def test_real_class_numbers_match_l1():
    for D in REAL:
        assert K.analytic_class_number(D) == pytest.approx(K.real_class_number(D), abs=1e-9)


@pytest.mark.parametrize("D, h, hp", [(5, 1, 1), (12, 1, 2), (136, 2, 4), (229, 3, 3),
                                      (401, 5, 5), (577, 7, 7), (1957, 3, 6)])
def test_known_real_class_numbers(D, h, hp):
    assert K.real_class_number(D) == h
    assert K.narrow_class_number(D) == hp


def test_cycles_are_reduced_and_closed():
    D = 1957
    forms = K.reduced_indefinite(D)
    assert sum(len(c) for c in K.cycles(D)) == len(forms)
    for a, b, c in forms:
        assert b * b - 4 * a * c == D
        assert K.rho((a, b, c), D) in forms


def test_unit_norm_and_regulator():
    assert K.unit_norm(5) == -1 and K.unit_norm(12) == 1
    assert K.regulator(5) == pytest.approx(math.log((1 + math.sqrt(5)) / 2))
    assert K.regulator(12) == pytest.approx(math.log(2 + math.sqrt(3)))


def test_zeta_minus_one():
    assert K.zeta_minus_one_siegel(5) == Fraction(1, 30)
    assert K.zeta_minus_one_siegel(8) == Fraction(1, 12)
    for D in REAL:
        assert K.zeta_minus_one_siegel(D) == K.zeta_minus_one_bernoulli(D)
    for D in (5, 13, 229):
        assert K.dedekind_zeta(-1.0, D).real == pytest.approx(float(K.zeta_minus_one_siegel(D)),
                                                             rel=1e-12)
    assert K.zeta_minus_one_bernoulli(-23) == 0


def test_superposition_law():
    s = np.linspace(0, 6, 6001)
    p = K.superposition_spacing(s)
    assert p[0] == pytest.approx(0.5, abs=1e-6)
    assert np.trapezoid(p, s) == pytest.approx(1, abs=1e-4)
    assert np.trapezoid(s * p, s) == pytest.approx(1, abs=1e-3)
    F = K.superposition_cdf(s)
    assert np.allclose(np.gradient(F, s)[1:-1], p[1:-1], atol=1e-4)
    assert np.allclose((1 - F)[1:-1], -np.gradient(K.superposition_gap(s), s)[1:-1], atol=1e-4)


def test_gaussian_field_zeros():
    gz, gl = K.zeros(-4, 200.0)
    assert gz[0] == pytest.approx(14.134725141734693)
    assert gl[0] == pytest.approx(6.020948904697597, abs=1e-9)
    u = np.sort(np.concatenate([K.unfold(-4, gz), K.unfold(-4, gl)]))
    assert abs(len(u) - K.unfold(-4, 200.0) - 1) < 2
    assert np.mean(np.diff(u)) == pytest.approx(1, abs=0.05)


def test_needs_fundamental():
    with pytest.raises(ValueError):
        K.ideal_counts(-12, 10)
