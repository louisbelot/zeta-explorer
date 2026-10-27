import math
from fractions import Fraction

import mpmath
import numpy as np
import pytest

from zeta_explorer.characters import (
    class_number,
    conrey,
    fundamental_unit,
    is_fundamental,
    primitive_characters,
    quadratic,
    unit_count,
)
from zeta_explorer.lfunction import (
    count_zeros_box,
    dirichlet_l,
    l_theta,
    l_z,
    l_zeros,
    verify_zeros,
)
from zeta_explorer.zeta import bernoulli, zeta

LMFDB_ZEROS = {
    (3, 2): [8.03973715568146668, 11.24920620777293525, 15.70461917672162557],
    (4, 3): [6.02094890469759665, 10.24377030416655455, 12.98809801231242251],
    (5, 2): [6.18357819545085391, 8.45722917442323072, 12.67494641701135578],
}


def bernoulli_poly(n, x):
    return sum(math.comb(n, k) * bernoulli(k) * x ** (n - k) for k in range(n + 1))


@pytest.mark.parametrize("q, n", [(5, 2), (7, 3), (8, 3), (9, 2), (12, 11), (13, 2), (15, 2),
                                  (21, 5), (32, 3), (40, 3)])
def test_l_matches_mpmath(q, n):
    chi = conrey(q, n)
    rng = np.random.default_rng(q * n)
    s = rng.uniform(-15, 15, 12) + 1j * rng.uniform(-80, 80, 12)
    got = dirichlet_l(s, chi)
    for si, g in zip(s, got, strict=True):
        ref = complex(mpmath.dirichlet(si, list(chi.values)))
        assert abs(g - ref) <= 1e-12 * abs(ref)


@pytest.mark.parametrize("t", [500.0, 2000.0])
def test_l_high_on_critical_line(t):
    chi = conrey(7, 3)
    ref = complex(mpmath.dirichlet(0.5 + 1j * t, list(chi.values)))
    assert abs(dirichlet_l(0.5 + 1j * t, chi) - ref) < 1e-10


def test_l_at_one():
    assert dirichlet_l(1, quadratic(-4)) == pytest.approx(math.pi / 4, rel=1e-15)
    assert dirichlet_l(1, quadratic(-3)) == pytest.approx(math.pi / (3 * math.sqrt(3)), rel=1e-15)
    assert dirichlet_l(1, quadratic(5)) == pytest.approx(
        2 * math.log((1 + math.sqrt(5)) / 2) / math.sqrt(5), rel=1e-15
    )


def test_l_near_one_has_no_cancellation():
    chi = quadratic(-4)
    h = 1e-9
    slope = (dirichlet_l(1 + h, chi) - dirichlet_l(1 - h, chi)) / (2 * h)
    ref = complex(mpmath.diff(lambda s: mpmath.dirichlet(s, [0, 1, 0, -1]), 1))
    assert abs(slope - ref) < 1e-6


@pytest.mark.parametrize("q, n", [(4, 3), (5, 2), (7, 3), (8, 5), (9, 4), (11, 2), (24, 5)])
@pytest.mark.parametrize("k", range(1, 8))
def test_generalized_bernoulli_numbers(q, n, k):
    chi = conrey(q, n)
    b = sum(chi(a) * float(bernoulli_poly(k, Fraction(a, q))) for a in range(1, q + 1))
    expected = -(q ** (k - 1)) * b / k
    assert abs(dirichlet_l(1 - k, chi) - expected) < 1e-12 * max(1, abs(expected))


@pytest.mark.parametrize("q, n", [(4, 3), (5, 4), (7, 3), (13, 5)])
def test_trivial_zeros(q, n):
    chi = conrey(q, n)
    s = -chi.parity - 2 * np.arange(6)
    assert np.all(dirichlet_l(s, chi) == 0)


def test_l_is_real_for_real_characters_on_the_real_axis():
    values = dirichlet_l(np.array([-3.5, -1.0, 0.25, 0.5, 2.0]), quadratic(-23))
    assert np.all(values.imag == 0)


def test_principal_character_is_zeta_with_euler_factors():
    s = np.array([2.0, 0.5 + 14j, -1.5 + 3j])
    expected = zeta(s) * (1 - 2.0**-s) * (1 - 3.0**-s)
    assert np.allclose(dirichlet_l(s, conrey(12, 1)), expected, rtol=1e-13)


def test_imprimitive_character_drops_euler_factors():
    chi = conrey(20, 13)
    assert chi.conductor == 5
    s = np.array([2.0, 0.5 + 10j, -2.5 + 1j])
    prim = chi.primitive()
    expected = dirichlet_l(s, prim) * (1 - prim(2) * 2.0**-s)
    assert np.allclose(dirichlet_l(s, chi), expected, rtol=1e-13)


@pytest.mark.parametrize("q", [5, 7, 8, 12, 13, 16, 21])
def test_functional_equation(q):
    rng = np.random.default_rng(q)
    s = rng.uniform(-3, 4, 6) + 1j * rng.uniform(-30, 30, 6)
    for _, chi in primitive_characters(q):
        k = chi.parity

        def completed(z, c, k=k):
            gamma = np.array([complex(mpmath.loggamma((x + k) / 2)) for x in z])
            return np.exp((z + k) / 2 * math.log(q / math.pi) + gamma) * dirichlet_l(z, c)

        rhs = chi.root_number() * completed(1 - s, chi.conj())
        assert np.allclose(completed(s, chi), rhs, rtol=1e-11)


@pytest.mark.parametrize("d", [d for d in range(-3, -400, -1) if is_fundamental(d)])
def test_class_number_formula_imaginary(d):
    L = dirichlet_l(1, quadratic(d))
    h = class_number(d)
    assert L == pytest.approx(2 * math.pi * h / (unit_count(d) * math.sqrt(-d)), rel=1e-13)


@pytest.mark.parametrize("d, h", [(5, 1), (40, 2), (60, 2), (229, 3), (316, 3), (328, 4)])
def test_class_number_formula_real(d, h):
    x, y = fundamental_unit(d)
    eps = (x + y * math.sqrt(d)) / 2
    L = dirichlet_l(1, quadratic(d))
    assert L * math.sqrt(d) / (2 * math.log(eps)) == pytest.approx(h, rel=1e-13)


@pytest.mark.parametrize("q, n", [(5, 2), (7, 3), (8, 3), (11, 2), (13, 2), (40, 3)])
def test_hardy_z_is_real_and_has_modulus_of_l(q, n):
    chi = conrey(q, n)
    t = np.linspace(-40, 40, 101)
    rotated = np.exp(1j * l_theta(t, chi)) * dirichlet_l(0.5 + 1j * t, chi)
    assert np.max(np.abs(rotated.imag)) < 1e-12
    assert np.allclose(np.abs(l_z(t, chi)), np.abs(dirichlet_l(0.5 + 1j * t, chi)), atol=1e-13)


@pytest.mark.parametrize("label", LMFDB_ZEROS)
def test_zeros_match_lmfdb(label):
    zeros = l_zeros(conrey(*label), 0, 16)[:3]
    assert np.allclose(zeros, LMFDB_ZEROS[label], rtol=0, atol=1e-13)


def test_conjugate_character_has_mirrored_zeros():
    chi = conrey(7, 3)
    assert np.allclose(l_zeros(chi.conj(), -30, 30), -l_zeros(chi, -30, 30)[::-1], atol=1e-12)


def test_real_character_zeros_are_symmetric():
    zeros = l_zeros(quadratic(-23), -30, 30)
    assert np.allclose(zeros, -zeros[::-1], atol=1e-12)


@pytest.mark.parametrize("q, n, T", [(4, 3, 40.0), (7, 3, 40.0), (12, 11, 30.0), (101, 2, 25.0)])
def test_verify_zeros(q, n, T):
    v = verify_zeros(conrey(q, n), T)
    assert v.ok
    assert v.in_box == len(v.zeros) > 0


def test_box_count_sees_zeros_off_the_line_too():
    chi = conrey(7, 3)
    assert count_zeros_box(chi, 20.0, sigma=(0.6, 2.0)) == 0
    assert count_zeros_box(chi, 20.0) == len(l_zeros(chi, -20, 20))


@pytest.mark.parametrize("q, n", [(4, 3), (5, 4)])
def test_box_count_skips_trivial_zeros(q, n):
    chi = conrey(q, n)
    assert count_zeros_box(chi, 15.0, sigma=(-4.5, 2.0)) == count_zeros_box(chi, 15.0)


def test_box_count_rejects_zeta():
    with pytest.raises(ValueError):
        count_zeros_box(conrey(9, 1), 10.0)
