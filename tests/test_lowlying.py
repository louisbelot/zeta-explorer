import math

import numpy as np
import pytest

from zeta_explorer.characters import conrey, quadratic
from zeta_explorer.lfunction import l_theta, l_zeros
from zeta_explorer.lowlying import (
    OneLevel,
    low_zeros,
    one_level,
    prime_family,
    quadratic_family,
    unfold_l,
    w_symplectic,
    w_unitary,
)


def test_unfolding_counts_zeros():
    chi = conrey(1009, 11)
    T = 60.0
    zeros = l_zeros(chi, 0, T)
    assert abs(len(zeros) - unfold_l(T, chi)) < 3


def test_unfolding_is_odd():
    chi = conrey(211, 2)
    t = np.array([0.3, 2.0, 9.0])
    assert np.allclose(unfold_l(-t, chi), -unfold_l(t, chi))


def test_unfolded_density_is_one_on_average():
    chi = conrey(1013, 3)
    x = unfold_l(l_zeros(chi, 0, 80), chi)
    gaps = np.diff(x)
    assert np.mean(gaps) == pytest.approx(1, abs=0.03)


def test_low_zeros_are_sorted_and_bounded():
    x = low_zeros(quadratic(-1003), 4.0)
    assert np.all(np.diff(x) > 0)
    assert 0 < x[0] and x[-1] <= 4.0


@pytest.mark.parametrize("q, n", [(5, 2), (7, 3), (101, 4), (211, 4)])
def test_low_zeros_refuse_small_conductors(q, n):
    with pytest.raises(ValueError):
        low_zeros(conrey(q, n), 3.0)


@pytest.mark.parametrize("q, n", [(11, 2), (227, 4), (1009, 4)])
def test_theta_increases_from_zero_above_the_threshold(q, n):
    chi = conrey(q, n)
    t = np.linspace(0, 5, 200)
    assert np.all(np.diff(l_theta(t, chi)) > 0)
    assert low_zeros(chi, 3.0)[0] > 0


def test_quadratic_family_counts_fundamental_discriminants():
    fam = quadratic_family(1, 30)
    assert len(fam) == 10 + 9
    assert all(chi.is_real and chi.is_primitive for chi in fam)


def test_prime_family_is_all_nontrivial_characters():
    fam = prime_family(101)
    assert len(fam) == 99
    assert sum(chi.is_real for chi in fam) == 1


def test_density_normalization():
    fam = OneLevel([np.array([0.5, 1.5]), np.array([0.25]), np.array([])])
    centers, d = fam.density(bins=4, x_max=2.0)
    assert np.allclose(centers, [0.25, 0.75, 1.25, 1.75])
    assert np.allclose(d, np.array([1, 1, 0, 1]) / (3 * 0.5))
    assert np.allclose(fam.first(), [0.5, 0.25])


def test_symplectic_density():
    assert w_symplectic(0.0) == 0
    assert w_symplectic(np.array([0.5, 1.0])) == pytest.approx([1, 1])
    assert w_symplectic(0.25) == pytest.approx(1 - 2 / math.pi)
    assert np.all(w_unitary([0.0, 2.0]) == 1)


def test_quadratic_family_avoids_the_central_point():
    quad = one_level(quadratic_family(1000, 1150), 1.0)
    unit = one_level(prime_family(1009)[:150], 1.0)
    _, dq = quad.density(bins=4, x_max=1.0)
    _, du = unit.density(bins=4, x_max=1.0)
    assert dq[0] < 0.45 < du[0]
