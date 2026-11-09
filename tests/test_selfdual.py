import math

import mpmath
import numpy as np
import pytest

from zeta_explorer import delta
from zeta_explorer.selfdual import (
    SelfDual,
    c0,
    multiplicative,
    s_integral_bound,
    turing,
)
from zeta_explorer.sympower import delta_power


def toy(mu=(0.0,), q=1.0):
    return SelfDual(lambda M: np.ones(M + 1), q, mu)


def test_multiplicative_builds_the_divisor_function():
    d = multiplicative(60, lambda p, K: list(range(1, K + 2)))
    assert d[1:13].tolist() == [1, 2, 2, 3, 2, 4, 2, 4, 3, 4, 2, 6]
    assert d[60] == 12
    assert d[49] == 3


def test_theta_matches_the_gamma_c_form():
    L = delta_power(1)
    t = np.array([0.0, 0.7, 9.0, 55.0])
    assert L.theta(t) == pytest.approx(delta.theta(t), abs=1e-12)
    assert L.density(55.0) == pytest.approx(
        (delta.theta(55.01) - delta.theta(54.99)) / 0.02 / math.pi, abs=1e-6
    )


def test_degree_two_agrees_with_the_incomplete_gamma_code():
    L = delta_power(1)
    for s in (0.5 + 10j, 1.2 + 3j, 0.5 + 0j, 0.2 - 6j):
        assert L(s) == pytest.approx(delta.l_delta(s), abs=1e-12)
    assert L.zeros(30) == pytest.approx(delta.zeros(30), abs=1e-10)


def test_hardy_is_the_whole_value_on_the_line():
    L = delta_power(2)
    t = np.array([3.0, 11.0, 26.5])
    assert np.abs(L.hardy(t)) == pytest.approx(np.abs(L(0.5 + 1j * t)), abs=1e-11)


def test_c0_is_booker_minus_one_log_zeta():
    assert c0(0.5) == pytest.approx(4.690288, abs=1e-6)
    assert c0(0.5) + float(mpmath.log(mpmath.zeta(1.5))) == pytest.approx(5.65055, abs=1e-5)
    assert c0(0.1) > c0(0.3) > c0(0.5)


def test_bound_at_half_has_bookers_coefficients():
    L = toy(mu=(0.0, 1.0), q=7.0)
    t1, t2 = 40.0, 55.0

    def log_q(s):
        return math.log(7) + sum(math.log(abs(s + m) / (2 * math.pi)) for m in L.mu)

    X = math.sqrt(t1 * t1 - 3.5**2)
    a = 2 * math.log(2)
    want = (log_q(3.5 + 1j * t2) / 4 + (math.log(2) - 0.5) * log_q(1.5 + 1j * t1)
            + 2 * c0(0.5) + 2 * (8.3 + 0.09 * a + 0.5) / X)
    assert s_integral_bound(L, t1, t2) == pytest.approx(want / math.pi, rel=1e-13)


def test_bound_needs_t1_past_the_gamma_shifts():
    L = toy(mu=(1.0, 11.0, 12.0))
    with pytest.raises(ValueError):
        s_integral_bound(L, 15.0, 30.0)
    with pytest.raises(ValueError):
        s_integral_bound(L, 30.0, 29.0)
    with pytest.raises(ValueError):
        s_integral_bound(L, 30.0, 40.0, eps=0.6)
    assert s_integral_bound(L, 16.0, 30.0) > s_integral_bound(L, 30.0, 40.0)


def test_turing_area_counts_zeros_inside_the_window():
    L = toy()
    full = turing(L, [1.0, 2.0, 7.5], 6.0, 10.0)
    less = turing(L, [2.0, 7.5], 6.0, 10.0)
    late = turing(L, [1.0, 2.0], 6.0, 10.0)
    assert full.found == 2
    assert full.area - less.area == pytest.approx(4.0, abs=1e-12)
    assert full.area - late.area == pytest.approx(2.5, abs=1e-12)
    assert less.missing - full.missing == pytest.approx(1.0, abs=1e-12)


def test_turing_certifies_sym2_delta_and_notices_a_dropped_zero():
    L = delta_power(2)
    z = L.zeros(40)
    r = turing(L, z, 30, 40)
    assert r.found == 13
    assert abs(r.found - r.expected) < 0.5
    assert r.complete
    assert not turing(L, np.delete(z, 4), 30, 40).complete
