import math

import mpmath
import numpy as np
import pytest

from zeta_explorer import delta
from zeta_explorer.elliptic import Curve
from zeta_explorer.mellin import MellinLine, rotation, terms, windows
from zeta_explorer.selfdual import turing
from zeta_explorer.sympower import delta_power
from zeta_explorer.symsquare import SymSquare

YOSHIDA_SYM4 = [
    2.3864500, 4.3752457, 6.0435487, 7.571907, 8.841633, 10.605890, 11.437474, 12.76622,
    13.76869, 15.2075, 15.6182, 16.9663, 18.0078, 18.874,
]
DIGITS = [7] * 3 + [6] * 4 + [5] * 2 + [4] * 4 + [3]


def test_rotation_tends_to_a_quarter_turn_per_gamma_factor():
    L = delta_power(3)
    h = rotation(L, [0.0, 5.0, 50.0, 5e5])
    assert h[0] == 0
    assert np.all(np.diff(h) > 0)
    assert h[-1] == pytest.approx(math.pi, abs=1e-4)


def test_windows_tile_the_line_with_growing_rotation():
    L = delta_power(4)
    w = windows(L, 80.0)
    assert w[0][0] == 0 and w[-1][1] == 80.0
    assert all(a[1] == b[0] for a, b in zip(w, w[1:], strict=False))
    thetas = [x[3] for x in w]
    assert thetas == sorted(thetas)
    assert thetas[-1] < 5 * math.pi / 4
    assert len(windows(L, 80.0, drop=4.0)) > len(w)


def test_window_weight_stays_within_the_drop():
    L = delta_power(4)

    def log_w(t, theta):
        return theta * t + sum(float(mpmath.loggamma((0.5 + m + 1j * t) / 2).real) for m in L.mu)

    for lo, hi, peak, theta in windows(L, 60.0, drop=6.0):
        low = min(log_w(t, theta) for t in np.linspace(lo, hi, 25))
        assert log_w(peak, theta) - low < 6.2


def test_more_rotation_needs_more_terms():
    L = delta_power(2)
    n = [terms(L, th) for th in (0.0, 1.0, 2.0, 2.3)]
    assert n == sorted(n)
    assert n[0] < 40 < n[-1]


@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_agrees_with_the_pointwise_functional_equation(m):
    L = delta_power(m)
    line = MellinLine(L, 45.0)
    t = np.array([0.4, 3.3, 11.0, 19.7, 28.1, 33.3, 44.5])
    assert line.hardy(t) == pytest.approx(L.hardy(t), abs=5e-11)
    assert isinstance(line.hardy(7.0), float)


def test_delta_against_the_incomplete_gamma_code():
    line = MellinLine(delta_power(1), 60.0)
    t = np.linspace(0.5, 59.5, 40)
    assert line.hardy(t) == pytest.approx(delta.hardy_z(t), abs=2e-11)
    assert line.zeros() == pytest.approx(delta.zeros(60.0), abs=1e-9)


def test_conductor_enters_through_the_scaling():
    L = SymSquare(Curve(0, -1, 1, 0, 0))
    line = MellinLine(L, 12.0)
    t = np.array([0.0, 2.2, 6.5, 11.9])
    assert line.hardy(t) == pytest.approx(L.hardy(t), abs=2e-11)
    assert line.zeros()[:3] == pytest.approx([3.899281495, 4.734595458, 6.189477388], abs=2e-9)


def test_sym4_zeros_against_yoshida():
    z = MellinLine(delta_power(4), 19.5).zeros()
    assert len(z) == len(YOSHIDA_SYM4)
    for got, want, d in zip(z, YOSHIDA_SYM4, DIGITS, strict=True):
        assert abs(got - want) < 1.01 * 10.0**-d


def test_zeros_pass_turing_at_a_height_the_pointwise_code_cannot_reach_quickly():
    L = delta_power(3)
    z = MellinLine(L, 100.0).zeros()
    r = turing(L, z, 88.0, 100.0)
    assert len(z) == 122
    assert r.complete
    assert not turing(L, np.delete(z, 60), 88.0, 100.0).complete
