import cmath
import math

import mpmath
import numpy as np
import pytest
from flint import arb

from zeta_explorer import certified as C
from zeta_explorer import davenport as D


def test_coefficients_come_from_the_two_characters():
    c = (1 - 1j * D.KAPPA) / 2
    a = [2 * (c * D.CHI(n)).real for n in range(5)]
    assert np.allclose(a, D.coefficients(), atol=1e-15)


def test_kappa_cancels_the_root_number():
    c = (1 - 1j * D.KAPPA) / 2
    assert abs(D.CHI.root_number() - c.conjugate() / c) < 1e-14
    assert abs(math.atan(D.KAPPA) - cmath.phase(D.CHI.root_number()) / 2) < 1e-14


@pytest.mark.parametrize("s", [0.5 + 20j, 0.3 + 5j, 0.81 + 85.7j, 1.3 + 400j, 3 + 1j,
                               -0.6 + 33j, -2.5 + 0.5j, 0.5 + 1500j, 2 - 10j])
def test_dh_matches_arb_hurwitz(s):
    ref = C.dh_ball(s, 128)
    ref = complex(float(ref.real.mid()), float(ref.imag.mid()))
    assert abs(D.dh(s) - ref) < 1e-11 * max(1.0, abs(ref))


def test_dh_is_the_dirichlet_series_for_large_sigma():
    s = 4 + 3j
    a = D.coefficients()
    n = np.arange(1, 200_001)
    assert abs(D.dh(s) - np.sum(a[n % 5] * n ** (-s))) < 1e-14


@pytest.mark.parametrize("s", [0.3 + 5j, 2 + 7j, -0.4 + 33j, 0.7 + 100j, 1.5 - 12j])
def test_completed_function_is_symmetric(s):
    def completed(z):
        w = (z + 1) / 2
        return mpmath.power(D.Q / mpmath.pi, w) * mpmath.gamma(w) * complex(D.dh(z))

    lhs, rhs = completed(s), completed(1 - s)
    assert abs(lhs - rhs) < 1e-11 * abs(lhs)


def test_dh_is_real_on_the_real_axis_with_trivial_zeros():
    x = np.linspace(-6, 4, 41)
    assert np.max(np.abs(np.imag(D.dh(x)))) < 1e-14
    assert np.max(np.abs(D.dh(np.array([-1.0, -3.0, -5.0])))) < 1e-12


def test_hardy_z_is_real():
    t = np.linspace(1, 200, 401)
    v = np.exp(1j * D.dh_theta(t)) * D.dh(0.5 + 1j * t)
    assert np.max(np.abs(v.imag) / np.maximum(1, np.abs(v))) < 1e-12
    assert np.allclose(D.dh_z(t), v.real)


def test_zero_free_sigma_balances_the_coefficients():
    s0 = D.zero_free_sigma()
    with C.precision(128):
        k = ((arb(10) - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
        h = [arb(s0).zeta(arb(a) / 5) for a in (1, 2, 3, 4)]
        total = arb(5) ** (-arb(s0)) * (h[0] + h[3] + k * (h[1] + h[2]))
    assert abs(float(total.mid()) - 2) < 1e-12
    assert 1.39 < s0 < 1.40


def test_no_zeros_right_of_sigma0():
    s0 = D.zero_free_sigma()
    t = np.linspace(0, 600, 6001)
    assert np.min(np.abs(D.dh(s0 + 1e-3 + 1j * t))) > 1e-4


def test_line_zeros_match_mpmath():
    zs = D.line_zeros(1, 15)
    a = [mpmath.mpf(x) for x in D.coefficients()]
    with mpmath.workdps(25):
        for z in zs:
            ref = mpmath.findroot(lambda t: mpmath.re(
                mpmath.exp(1j * _theta_mp(t)) * mpmath.dirichlet(0.5 + 1j * t, a)), z)
            assert abs(z - float(ref)) < 1e-11


def _theta_mp(t):
    return t / 2 * mpmath.log(5 / mpmath.pi) + mpmath.im(mpmath.loggamma((1.5 + 1j * t) / 2))


def test_spira_zeros_and_their_mirrors():
    roots, ok = D.newton(list(D.SPIRA))
    assert ok.all()
    assert np.max(np.abs(roots - np.array(D.SPIRA))) < 1e-6
    assert np.max(np.abs(D.dh(roots))) < 1e-12
    assert np.max(np.abs(D.dh(1 - np.conj(roots)))) < 1e-12


def test_spira_zeros_certified_off_the_line():
    roots, _ = D.newton(list(D.SPIRA))
    for r in roots:
        box = C.krawczyk(complex(r), 1e-10)
        assert box is not None
        assert box.real > 0.57
        assert box.rad() < 1e-15
        assert abs(complex(float(box.real.mid()), float(box.imag.mid())) - r) < 1e-12


def test_krawczyk_rejects_a_box_without_a_zero():
    assert C.krawczyk(0.8 + 90j, 1e-3) is None


def test_count_below_the_first_off_line_zero_is_all_on_the_line():
    assert D.count_zeros(1, 80) == len(D.line_zeros(1, 80, 32))


def test_window_with_spira_pair():
    w = D.zeros_in_window(80, 90)
    assert w.total == len(w.line) + 2
    assert np.allclose(sorted(w.off.real), [1 - 0.808517182, 0.808517182], atol=1e-8)


def test_window_with_two_pairs():
    w = D.zeros_in_window(160, 180)
    assert len(w.off) == 4
    right = np.sort(w.off[w.off.real > 0.5].imag)
    assert np.allclose(right, [166.479305913, 176.702461243], atol=1e-8)
