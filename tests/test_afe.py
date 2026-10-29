import cmath
import math

import mpmath
import numpy as np
import pytest

from zeta_explorer import afe
from zeta_explorer.afe import afe_length, gammaincc, l_afe
from zeta_explorer.characters import conrey, quadratic


def test_gammaincc_matches_mpmath():
    rng = np.random.default_rng(7)
    for _ in range(80):
        z = complex(rng.uniform(0.05, 3), rng.uniform(-300, 300))
        w = 10 ** rng.uniform(-4, 3) * cmath.exp(1j * math.copysign(rng.uniform(0, 1.56), z.imag))
        ref = complex(mpmath.gammainc(z, w, regularized=True))
        assert abs(gammaincc(z, w) - ref) <= 1e-11 * max(1.0, abs(ref))


def test_gammaincc_order_one_is_exponential():
    w = np.array([1e-3, 0.7 + 2j, 5 - 30j, 80 + 400j])
    assert np.allclose(gammaincc(1.0, w), np.exp(-w), rtol=1e-13, atol=0)


def test_gammaincc_recurrence():
    w = np.array([0.02 + 0.1j, 1.5 + 3j, 12 + 40j, 60 + 0.5j])
    for z in (0.3 + 4j, 1.7 - 20j):
        step = np.exp(z * np.log(w) - w) / complex(mpmath.gamma(z + 1))
        assert np.allclose(gammaincc(z + 1, w) - gammaincc(z, w), step, rtol=1e-12, atol=1e-14)


def test_gammaincc_rejects_negative_axis():
    with pytest.raises(ValueError):
        gammaincc(0.5, -2.0)
    with pytest.raises(ValueError):
        gammaincc(-0.5, 2.0)


@pytest.mark.parametrize("q, n", [(5, 2), (7, 3), (8, 3), (12, 11), (13, 2), (21, 5), (32, 3),
                                  (40, 3), (97, 5)])
def test_afe_matches_mpmath(q, n):
    chi = conrey(q, n)
    rng = np.random.default_rng(q + n)
    s = rng.uniform(0.02, 0.98, 8) + 1j * rng.uniform(-200, 200, 8)
    for si, got in zip(s, l_afe(s, chi), strict=True):
        ref = complex(mpmath.dirichlet(si, list(chi.values)))
        assert abs(got - ref) <= 1e-12 * abs(ref)


@pytest.mark.parametrize("t", [0.0, 2.0, 40.0, 400.0])
def test_afe_independent_of_rotation(monkeypatch, t):
    chi = conrey(1009, 11)
    s = 0.37 + 1j * t
    values = []
    for loss in (1.5, 3.0, 6.0):
        monkeypatch.setattr(afe, "_LOSS", loss)
        values.append(l_afe(s, chi))
    assert max(abs(v - values[0]) for v in values) <= 2e-12 * abs(values[0])


def test_real_character_is_real_on_the_axis():
    assert l_afe(0.5, quadratic(-1003)).imag == pytest.approx(0, abs=1e-14)
    assert l_afe(0.5, quadratic(1005)).imag == pytest.approx(0, abs=1e-14)


def test_length_grows_like_sqrt_qt():
    assert afe_length(100 * 1009, 1.0) / afe_length(1009, 1.0) == pytest.approx(10, rel=0.01)
    assert afe_length(1009, 4000.0) / afe_length(1009, 40.0) == pytest.approx(10, rel=0.01)
    assert afe_length(1009, 0.0) == afe_length(1009, 3.5)


def test_afe_rejects_bad_input():
    with pytest.raises(ValueError):
        l_afe(0.5, conrey(7, 1))
    with pytest.raises(ValueError):
        l_afe(0.5, conrey(21, 4))
    with pytest.raises(ValueError):
        l_afe(1.2 + 3j, conrey(7, 3))
