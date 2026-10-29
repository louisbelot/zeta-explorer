import cmath
import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.afe import gammaincc


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
