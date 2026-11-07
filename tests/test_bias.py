
import mpmath
import numpy as np

from zeta_explorer.bias import (
    bessel_j0,
    square_roots,
)


def test_j0_matches_mpmath():
    z = np.random.default_rng(3).uniform(0, 80, 400)
    ref = np.array([float(mpmath.besselj(0, v)) for v in z])
    assert np.max(np.abs(bessel_j0(z) - ref)) < 1e-14


def test_square_roots():
    assert [square_roots(8, a) for a in (1, 3, 5, 7)] == [4, 0, 0, 0]
    assert [square_roots(7, a) for a in range(1, 7)] == [2, 2, 0, 2, 0, 0]