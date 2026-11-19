
import numpy as np
import pytest

from zeta_explorer.onelevel import (
    power_sums,
    rmt_trace,
)


def test_power_sums():
    assert power_sums((1, -5, 6), 3) == [5, 13, 35]


def _haar(dim, size, rng):
    Z = rng.standard_normal((size, dim, dim))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]


@pytest.mark.parametrize("dim", [6, 7])
def test_rmt_trace_against_haar(dim):
    Q = _haar(dim, 40_000, np.random.default_rng(dim))
    lam = np.linalg.eigvals(Q)
    det = np.sign(np.linalg.det(Q))
    for n in range(1, dim + 3):
        tr = (lam**n).sum(axis=1).real
        for s in (1, -1):
            x = tr[det == s]
            assert abs(x.mean() - rmt_trace(dim, s, n)) < 5 * x.std() / np.sqrt(len(x)) + 1e-9


def test_rmt_trace_needs_positive_n():
    with pytest.raises(ValueError):
        rmt_trace(4, 1, 0)
