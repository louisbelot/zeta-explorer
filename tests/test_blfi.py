import math

import numpy as np
import pytest

from zeta_explorer.blfi import kernel, subdivide


@pytest.mark.parametrize("n2", [1, 2, 7, 64, 65, 127, 1000, 39894, 398942])
@pytest.mark.parametrize("M", [1, 3, 16, 100])
def test_blocks_tile_the_main_sum(n2, M):
    if M > n2:
        with pytest.raises(ValueError):
            subdivide(n2, M)
        return
    v, K = subdivide(n2, M)
    assert v[0] == M
    assert np.all(v[1:] == v[:-1] + K[:-1])
    assert v[-1] + K[-1] - 1 == n2
    assert np.all((K - 1) * M < v)
    assert np.all(K >= 1)
    m = (n2 // M).bit_length() - 1
    assert len(v) <= M * (m + 1) + 1
    assert set(K.tolist()) <= {1 << j for j in range(m + 1)} | {int(K[-1])}


def test_kernel_rebuilds_a_band_limited_sum():
    rng = np.random.default_rng(5)
    tau = 0.01
    freq = rng.uniform(-tau, tau, 200)
    amp = rng.standard_normal(200) + 1j * rng.standard_normal(200)

    def G(t):
        return np.exp(1j * np.multiply.outer(t, freq)) @ amp

    sup = np.abs(amp).sum()
    h = math.pi / (3 * tau)
    for c in (10.0, 20.0, 30.0):
        grid = np.arange(-120, 121) * h
        t = rng.uniform(-3000, 3000, 50)
        approx = kernel(t[:, None] - grid, tau, c) @ G(grid)
        assert np.abs(approx - G(t)).max() < 6 * sup * math.exp(-c)


def test_kernel_has_finite_reach():
    tau, c = 0.01, 25.0
    assert kernel(0.0, tau, c) == pytest.approx(2 / 3)
    assert kernel(np.array([c / tau, -c / tau, 1.01 * c / tau]), tau, c) == pytest.approx(0)
