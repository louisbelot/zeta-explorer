import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.blfi import BLFIWindow, kernel, log_split, spin, subdivide
from zeta_explorer.multieval import ZWindow


def _siegelz(t0, u):
    with mpmath.workdps(35):
        return float(mpmath.siegelz(mpmath.mpf(t0) + mpmath.mpf(u)))


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


def test_spin_keeps_full_precision_far_from_t0():
    rng = np.random.default_rng(1)
    n = rng.integers(2, 10**7, 60)
    u = rng.uniform(-60000, 60000, 60)
    phi = rng.uniform(0, 2 * math.pi, 60)
    got = spin(phi, u, *log_split(n))
    with mpmath.workdps(40):
        exact = [
            complex(mpmath.expj(mpmath.mpf(f) + mpmath.mpf(x) * mpmath.log(int(k))))
            for f, x, k in zip(phi, u, n, strict=True)
        ]
    assert np.abs(got - exact).max() < 1e-14


def test_log_split_puts_the_head_on_a_grid():
    hi, lo = log_split([2, 3, 10**6 + 3])
    assert np.all(hi * 2**16 == np.rint(hi * 2**16))
    assert np.all(np.abs(lo) <= 2**-17)
    with mpmath.workdps(40):
        for k, h, low in zip((2, 3, 10**6 + 3), hi, lo, strict=True):
            assert abs(mpmath.mpf(h) + mpmath.mpf(low) - mpmath.log(k)) < 1e-21


@pytest.mark.parametrize("T", [1e4, 1e8, 1e10])
def test_window_matches_mpmath(T):
    win = BLFIWindow(T, T + T**0.25)
    for u in np.random.default_rng(int(math.log10(T))).uniform(win.u_lo, win.u_hi, 3):
        assert win.z(u)[0] == pytest.approx(_siegelz(win.t0, u), abs=2e-12)


def test_main_sum_agrees_with_the_nufft_window():
    T = 3e9
    win = BLFIWindow(T, T + 300)
    other = ZWindow(T, T + 300)
    u = np.linspace(win.u_lo, win.u_hi, 777)
    assert np.abs(win.main_sum(u) - np.conj(other.main_sum(u))).max() < 1e-11


def test_tolerance_bounds_the_error():
    T = 1e8
    ref = ZWindow(T, T + 100)
    u = np.linspace(ref.u_lo, ref.u_hi, 200)
    exact = ref.z(u)
    errors = []
    for tol in (1e-3, 1e-6, 1e-9):
        errors.append(np.abs(BLFIWindow(T, T + 100, tol=tol).z(u) - exact).max())
        assert errors[-1] < tol
    assert errors[0] > errors[1] > errors[2]


def test_zeta_matches_mpmath():
    T = 1e6
    win = BLFIWindow(T, T + 40)
    for u in (-17.1, 0.3, 19.9):
        with mpmath.workdps(30):
            exact = complex(mpmath.zeta(mpmath.mpc(0.5, mpmath.mpf(win.t0) + mpmath.mpf(u))))
        assert win.zeta(u)[0] == pytest.approx(exact, abs=1e-12)


def test_window_across_a_change_in_the_number_of_terms():
    k = 200
    edge = 2 * math.pi * k * k
    win = BLFIWindow(edge - 300, edge + 300)
    assert (win.n2, win.n) == (k - 1, k)
    for u in edge - win.t0 + np.array([-200.0, -0.01, 0.01, 150.0]):
        assert win.z(u)[0] == pytest.approx(_siegelz(win.t0, u), abs=1e-12)


def test_order_of_points_does_not_matter():
    win = BLFIWindow(1e7, 1e7 + 60)
    u = np.random.default_rng(3).uniform(win.u_lo, win.u_hi, 300)
    assert win.z(u) == pytest.approx(win.z(np.sort(u))[np.argsort(np.argsort(u))], abs=1e-15)


def test_block_count_grows_like_the_fourth_root():
    counts = [len(BLFIWindow(T, T + 1).v) for T in (1e6, 1e8, 1e10)]
    for small, big in zip(counts, counts[1:], strict=False):
        assert 10**0.5 < big / small < 10**0.5 * 2
    assert counts[-1] < math.sqrt(1e10 / (2 * math.pi)) / 15


def test_window_rejects_bad_input():
    win = BLFIWindow(1000, 1010)
    with pytest.raises(ValueError):
        win.z(win.u_hi + 1)
    with pytest.raises(ValueError):
        BLFIWindow(10, 20)
    with pytest.raises(ValueError):
        BLFIWindow(1000, 1010, M=0)
    with pytest.raises(ValueError):
        BLFIWindow(1000, 1010, tol=0)
    with pytest.raises(ValueError):
        BLFIWindow(1e12, 1e12 + 2**17)
