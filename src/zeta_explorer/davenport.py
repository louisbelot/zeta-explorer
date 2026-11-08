from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import conrey
from zeta_explorer.lfunction import arg_change, dirichlet_l, mean_spacing
from zeta_explorer.zeros import illinois_many
from zeta_explorer.zeta import _out, loggamma

Q = 5
KAPPA = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
CHI = conrey(5, 2)
_C = (1 - 1j * KAPPA) / 2

SPIRA = (
    0.808517 + 85.699348j,
    0.650830 + 114.163343j,
    0.574356 + 166.479306j,
    0.724258 + 176.702461j,
)


def coefficients() -> np.ndarray:
    return np.array([0.0, 1.0, KAPPA, -KAPPA, -1.0])


def dh(s: ArrayLike) -> complex | np.ndarray:
    z = np.array(s, dtype=complex, ndmin=1).ravel()
    out = _C * dirichlet_l(z, CHI) + np.conj(_C * dirichlet_l(np.conj(z), CHI))
    return _out(out.reshape(np.shape(s)), s)


def dh_theta(t: ArrayLike) -> float | np.ndarray:
    ts = np.asarray(t, dtype=float)
    out = ts / 2 * math.log(Q / math.pi) + np.imag(loggamma((1.5 + 1j * ts) / 2))
    return float(out) if np.ndim(out) == 0 else out


def dh_z(t: ArrayLike) -> float | np.ndarray:
    ts = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.real(np.exp(1j * dh_theta(ts)) * dh(0.5 + 1j * ts))
    return float(out[0]) if np.ndim(t) == 0 else out


@cache
def zero_free_sigma(dps: int = 20) -> float:
    k = mpmath.mpf(KAPPA)

    def excess(sigma):
        h = [mpmath.zeta(sigma, mpmath.mpf(a) / Q) for a in (1, 2, 3, 4)]
        return Q ** (-sigma) * (h[0] + h[3] + k * (h[1] + h[2])) - 2

    with mpmath.workdps(dps):
        return float(mpmath.findroot(excess, (1.2, 3.0), solver="anderson"))


def line_zeros(t_lo: float, t_hi: float, per_gap: int = 8) -> np.ndarray:
    step = mean_spacing(t_hi, Q) / per_gap
    ts = np.linspace(t_lo, t_hi, max(2, math.ceil((t_hi - t_lo) / step)) + 1)
    z = dh_z(ts)
    i = np.flatnonzero(z[:-1] * z[1:] < 0)
    return illinois_many(dh_z, ts[i], ts[i + 1], z[i], z[i + 1], 1e-14, atol=1e-13)


def count_zeros(t_lo: float, t_hi: float, sigma: tuple[float, float] | None = None) -> int:
    if sigma is None:
        hi = zero_free_sigma() + 0.25
        sigma = (1 - hi, hi)
    lo, hi = sigma
    corners = [lo + 1j * t_lo, hi + 1j * t_lo, hi + 1j * t_hi, lo + 1j * t_hi, lo + 1j * t_lo]
    winding = sum(arg_change(dh, a, b) for a, b in zip(corners[:-1], corners[1:], strict=True))
    n = winding / (2 * math.pi)
    if abs(n - round(n)) > 0.05:
        raise RuntimeError(f"winding number {n:.3f} is not close to an integer")
    return round(n)


def newton(s: ArrayLike, tol: float = 1e-13, max_iter: int = 40) -> tuple[np.ndarray, np.ndarray]:
    z = np.array(s, dtype=complex, ndmin=1)
    ok = np.zeros(len(z), dtype=bool)
    h = 1e-6
    for _ in range(max_iter):
        live = np.flatnonzero(~ok)
        if len(live) == 0:
            break
        w = z[live]
        step = dh(w) / ((dh(w + h) - dh(w - h)) / (2 * h))
        z[live] = w - step
        ok[live] = np.abs(step) < tol * np.maximum(1.0, np.abs(w))
    return z, ok & np.isfinite(z)


def _distinct(z: np.ndarray, tol: float = 1e-7) -> np.ndarray:
    z = np.sort_complex(z)
    keep = np.ones(len(z), dtype=bool)
    for i in range(1, len(z)):
        keep[i] = np.min(np.abs(z[i] - z[:i][keep[:i]])) > tol
    return z[keep]


def off_line_zeros(
    t_lo: float, t_hi: float, sigma_hi: float, d_sigma: float = 0.04, d_t: float = 0.1
) -> np.ndarray:
    sig = np.arange(0.5 - d_sigma, sigma_hi + d_sigma, d_sigma)
    ts = np.arange(t_lo - d_t, t_hi + 2 * d_t, d_t)
    grid = sig[None, :] + 1j * ts[:, None]
    a = np.abs(dh(grid.ravel())).reshape(grid.shape)
    inner = a[1:-1, 1:-1]
    is_min = np.ones_like(inner, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                is_min &= inner <= a[1 + di : a.shape[0] - 1 + di, 1 + dj : a.shape[1] - 1 + dj]
    roots, ok = newton(grid[1:-1, 1:-1][is_min])
    keep = ok & (roots.real > 0.5 + 1e-8) & (roots.imag >= t_lo) & (roots.imag < t_hi)
    right = _distinct(roots[keep])
    return np.sort_complex(np.concatenate([right, 1 - np.conj(right)]))


@dataclass(frozen=True)
class Window:
    t_lo: float
    t_hi: float
    total: int
    line: np.ndarray
    off: np.ndarray


def zeros_in_window(t_lo: float, t_hi: float, levels: int = 4) -> Window:
    sigma_hi = zero_free_sigma() + 0.25
    total = count_zeros(t_lo, t_hi, (1 - sigma_hi, sigma_hi))
    off = off_line_zeros(t_lo, t_hi, sigma_hi)
    for k in range(levels):
        line = line_zeros(t_lo, t_hi, 8 << k)
        if len(line) + len(off) == total:
            return Window(t_lo, t_hi, total, line, off)
    d = 0.04
    for _ in range(levels - 1):
        d /= 2
        off = off_line_zeros(t_lo, t_hi, sigma_hi, d, 2.5 * d)
        if len(line) + len(off) == total:
            return Window(t_lo, t_hi, total, line, off)
    raise RuntimeError(
        f"found {len(line)} + {len(off)} zeros in [{t_lo}, {t_hi}], the box holds {total}"
    )
