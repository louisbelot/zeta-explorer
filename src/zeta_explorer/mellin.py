from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.selfdual import SelfDual, _log_gamma_factor, _psi_gamma_factor, scan_zeros

_C = 1.0
_TAIL = 46.0


def _log_weight(L: SelfDual, t: np.ndarray, theta: float, re: float = 0.5) -> np.ndarray:
    return theta * t + _log_gamma_factor(re + 1j * t, L.mu).real


def rotation(L: SelfDual, t: ArrayLike) -> np.ndarray:
    return _psi_gamma_factor(0.5 + 1j * np.asarray(t, dtype=float), L.mu).imag


def windows(L: SelfDual, T: float, drop: float = 8.0) -> list[tuple[float, float, float, float]]:
    grid = np.linspace(0.0, 2 * T + 40, 1 + int(40 * (T + 20)))
    ell = _log_gamma_factor(0.5 + 1j * grid, L.mu).real
    h = rotation(L, grid)
    out = []
    a = 0
    while grid[a] < T:
        rise = h[a:] * (grid[a:] - grid[a]) + ell[a:] - ell[a]
        peak = a + max(1, int(np.searchsorted(rise, drop, side="right")) - 1)
        fall = h[peak] * (grid[peak] - grid[peak:]) + ell[peak] - ell[peak:]
        b = peak + int(np.searchsorted(fall, drop, side="right")) - 1
        if grid[b] >= T:
            b = int(np.searchsorted(grid, T))
            c = np.arange(a + 1, min(peak, b) + 1)
            ok = h[c] * (grid[c] - grid[b]) + ell[c] - ell[b] <= drop
            peak = int(c[np.argmax(ok)])
        out.append((float(grid[a]), float(min(grid[b], T)), float(grid[peak]), float(h[peak])))
        a = b
    return out


def _extent(f: Callable[[np.ndarray], np.ndarray], floor: float, R: float) -> tuple[float, float]:
    while f(np.array([-R, R])).max() > floor:
        R *= 2
    x = np.linspace(-R, R, 8193)
    keep = np.flatnonzero(f(x) > floor)
    step = x[1] - x[0]
    return float(x[keep[0]] - step), float(x[keep[-1]] + step)


def terms(L: SelfDual, theta: float) -> int:
    d = len(L.mu)
    y = np.geomspace(1e-3, 1e15, 7201)
    f = (sum(L.mu) - d + 1) / d * np.log(y) - d * math.pi * math.cos(2 * theta / d) * y ** (2 / d)
    i = int(np.argmax(f))
    j = i + int(np.argmax(f[i:] < f[i] - _TAIL - np.log(y[i:] / y[i])))
    return max(2, math.ceil(math.sqrt(L.q) * y[j]))


@dataclass(frozen=True)
class Band:
    lo: float
    hi: float
    theta: float
    ref: float
    u: np.ndarray
    g: np.ndarray
    terms: int

    def integral(self, t: np.ndarray) -> np.ndarray:
        out = np.empty(len(t), dtype=complex)
        for i in range(0, len(t), 512):
            out[i : i + 512] = np.exp(1j * np.outer(t[i : i + 512], self.u)) @ self.g
        return out


def band(L: SelfDual, lo: float, hi: float, peak: float, theta: float) -> Band:
    M = terms(L, theta)
    U = math.log(M)
    h_tau = 2 * math.pi / (U + _TAIL / _C)
    ref = float(_log_weight(L, np.array([peak]), theta)[0])
    R = 2 * peak + 32

    edge = _log_weight(L, np.array([lo, hi]), theta).min()
    t_lo, t_hi = _extent(lambda t: _log_weight(L, t, theta), edge - _TAIL, R)

    def log_f(tau: np.ndarray) -> np.ndarray:
        return _log_weight(L, tau, theta, _C)

    top = log_f(np.linspace(-R, R, 2049)).max()
    tau_lo, tau_hi = _extent(log_f, top - _TAIL, R)

    k = np.arange(math.floor(tau_lo / h_tau), math.ceil(tau_hi / h_tau) + 1)
    need = max(2 * np.abs(k).max() + 2, (t_hi - t_lo) / h_tau)
    N = 1 << math.ceil(math.log2(need))
    tau = k * h_tau
    s = _C + 1j * tau

    n = np.arange(1, M + 1)
    a = L.coeffs(1 << max(10, (M - 1).bit_length()))[1 : M + 1] * n ** -_C
    ln = np.log(n)
    D = np.zeros(len(tau), dtype=complex)
    for i in range(0, M, 1024):
        D += np.exp(-1j * np.outer(tau, ln[i : i + 1024])) @ a[i : i + 1024]

    log_a = (_log_gamma_factor(s, L.mu) + s * math.log(L.q) / 2
             + theta * tau - 1j * theta * (_C - 0.5) - ref)
    padded = np.zeros(N, dtype=complex)
    padded[k % N] = np.exp(log_a) * D

    h_u = 2 * math.pi / (N * h_tau)
    J = math.ceil(U / h_u) + 1
    u = np.arange(J) * h_u
    g = h_u * h_tau / (2 * math.pi) * np.exp(-u * (_C - 0.5)) * np.fft.fft(padded)[:J]
    g[0] /= 2
    return Band(lo, hi, theta, ref, u, g, M)


class MellinLine:
    def __init__(self, L: SelfDual, T: float, drop: float = 8.0):
        self.L = L
        self.T = float(T)
        self.bands = [band(L, *w) for w in windows(L, T, drop)]
        self._edges = np.array([b.hi for b in self.bands[:-1]])

    def hardy(self, t: ArrayLike) -> float | np.ndarray:
        ts = np.atleast_1d(np.asarray(t, dtype=float)).ravel()
        which = np.searchsorted(self._edges, ts, side="left")
        out = np.empty(len(ts))
        for i, b in enumerate(self.bands):
            m = which == i
            if not m.any():
                continue
            w = b.integral(ts[m])
            part = w.real if self.L.eps == 1 else w.imag
            scale = b.ref - _log_weight(self.L, ts[m], b.theta) - math.log(self.L.q) / 4
            out[m] = 2 * part * np.exp(scale)
        return float(out[0]) if np.ndim(t) == 0 else out.reshape(np.shape(t))

    def zeros(self, per_gap: int = 8, rtol: float = 1e-12) -> np.ndarray:
        return scan_zeros(self.hardy, self.L.density, self.T, per_gap, rtol)
