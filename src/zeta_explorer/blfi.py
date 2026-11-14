from __future__ import annotations

import math
import time

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.multieval import ZWindow, log_phases, theta_reference, theta_shift
from zeta_explorer.riemann_siegel import TWO_PI, rs_remainder, z_rs

_CHUNK = 1 << 20

_GRID = 2.0**16
_TWO_PI_HI = round(TWO_PI * 2**31) / 2**31
with mpmath.workdps(40):
    _TWO_PI_LO = float(2 * mpmath.pi - _TWO_PI_HI)


def subdivide(n2: int, M: int) -> tuple[np.ndarray, np.ndarray]:
    if not 1 <= M <= n2:
        raise ValueError("need 1 <= M <= n2")
    m = (n2 // M).bit_length() - 1
    starts = [(M << j) + (np.arange(M, dtype=np.int64) << j) for j in range(m)]
    lengths = [np.full(M, 1 << j, dtype=np.int64) for j in range(m)]
    base = M << m
    full = (n2 - base + 1) >> m
    starts.append(base + (np.arange(full, dtype=np.int64) << m))
    lengths.append(np.full(full, 1 << m, dtype=np.int64))
    rest = n2 - base + 1 - (full << m)
    if rest:
        starts.append(np.array([base + (full << m)], dtype=np.int64))
        lengths.append(np.array([rest], dtype=np.int64))
    return np.concatenate(starts), np.concatenate(lengths)


def log_split(n: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    with mpmath.workdps(40):
        logs = [mpmath.log(int(k)) for k in np.asarray(n, dtype=np.int64)]
        hi = [float(mpmath.nint(x * _GRID)) / _GRID for x in logs]
        lo = [float(x - h) for x, h in zip(logs, hi, strict=True)]
    return np.array(hi), np.array(lo)


def cis(x: np.ndarray) -> np.ndarray:
    out = np.empty(np.shape(x), dtype=complex)
    np.cos(x, out=out.real)
    np.sin(x, out=out.imag)
    return out


def spin(phi: ArrayLike, u: ArrayLike, log_hi: ArrayLike, log_lo: ArrayLike) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    u_hi = np.rint(u * _GRID) / _GRID
    p = u_hi * log_hi
    k = np.rint(p / TWO_PI)
    small = u_hi * log_lo + (u - u_hi) * (log_hi + log_lo) + phi
    return cis(p - k * _TWO_PI_HI - k * _TWO_PI_LO + small)


def kernel(d: ArrayLike, tau: float, c: float) -> np.ndarray:
    d = np.asarray(d, dtype=float)
    beta, lam, eps = 3 * tau, 2 * tau, tau
    x = c * c - (eps * d) ** 2
    inside = x > 0
    r = np.sqrt(np.where(inside, x, 1.0))
    w = lam / beta * np.sinc(lam * d / math.pi) * np.sinh(r) / r * (c / math.sinh(c))
    return np.where(inside, w, 0.0)


class BLFIWindow:
    def __init__(self, t_lo: float, t_hi: float, M: int | None = None, tol: float = 1e-10):
        if not 8 * math.pi <= t_lo < t_hi:
            raise ValueError("need 8π <= t_lo < t_hi")
        if t_hi - t_lo >= _GRID * 2:
            raise ValueError("window must be shorter than 2^17")
        if not 0 < tol < 1:
            raise ValueError("tol must be in (0, 1)")
        self.t0 = float(round((t_lo + t_hi) / 2))
        self.u_lo, self.u_hi = t_lo - self.t0, t_hi - self.t0
        self.n2 = n2 = math.isqrt(int(t_lo / TWO_PI))
        self.n = n = math.isqrt(int(t_hi / TWO_PI))
        if M is None:
            M = math.ceil(t_lo**0.25)
        if M < 1:
            raise ValueError("M must be positive")
        self.M = M = min(M, n2)
        self.n_ref, self.theta_base = theta_reference(self.t0)
        self.phase = log_phases(self.t0, n)
        self.v, K = subdivide(n2, M)
        self.head = np.arange(1, M)
        self.extra = np.arange(n2 + 1, n + 1)
        self.head_log, self.v_log, self.extra_log = (
            log_split(self.head),
            log_split(self.v),
            log_split(self.extra),
        )
        self.tau = 1 / (2 * M)
        self.c = math.log(24 * math.sqrt(n2) / tol)
        self.reach = self.c / self.tau
        h = math.pi / (3 * self.tau)
        j = np.arange(
            math.floor((self.u_lo - self.reach) / h), math.ceil((self.u_hi + self.reach) / h) + 1
        )
        self.grid = j * h
        self.table = self._tabulate(K)

    def _tabulate(self, K: np.ndarray) -> np.ndarray:
        first = np.repeat(self.v, K)
        n = np.arange(self.M, self.n2 + 1)
        amp = cis(self.phase[n - 1] - self.phase[first - 1]) / np.sqrt(n)
        freq = np.log1p((n - first) / first) - self.tau
        cuts = np.concatenate(([0], np.cumsum(K)[:-1]))
        return np.stack([np.add.reduceat(amp * cis(u * freq), cuts) for u in self.grid], axis=1)

    def _initial(self, u: np.ndarray) -> np.ndarray:
        k = self.head
        return (spin(self.phase[k - 1], u[:, None], *self.head_log) / np.sqrt(k)).sum(axis=1)

    def _blocks(self, u: np.ndarray) -> np.ndarray:
        lo = np.searchsorted(self.grid, u[0] - self.reach)
        hi = np.searchsorted(self.grid, u[-1] + self.reach, side="right")
        A = spin(self.phase[self.v - 1], u[:, None], *self.v_log)
        W = kernel(u[:, None] - self.grid[lo:hi], self.tau, self.c)
        return ((A @ self.table[:, lo:hi]) * W).sum(axis=1) * cis(self.tau * u)

    def main_sum(self, u: ArrayLike) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, dtype=float))
        if np.any((u < self.u_lo) | (u > self.u_hi)):
            raise ValueError("u is outside the window")
        order = np.argsort(u)
        us = u[order]
        out = np.empty(len(u), dtype=complex)
        rows = max(1, _CHUNK // max(len(self.v), self.M))
        for i in range(0, len(us), rows):
            x = us[i : i + rows]
            out[order[i : i + rows]] = self._initial(x) + self._blocks(x)
        N = np.sqrt((self.t0 + u) / TWO_PI).astype(np.int64)
        for k, hi, lo in zip(self.extra, *self.extra_log, strict=True):
            more = N >= k
            out[more] += spin(self.phase[k - 1], u[more], hi, lo) / math.sqrt(k)
        return out

    def _theta(self, u: np.ndarray) -> np.ndarray:
        return theta_shift(self.t0, self.theta_base, u) + math.pi * (self.n_ref % 2)

    def z(self, u: ArrayLike) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, dtype=float))
        main = 2 * (np.exp(-1j * self._theta(u)) * self.main_sum(u)).real
        return main + rs_remainder(self.t0 + u)

    def zeta(self, u: ArrayLike) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, dtype=float))
        return np.exp(-1j * self._theta(u)) * self.z(u)


def compare(T: float, points: int | None = None) -> dict[str, tuple[float, float]]:
    L = T**0.25
    P = points or int(L)
    out = {}
    for name, build in (("blfi", BLFIWindow), ("nufft", ZWindow)):
        start = time.perf_counter()
        win = build(T, T + L)
        setup = time.perf_counter() - start
        u = np.linspace(win.u_lo, win.u_hi, P)
        win.z(u[:2])
        start = time.perf_counter()
        win.z(u)
        out[name] = (setup, (time.perf_counter() - start) / P)
    few = win.t0 + u[: min(P, 40)]
    start = time.perf_counter()
    z_rs(few)
    out["direct"] = (0.0, (time.perf_counter() - start) / len(few))
    return out
