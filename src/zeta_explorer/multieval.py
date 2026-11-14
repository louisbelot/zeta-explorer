from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.riemann_siegel import TWO_PI, rs_remainder, theta_prime, theta_tail
from zeta_explorer.spacing import gram_block_zeros
from zeta_explorer.zeros import gram_points

_CHUNK = 1 << 20


def nufft(x: ArrayLike, c: ArrayLike, m: int, tol: float = 1e-13) -> np.ndarray:
    if m < 2 or m % 2:
        raise ValueError("m must be even and positive")
    x = np.mod(np.asarray(x, dtype=float), TWO_PI)
    c = np.asarray(c, dtype=complex)
    mr = 2 * m
    sp = math.ceil(1.5 * math.log(1 / tol) / math.pi)
    tau = math.pi * sp / (3 * m * m)
    h = TWO_PI / mr
    offsets = np.arange(-sp, sp + 1)
    grid = np.zeros(mr, dtype=complex)
    rows = max(1, _CHUNK // len(offsets))
    for i in range(0, len(x), rows):
        xs, cs = x[i : i + rows], c[i : i + rows]
        idx = np.rint(xs / h).astype(np.int64)[:, None] + offsets
        w = (np.exp(-((idx * h - xs[:, None]) ** 2) / (4 * tau)) * cs[:, None]).ravel()
        idx = (idx % mr).ravel()
        grid += np.bincount(idx, w.real, mr) + 1j * np.bincount(idx, w.imag, mr)
    k = np.arange(-m // 2, m // 2)
    return np.fft.fft(grid)[k % mr] * np.exp(tau * k**2) / (mr * math.sqrt(tau / math.pi))


def fft_size(n: int) -> int:
    best = 2 ** math.ceil(math.log2(n))
    p5 = 1
    while p5 < best:
        p35 = p5
        while p35 < best:
            p = p35
            while p < n:
                p *= 2
            best = min(best, p)
            p35 *= 3
        p5 *= 5
    return best


class BandSum:
    def __init__(self, a: ArrayLike, u_lo: float, u_hi: float, oversample: float = 2.0,
                 tol: float = 1e-12):
        if oversample <= 1:
            raise ValueError("oversample must exceed 1")
        a = np.asarray(a, dtype=complex)
        n = len(a)
        self.u_lo, self.u_hi = u_lo, u_hi
        self.log_n = np.log(np.arange(1, n + 1, dtype=float))
        self.beta = math.log(n) / 2
        alpha = oversample * self.beta
        self.delta = math.pi / alpha
        self.taps = math.ceil(2 * math.log(1 / tol) / (math.pi * (1 - 1 / oversample)))
        self.eps2 = (alpha - self.beta) / (self.taps * self.delta)
        half = math.ceil(max(-u_lo, u_hi) / self.delta) + self.taps + 1
        self.m = 2 * fft_size(half)
        j = np.arange(-self.m // 2, self.m // 2)
        F = nufft(self.delta * self.log_n, a, self.m, tol / 10)
        self.samples = F * np.exp(1j * self.beta * self.delta * j)

    def __call__(self, u: ArrayLike) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, dtype=float))
        if np.any((u < self.u_lo) | (u > self.u_hi)):
            raise ValueError("u is outside the window")
        out = np.empty(len(u), dtype=complex)
        taps = np.arange(-self.taps, self.taps + 1)
        rows = max(1, _CHUNK // len(taps))
        for i in range(0, len(u), rows):
            us = u[i : i + rows]
            js = np.rint(us / self.delta).astype(np.int64)[:, None] + taps
            d = us[:, None] - js * self.delta
            w = np.sinc(d / self.delta) * np.exp(-self.eps2 * d * d / 2)
            out[i : i + rows] = (self.samples[js + self.m // 2] * w).sum(axis=1)
        return out * np.exp(-1j * self.beta * u)


def smallest_factor(n: int) -> np.ndarray:
    spf = np.zeros(n + 1, dtype=np.int64)
    for p in range(2, math.isqrt(n) + 1):
        if spf[p] == 0:
            multiples = spf[p * p :: p]
            multiples[multiples == 0] = p
    ks = np.arange(n + 1)
    spf[spf == 0] = ks[spf == 0]
    return spf


def _two_sum(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)


def _dd_add(ah, al, bh, bl):
    s, e = _two_sum(ah, bh)
    e = e + (al + bl)
    hi = s + e
    return hi, e - (hi - s)


def log_phases(t0: float, n: int) -> np.ndarray:
    spf = smallest_factor(n)
    ps = np.flatnonzero(spf[2:] == np.arange(2, n + 1)) + 2
    hi, lo = np.zeros(n + 1), np.zeros(n + 1)
    with mpmath.workdps(25 + int(math.log10(t0))):
        t = mpmath.mpf(t0)
        two_pi = 2 * mpmath.pi
        for p in ps:
            x = t * mpmath.log(int(p)) % two_pi
            hi[p] = float(x)
            lo[p] = float(x - hi[p])
        tp_hi = float(two_pi)
        tp_lo = float(two_pi - tp_hi)
    ks = np.arange(n + 1)
    q = ks // np.maximum(spf, 1)
    done = np.zeros(n + 1, dtype=bool)
    done[:2] = True
    done[ps] = True
    while not done.all():
        k = np.flatnonzero(~done & done[q])
        p = spf[k]
        s_hi, s_lo = _dd_add(hi[p], lo[p], hi[q[k]], lo[q[k]])
        wrap = (s_hi > tp_hi) | ((s_hi == tp_hi) & (s_lo >= tp_lo))
        w_hi, w_lo = _dd_add(s_hi, s_lo, np.full(len(k), -tp_hi), np.full(len(k), -tp_lo))
        hi[k], lo[k] = np.where(wrap, w_hi, s_hi), np.where(wrap, w_lo, s_lo)
        done[k] = True
    return hi[1:]


def theta_reference(t0: float) -> tuple[int, float]:
    with mpmath.workdps(25 + int(math.log10(t0))):
        th = mpmath.siegeltheta(mpmath.mpf(t0))
        n = int(mpmath.nint(th / mpmath.pi))
        return n, float(th - n * mpmath.pi)


def theta_shift(t0: float, base: float, u: ArrayLike) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    t = t0 + u
    return (
        base
        + u / 2 * math.log(t0 / TWO_PI)
        + t / 2 * np.log1p(u / t0)
        - u / 2
        + theta_tail(t)
        - theta_tail(np.float64(t0))
    )


class ZWindow:
    def __init__(self, t_lo: float, t_hi: float, oversample: float = 2.0, tol: float = 1e-12):
        if not 8 * math.pi <= t_lo < t_hi:
            raise ValueError("need 8π <= t_lo < t_hi")
        self.t0 = float(round((t_lo + t_hi) / 2))
        self.u_lo, self.u_hi = t_lo - self.t0, t_hi - self.t0
        self.n_min = math.isqrt(int(t_lo / TWO_PI))
        self.n = n = math.isqrt(int(t_hi / TWO_PI))
        self.n_ref, self.theta_base = theta_reference(self.t0)
        ks = np.arange(1, n + 1)
        self.log_n = np.log(ks)
        self.b = np.exp(-1j * log_phases(self.t0, n)) / np.sqrt(ks)
        self.sum = BandSum(self.b, self.u_lo, self.u_hi, oversample, tol)

    def theta_offset(self, u: ArrayLike) -> np.ndarray:
        return theta_shift(self.t0, self.theta_base, u)

    def main_sum(self, u: ArrayLike) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, dtype=float))
        out = self.sum(u)
        if self.n_min < self.n:
            N = np.sqrt((self.t0 + u) / TWO_PI).astype(np.int64)
            for k in range(self.n_min + 1, self.n + 1):
                extra = N < k
                out[extra] -= self.b[k - 1] * np.exp(-1j * u[extra] * self.log_n[k - 1])
        return out

    def z(self, u: ArrayLike) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, dtype=float))
        th = self.theta_offset(u) + math.pi * (self.n_ref % 2)
        main = 2 * (np.exp(1j * th) * self.main_sum(u)).real
        return main + rs_remainder(self.t0 + u)

    def unfold(self, u: ArrayLike) -> np.ndarray:
        return self.theta_offset(u) / math.pi

    def gram(self, ns: ArrayLike, tol: float = 1e-13) -> np.ndarray:
        ns = np.asarray(ns)
        target = (ns - self.n_ref) * math.pi
        u = (target - self.theta_base) / theta_prime(self.t0)
        for _ in range(50):
            du = (self.theta_offset(u) - target) / theta_prime(self.t0 + u)
            u -= du
            if np.all(np.abs(du) < tol * max(1.0, float(np.max(np.abs(u))))):
                return u
        raise RuntimeError("Newton didn't converge for some Gram points")


def zeros_fast(
    n_lo: int, n_hi: int, atol: float = 1e-10, **window
) -> tuple[ZWindow, np.ndarray]:
    g = gram_points([n_lo, n_hi])
    win = ZWindow(g[0] - 1, g[1] + 1, **window)
    u = gram_block_zeros(win.gram(np.arange(n_lo, n_hi + 1)), n_lo, win.z, win.z, atol=atol)
    return win, u
