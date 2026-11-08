from __future__ import annotations

import math
from functools import partial

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.multieval import ZWindow
from zeta_explorer.riemann_siegel import TWO_PI
from zeta_explorer.values import INV_PHI


def bessel_k(nu: float, x: ArrayLike, h: float = 0.05) -> np.ndarray:
    x = np.atleast_1d(np.asarray(x, dtype=float))
    s_max = math.acosh(800 / float(x.min()) + 1)
    s = np.arange(0, s_max + h, h)
    w = np.full(len(s), h)
    w[0] = h / 2
    out = np.empty(len(x))
    for i in range(0, len(x), 4096):
        xs = x[i : i + 4096, None]
        out[i : i + 4096] = (np.exp(-xs * np.cosh(s)) * np.cosh(nu * s) * w).sum(axis=1)
    return out


def fhk_density(y: ArrayLike) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    return 4 * np.exp(-2 * y) * bessel_k(0, 2 * np.exp(-y)).reshape(y.shape)


def fhk_cdf(y: ArrayLike) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    return 2 * np.exp(-y) * bessel_k(1, 2 * np.exp(-y)).reshape(y.shape)


def fhk_centre(N: ArrayLike) -> np.ndarray:
    N = np.asarray(N, dtype=float)
    return np.log(N) - 0.75 * np.log(np.log(N))


def _golden_max(f, a: np.ndarray, b: np.ndarray, iters: int) -> tuple[np.ndarray, np.ndarray]:
    c = b - INV_PHI * (b - a)
    d = a + INV_PHI * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        left = fc > fd
        b = np.where(left, d, b)
        a = np.where(left, a, c)
        c, d = np.where(left, b - INV_PHI * (b - a), d), np.where(left, c, a + INV_PHI * (b - a))
        f_new = f(np.where(left, c, d))
        fc, fd = np.where(left, f_new, fd), np.where(left, fc, f_new)
    x = (a + b) / 2
    return x, f(x)


def _refine(f, a: np.ndarray, b: np.ndarray, va: np.ndarray, vb: np.ndarray, k: np.ndarray,
            best: np.ndarray, margin: float, iters: int) -> np.ndarray:
    pick = np.maximum(va, vb) > best[k] - margin
    _, fx = _golden_max(f, a[pick], b[pick], iters)
    np.maximum.at(best, k[pick], fx)
    return best


def _log_z(win: ZWindow, t: np.ndarray) -> np.ndarray:
    return np.log(np.abs(win.z(t)))


def zeta_maxima(
    T: float,
    count: int,
    length: float = TWO_PI,
    per_zero: int = 16,
    margin: float = 0.02,
    chunk: int = 4000,
    iters: int = 32,
) -> np.ndarray:
    out = []
    for c0 in range(0, count, chunk):
        m = min(chunk, count - c0)
        lo = T + length * c0
        hi = lo + length * m
        win = ZWindow(lo - 1, hi + 1)
        steps = math.ceil(length * math.log(hi / TWO_PI) / TWO_PI * per_zero)
        u = lo - win.t0 + length * np.arange(m * steps + 1) / steps
        v = np.log(np.abs(win.z(u)))
        best = np.maximum(v[:-1].reshape(m, steps).max(axis=1), v[steps::steps])
        k = np.arange(m * steps) // steps
        out.append(_refine(partial(_log_z, win), u[:-1], u[1:], v[:-1], v[1:], k, best,
                           margin, iters))
    return np.concatenate(out)
