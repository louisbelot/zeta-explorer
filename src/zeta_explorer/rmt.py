from __future__ import annotations

import functools
import math
from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike

Kernel = Callable[[np.ndarray, np.ndarray], np.ndarray]


def sine(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sinc(x - y)


def symplectic(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sinc(x - y) - np.sinc(x + y)


def even_orthogonal(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sinc(x - y) + np.sinc(x + y)


def _dirichlet(z: np.ndarray, m: float, n: float) -> np.ndarray:
    a = math.pi * z / (2 * n)
    den = np.sin(a)
    small = np.abs(den) < 1e-14
    out = np.sin(m * a) / (2 * n * np.where(small, 1.0, den))
    return np.where(small, m / (2 * n) * np.cos(m * a) / np.cos(a), out)


def symplectic_n(n: float) -> Kernel:
    if n < 1:
        raise ValueError("n must be >= 1")
    m = 2 * n + 1
    return lambda x, y: _dirichlet(x - y, m, n) - _dirichlet(x + y, m, n)


def even_orthogonal_n(n: float) -> Kernel:
    if n < 1:
        raise ValueError("n must be >= 1")
    m = 2 * n - 1
    return lambda x, y: _dirichlet(x - y, m, n) + _dirichlet(x + y, m, n)


@functools.cache
def _legendre(m: int) -> tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(m)


def _nystrom(
    kernel: Kernel, a: float, b: float, m: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, w = _legendre(m)
    x = (b - a) / 2 * u + (a + b) / 2
    sw = np.sqrt(w * (b - a) / 2)
    return x, sw, np.eye(m) - sw[:, None] * kernel(x[:, None], x[None, :]) * sw[None, :]


def fredholm_det(kernel: Kernel, a: float, b: float, m: int = 40) -> float:
    return float(np.linalg.det(_nystrom(kernel, a, b, m)[2]))


def _gap_and_density(kernel: Kernel, s: float, m: int) -> tuple[float, float]:
    if s == 0:
        return 1.0, float(kernel(np.zeros(1), np.zeros(1))[0])
    x, sw, I_A = _nystrom(kernel, 0.0, s, m)
    edge = kernel(x, np.full(m, s))
    rho = np.linalg.solve(I_A, sw * edge)
    resolvent = float(kernel(np.array([s]), np.array([s]))[0]) + sw @ (edge * rho)
    E = float(np.linalg.det(I_A))
    return E, E * resolvent


def gap(kernel: Kernel, s: ArrayLike, m: int = 40) -> float | np.ndarray:
    ss = np.atleast_1d(np.asarray(s, dtype=float))
    out = np.array([_gap_and_density(kernel, si, m)[0] for si in ss])
    return float(out[0]) if np.ndim(s) == 0 else out


def first_density(kernel: Kernel, s: ArrayLike, m: int = 40) -> float | np.ndarray:
    ss = np.atleast_1d(np.asarray(s, dtype=float))
    out = np.array([_gap_and_density(kernel, si, m)[1] for si in ss])
    return float(out[0]) if np.ndim(s) == 0 else out


def fit_symplectic_n(
    first: ArrayLike, grid: ArrayLike, m: int = 40
) -> tuple[float, float, float, np.ndarray]:
    x = np.asarray(first, dtype=float)
    ns = np.asarray(grid, dtype=float)
    s = np.linspace(0, x.max() * 1.001, 200)
    loglik = np.full(len(ns), -np.inf)
    for i, n in enumerate(ns):
        if n <= s[-1]:
            continue
        p = np.interp(x, s, first_density(symplectic_n(n), s, m))
        if np.all(p > 0):
            loglik[i] = np.sum(np.log(p))
    best = int(np.argmax(loglik))
    inside = ns[loglik >= loglik[best] - 1.92]
    return float(ns[best]), float(inside.min()), float(inside.max()), loglik
