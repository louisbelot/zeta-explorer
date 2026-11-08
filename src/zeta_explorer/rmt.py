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


def odd_orthogonal_n(n: float) -> Kernel:
    if n < 1:
        raise ValueError("n must be >= 1")
    m = 2 * n
    return lambda x, y: _dirichlet(x - y, m, n) - _dirichlet(x + y, m, n)


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


def fit_n(
    first: ArrayLike, grid: ArrayLike, family: Callable[[float], Kernel] = symplectic_n,
    m: int = 40,
) -> tuple[float, float, float, np.ndarray]:
    x = np.asarray(first, dtype=float)
    ns = np.asarray(grid, dtype=float)
    s = np.linspace(0, x.max() * 1.001, 200)
    loglik = np.full(len(ns), -np.inf)
    for i, n in enumerate(ns):
        if n <= s[-1]:
            continue
        p = np.interp(x, s, first_density(family(n), s, m))
        if np.all(p > 0):
            loglik[i] = np.sum(np.log(p))
    best = int(np.argmax(loglik))
    inside = ns[loglik >= loglik[best] - 1.92]
    return float(ns[best]), float(inside.min()), float(inside.max()), loglik


def so_angles(
    n: int, size: int, seed: int = 0, batch: int = 10_000, odd: bool = False
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    dim = 2 * n + odd
    out = []
    for k in range(0, size, batch):
        Z = rng.standard_normal((min(batch, size - k), dim, dim))
        Q, R = np.linalg.qr(Z)
        Q *= np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]
        Q[np.linalg.det(Q) < 0, :, 0] *= -1
        theta = np.sort(np.abs(np.angle(np.linalg.eigvals(Q))), axis=1)
        out.append(theta[:, int(odd)::2])
    return np.concatenate(out)


def char_poly_at_one(theta: np.ndarray) -> np.ndarray:
    return np.prod(2 - 2 * np.cos(theta), axis=-1)


def excised_first(
    n: int, cutoff: float, size: int, seed: int = 0, odd: bool = False
) -> np.ndarray:
    theta = so_angles(n, size, seed, odd=odd)
    keep = char_poly_at_one(theta) >= cutoff
    return n * theta[keep, 0] / math.pi

