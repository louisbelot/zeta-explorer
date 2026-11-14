from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike


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


def kernel(d: ArrayLike, tau: float, c: float) -> np.ndarray:
    d = np.asarray(d, dtype=float)
    beta, lam, eps = 3 * tau, 2 * tau, tau
    x = c * c - (eps * d) ** 2
    inside = x > 0
    r = np.sqrt(np.where(inside, x, 1.0))
    w = lam / beta * np.sinc(lam * d / math.pi) * np.sinh(r) / r * (c / math.sinh(c))
    return np.where(inside, w, 0.0)
