from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.explicit import _CHUNK


def bessel_j0(z: ArrayLike) -> np.ndarray:
    z = np.abs(np.asarray(z, dtype=float))
    top = float(z.max(initial=0.0))
    m = math.ceil((top + 10 * top ** (1 / 3) + 30) / 2)
    theta = np.pi * (np.arange(m) + 0.5) / m
    flat = z.ravel()
    out = np.empty_like(flat)
    rows = max(1, _CHUNK // m)
    for i in range(0, len(flat), rows):
        s = slice(i, i + rows)
        out[s] = np.cos(flat[s, None] * np.sin(theta)).mean(axis=1)
    return out.reshape(z.shape)


def square_roots(q: int, a: int) -> int:
    x = np.arange(q, dtype=np.int64)
    return int(np.count_nonzero(x * x % q == a % q))