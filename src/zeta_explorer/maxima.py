from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike


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
