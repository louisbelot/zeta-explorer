from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

BLOCK = 1 << 22


def lehman_sum(omega: ArrayLike, gamma: ArrayLike, alpha: float | None = None) -> np.ndarray:
    w = np.atleast_1d(np.asarray(omega, dtype=float))
    g = np.asarray(gamma, dtype=float)
    weight = np.ones_like(g) if alpha is None else np.exp(-(g**2) / (2 * alpha))
    c1 = weight / (0.25 + g**2)
    c2 = 2 * g * c1
    out = np.zeros_like(w)
    step = max(1, BLOCK // max(len(w), 1))
    for i in range(0, len(g), step):
        phase = np.outer(w, g[i : i + step])
        out -= np.cos(phase) @ c1[i : i + step] + np.sin(phase) @ c2[i : i + step]
    return out


def lehman_errors(
    omega: float, eta: float, alpha: float, T: float, A: float, improved: bool = True
) -> dict[str, mpmath.mpf]:
    lehman = 4 * A / omega <= alpha <= A**2 and 2 * A / alpha <= eta <= omega / 2
    if not (omega - eta > 1 and lehman):
        raise ValueError("parameters outside Lehman's conditions")
    if not 2 * math.pi * math.e < T <= A:
        raise ValueError("need 2 pi e < T <= A")
    if improved and omega - eta <= 437.5:
        raise ValueError("the improved s1 needs omega - eta > 437.5")
    with mpmath.workdps(30):
        lo, hi = mpmath.mpf(omega) - eta, mpmath.mpf(omega) + eta
        alpha, eta, T, A = (mpmath.mpf(v) for v in (alpha, eta, T, A))
        gauss = mpmath.exp(-alpha * eta**2 / 2)
        return {
            "s1": (mpmath.mpf("2.1457") if improved else mpmath.mpf("3.05")) / lo,
            "s2": 4 * hi * mpmath.exp(-lo / 6),
            "s3": 2 * gauss / (mpmath.sqrt(2 * mpmath.pi * alpha) * eta),
            "s4": mpmath.mpf("0.08") * mpmath.sqrt(alpha) * gauss,
            "s5": mpmath.exp(-(T**2) / (2 * alpha))
            * (alpha / (mpmath.pi * T**2) * mpmath.log(T / (2 * mpmath.pi)) + 8 * mpmath.log(T) / T
               + 4 * alpha / T**3),
            "s6": A * mpmath.log(A) * mpmath.exp(-(A**2) / (2 * alpha) + hi / 2)
            * (4 / mpmath.sqrt(alpha) + 15 * eta),
        }
