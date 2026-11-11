from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

_CHUNK = 1 << 22


def _ranges(lo: int, hi: int):
    for a in range(lo, hi + 1, _CHUNK):
        yield np.arange(a, min(hi, a + _CHUNK - 1) + 1, dtype=np.int64)


def log_b(n: ArrayLike, t: float) -> np.ndarray:
    return t / 4 * np.log(n) ** 2


def bsum(M: int, sigma: float, t: float) -> float:
    return sum(float(np.exp(log_b(n, t) - sigma * np.log(n)).sum()) for n in _ranges(1, M))


def lemma_tail(M: int, N: float, sigma: float, t: float) -> float:
    if sigma <= t / 2 * math.log(N):
        raise ValueError("Lemma 8.2 needs sigma > (t/2) log N")
    ends = [(1 - sigma) * math.log(v) + t / 4 * math.log(v) ** 2 for v in (M, N)]
    return math.exp(max(ends)) * math.log(N / M)


def erfi_tail(M: int, N: float, sigma: float, t: float) -> float:
    with mpmath.workdps(30):
        r = mpmath.sqrt(t)

        def g(v):
            return mpmath.erfi(r / 2 * mpmath.log(v) - (sigma - 1) / r)

        return float(mpmath.sqrt(mpmath.pi / t) * mpmath.exp(-((sigma - 1) ** 2) / t)
                     * (g(N) - g(M)))


def x_of(N: float, t: float) -> float:
    return 4 * math.pi * (N * N - t / 16)


def sigma_floor(N: float, y: float, t: float) -> float:
    xN = x_of(N, t)
    return (1 + y) / 2 + t / 4 * math.log(xN / (4 * math.pi)) - t / (2 * xN * xN)


def gamma_constant(N: float, y: float, t: float, n_max: float) -> float:
    return (math.exp(0.02 * y) * (1 - t / (16 * N * N)) ** (-y / 2)
            * n_max ** (t / (2 * (x_of(N, t) - 6))))


def far_bound(N: float, t: float, y: float, M: int) -> tuple[float, float]:
    a = sigma_floor(N, y, t)
    c = gamma_constant(N, y, t, N) * N**-y
    A = bsum(M, a, t) + c * bsum(M, a - y, t)
    B = lemma_tail(M, N, a, t) + c * lemma_tail(M, N, a - y, t)
    return A, B


def far_sup(N_lo: float, N_hi: float, t: float, y: float, M: int, samples: int = 400) -> float:
    return max(sum(far_bound(N, t, y, M)) for N in np.geomspace(N_lo, N_hi, samples))


def top_bound(N0: int, N1: int, t: float, y: float = 1.0) -> float:
    a = sigma_floor(N0, y, t)
    c = gamma_constant(N0, y, t, N1)
    head = bsum(N1, a, t)
    return head - 1 + c * (bsum(N0, a - y, t) * N0**-y + head - bsum(N0, a, t))
