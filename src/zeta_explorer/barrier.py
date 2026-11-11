from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

PRIMES = (2, 3, 5)
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


def mollifier_terms(t: float) -> list[tuple[int, float]]:
    out = []
    for mask in range(1 << len(PRIMES)):
        ps = [p for i, p in enumerate(PRIMES) if mask >> i & 1]
        out.append((math.prod(ps), math.prod(-math.exp(log_b(p, t)) for p in ps)))
    return sorted(out)


def mollifier_range(sigma: float, t: float) -> tuple[float, float, float]:
    r = [math.exp(log_b(p, t)) / p**sigma for p in PRIMES]
    return (math.prod(1 - v for v in r), math.prod(1 + v for v in r),
            sum(math.asin(v) for v in r))


def coefficients(n: np.ndarray, cap: ArrayLike, t: float,
                 y: float) -> tuple[np.ndarray, np.ndarray]:
    beta = np.zeros(len(n))
    lam = np.zeros(len(n))
    cap = np.broadcast_to(cap, n.shape)
    for d, a in mollifier_terms(t):
        m = (n % d == 0) & (n <= d * cap)
        q = (n[m] // d).astype(float)
        bq = a * np.exp(log_b(q, t))
        beta[m] += bq
        lam[m] += bq * q**y
    return beta, lam


def lower_edge(N_lo: int, N_hi: int, t: float = 0.2, y: float = 0.2, improved: bool = True,
               every_n: bool = True, c: float | None = None) -> float:
    c = gamma_constant(N_lo, y, t, N_hi) if c is None else c
    r = c * N_lo**-y
    k = (1 - r) / (1 + r)
    sigma = sigma_floor(N_lo, y, t)
    D = math.prod(PRIMES)
    divisors = [d for d, _ in mollifier_terms(t)]
    total = 0.0
    for n in _ranges(2, D * N_hi):
        caps = [float(N_lo), *(n / d for d in divisors)] if every_n else [float(N_hi)]
        worst = np.zeros(len(n))
        for cap in caps:
            beta, lam = coefficients(n, cap, t, y)
            if improved:
                term = np.maximum(np.abs(beta - r * lam), k * np.abs(beta + r * lam))
            else:
                term = np.abs(beta) + r * np.abs(lam)
            if every_n:
                term = np.where((cap >= N_lo) & (cap <= N_hi), term, 0.0)
            worst = np.maximum(worst, term)
        total += float((worst * np.exp(-sigma * np.log(n))).sum())
    return 1 - r - total


def polar_distance(mag: float, arg: float) -> float:
    if arg >= math.pi:
        return 0.0
    return mag if arg <= math.pi / 2 else mag * math.sin(arg)
