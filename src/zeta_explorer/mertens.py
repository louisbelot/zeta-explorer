from __future__ import annotations

import math

import numpy as np


def mobius_sieve(n: int) -> np.ndarray:
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, math.isqrt(n) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    for p in np.flatnonzero(is_prime):
        mu[p::p] *= -1
    for p in np.flatnonzero(is_prime[: math.isqrt(n) + 1]):
        mu[p * p :: p * p] = 0
    return mu


def mertens_table(n: int) -> np.ndarray:
    return np.cumsum(mobius_sieve(n), dtype=np.int64)


def mertens(x: int, u: int | None = None) -> int:
    if x < 1:
        return 0
    if u is None:
        u = max(math.isqrt(x) + 1, int(x ** (2 / 3) / 4))
    small = mertens_table(max(u, math.isqrt(x) + 1))
    if x < len(small):
        return int(small[x])
    K = x // (len(small))
    L = np.zeros(K + 1, dtype=np.int64)
    for k in range(K, 0, -1):
        v = x // k
        r = math.isqrt(v)
        s = v // (r + 1)
        total = 1
        nb = min(r, K // k)
        if nb >= 2:
            total -= int(L[k * np.arange(2, nb + 1)].sum())
        d = np.arange(max(nb, 1) + 1, r + 1, dtype=np.int64)
        total -= int(small[v // d].sum())
        q = np.arange(1, s + 1, dtype=np.int64)
        counts = v // q - v // (q + 1)
        counts[-1] = v // s - r
        total -= int((counts * small[q]).sum())
        L[k] = total
    return int(L[1])
