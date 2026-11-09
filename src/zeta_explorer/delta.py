from __future__ import annotations

import functools
import math

import numpy as np

from zeta_explorer.explicit import primes

KAPPA = 5.5


def _moduli(N: int) -> list[int]:
    small = primes(46341)
    out, M, p = [], 1, (1 << 31) - 1
    while M <= 4 * N**6:
        if all(p % q for q in small):
            out.append(p)
            M *= p
        p -= 2
    return out


def _residues(N: int, ps: list[int]) -> np.ndarray:
    k = np.arange(math.isqrt(2 * N) + 2)
    e = k * (k + 1) // 2
    c = np.where(k % 2, -1, 1) * (2 * k + 1)
    e, c = e[e < N], c[e < N]
    P = np.array(ps, dtype=np.int64)[:, None]
    pair = (e[:, None] + e[None, :]).ravel()
    keep = pair < N
    f = np.zeros(N, dtype=np.int64)
    np.add.at(f, pair[keep], np.outer(c, c).ravel()[keep])
    f = f[None, :] % P
    for _ in range(6):
        g = np.zeros_like(f)
        for s, ck in zip(e, c, strict=True):
            g[:, s:] += ck * f[:, : N - s]
        f = g % P
    return f


def _garner(r: np.ndarray, ps: list[int]) -> np.ndarray:
    v = []
    for j, p in enumerate(ps):
        acc = np.zeros(r.shape[1], dtype=np.int64)
        m = 1
        for i in range(j):
            acc = (acc + v[i] * m) % p
            m = m * ps[i] % p
        v.append((r[j] - acc) * pow(m, -1, p) % p)
    x = v[-1].astype(object)
    for i in range(len(ps) - 2, -1, -1):
        x = x * ps[i] + v[i].astype(object)
    M = math.prod(ps)
    return np.where(x > M // 2, x - M, x)


@functools.lru_cache(maxsize=4)
def _tau(N: int) -> np.ndarray:
    ps = _moduli(N)
    return np.concatenate([[0], _garner(_residues(N, ps), ps)]).astype(object)


def _size(N: int) -> int:
    return 1 << max(10, (N - 1).bit_length())


def tau(N: int) -> np.ndarray:
    if N < 1:
        raise ValueError("need N >= 1")
    return _tau(_size(N))[: N + 1].copy()


@functools.lru_cache(maxsize=4)
def _lam(N: int) -> np.ndarray:
    t = _tau(N)
    n = np.arange(1, N + 1, dtype=float)
    out = np.zeros(N + 1)
    out[1:] = np.array([float(x) for x in t[1:]]) * n**-KAPPA
    return out


def lam(N: int) -> np.ndarray:
    return _lam(_size(N))[: N + 1]
