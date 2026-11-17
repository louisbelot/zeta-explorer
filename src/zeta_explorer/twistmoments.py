from __future__ import annotations

from fractions import Fraction

import numpy as np

from zeta_explorer.funcfield import Twist, legendre, traces
from zeta_explorer.gf import field
from zeta_explorer.rmt import so_moment


def place_traces(q: int, d: int) -> np.ndarray:
    F = field(q, d)
    t = F.orbits[0]
    keep = (F.degrees(t) == d) & (t != 0) & (t != 1)
    return traces(legendre(q), d)[t[keep]].astype(np.int64)


def _mul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    n, m = A.shape
    out = np.zeros_like(A)
    for i, j in zip(*np.nonzero(A), strict=True):
        out[i:, j:] += A[i, j] * B[: n - i, : m - j]
    return out


def _series(n: int) -> np.ndarray:
    return np.zeros((n + 1, n + 1), dtype=object)


def coefficients(a: int, Q: int, m: int) -> list[int]:
    b = [1, a]
    while len(b) <= m:
        b.append(a * b[-1] - Q * b[-2])
    return b[: m + 1]


def squarefree_counts(q: int, n: int) -> list[int]:
    S = [1, q] + [q**k - q ** (k - 1) for k in range(2, n + 1)]
    for _ in range(2):
        S = [sum((-1) ** (k - j) * S[j] for j in range(k + 1)) for k in range(n + 1)]
    return S


def diagonal(q: int, n: int) -> np.ndarray:
    a0, a1 = (int(x) for x in traces(legendre(q), 1)[[0, 1]])
    G = _series(n)
    for k in range(n // 2 + 1):
        G[2 * k, 0] = (k + 1) + a0 * a1 * k
    S = _series(n)
    S[0, :] = squarefree_counts(q, n)
    G = _mul(G, S)
    for d in range(1, n // 2 + 1):
        for a in place_traces(q, d):
            b = coefficients(int(a), q**d, n // d)
            local = _series(n)
            local[0, 0] = 1
            for k in range(1, n // (2 * d) + 1):
                for j in range(0, n + 1, d):
                    local[2 * k * d, j] = b[2 * k] * (-1) ** (j // d)
            G = _mul(G, local)
    return G[:, n]


def recipe_first(q: int, n: int) -> Fraction:
    if n % 2:
        raise ValueError("even degree only")
    c = diagonal(q, n)
    total = sum((Fraction(int(c[j]), q**j) for j in range(n)), Fraction(0))
    total += Fraction(int(c[n]), 2 * q**n)
    return total / squarefree_counts(q, n)[n]


def first_moment(fam: list[Twist]) -> Fraction:
    return sum((t.L.central for t in fam), Fraction(0)) / len(fam)


def moments(fam: list[Twist], ks: list[int]) -> np.ndarray:
    v = np.array([float(t.L.central) for t in fam if t.L.eps == 1])
    return np.array([np.mean(v**k) for k in ks])


def so_ratios(fam: list[Twist], ks: list[int]) -> np.ndarray:
    return moments(fam, ks) / so_moment(fam[0].L.N // 2, ks)
