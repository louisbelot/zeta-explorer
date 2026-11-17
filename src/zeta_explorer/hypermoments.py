from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator

import numpy as np

from zeta_explorer.funcfield import _newton
from zeta_explorer.gf import Poly, character_sums, field, pdivmod, pmul, ppowmod


def irreducibles(q: int, d: int) -> list[Poly]:
    F = field(q, d)
    reps, inv = F.orbits
    alpha = reps[np.bincount(inv) == d]
    c = np.zeros((len(alpha), d + 1), dtype=np.int64)
    c[:, 0] = 1
    beta = alpha
    for _ in range(d):
        shifted = np.zeros_like(c)
        shifted[:, 1:] = c[:, :-1]
        c = F.add(shifted, F.neg(F.mul(c, beta[:, None])))
        beta = F.frobenius(beta)
    if np.any(c >= q):
        raise ArithmeticError("minimal polynomial outside the prime field")
    out = [tuple(int(x) for x in row) for row in c]
    return sorted(out, key=lambda P: P[::-1])


def squarefrees(q: int, n: int) -> Iterator[tuple[Poly, tuple[int, ...]]]:
    primes = [P for d in range(1, n + 1) for P in irreducibles(q, d)]

    def grow(f: Poly, degs: tuple[int, ...], start: int):
        yield f, degs
        for i in range(start, len(primes)):
            P = primes[i]
            if len(f) - 1 + len(P) - 1 > n:
                break
            yield from grow(pmul(f, P, q), degs + (len(P) - 1,), i + 1)

    yield from grow((1,), (), 0)


def squares(q: int, n: int) -> list[tuple[Poly, list[Poly]]]:
    primes = [P for d in range(1, n // 2 + 1) for P in irreducibles(q, d)]
    out = []

    def grow(h: Poly, support: list[Poly], start: int):
        out.append((h, support))
        for i in range(start, len(primes)):
            P = primes[i]
            if 2 * (len(h) + len(P) - 2) > n:
                break
            hk = pmul(h, P, q)
            while 2 * (len(hk) - 1) <= n:
                grow(hk, support + [P], i + 1)
                hk = pmul(hk, P, q)

    grow((1,), [], 0)
    return out


def lpoly(q: int, fs: list[Poly]) -> np.ndarray:
    d = len(fs[0]) - 1
    if d == 0:
        raise ValueError("the trivial character has no L-polynomial")
    N = d - 1 if d % 2 else d - 2
    S = [character_sums(field(q, m), fs, np.bincount(field(q, m).orbits[1]))
         for m in range(1, N // 2 + 1)]
    S = np.array(S, dtype=np.int64).T if S else np.zeros((len(fs), 0), dtype=np.int64)
    s = -1 if (q - 1) // 2 * d % 2 else 1
    out = np.zeros((len(fs), d), dtype=object)
    for i, row in enumerate(S):
        c = _newton([int(x) for x in row])
        if d % 2 == 0:
            c = [sum(c[: k + 1]) for k in range(len(c))]
        full = c + [0] * (N + 1 - len(c))
        for k in range(N // 2):
            full[N - k] = q ** (N // 2 - k) * c[k]
        if d % 2 == 0:
            full = [a - b for a, b in zip(full + [0], [0] + full, strict=True)]
        out[i] = [int(a) * s**k for k, a in enumerate(full)]
    return out


def _symbol(f: Poly, P: Poly, q: int) -> int:
    r = pdivmod(f, P, q)[1]
    if not r:
        return 0
    v = ppowmod(r, (q ** (len(P) - 1) - 1) // 2, P, q)
    return 1 if v == (1,) else -1


def dual_coefficients(q: int, g: int) -> list[int]:
    n = 2 * g + 1
    hs = sorted(squares(q, g), key=lambda h: len(h[0]))
    upto = [sum(2 * (len(h) - 1) <= m for h, _ in hs) for m in range(g + 1)]
    by_deg: dict[int, list[tuple[Poly, tuple[int, ...]]]] = defaultdict(list)
    for f, degs in squarefrees(q, g):
        by_deg[len(f) - 1].append((f, degs))
    T = [0] * (g + 1)
    for d, group in sorted(by_deg.items()):
        Ls = lpoly(q, [f for f, _ in group]) if d else [None]
        s = -1 if (q - 1) // 2 * d % 2 else 1
        for (f, degs), L in zip(group, Ls, strict=True):
            if d == 0:
                A = [1, 0, -q]
                b0 = [q**k for k in range(n + 1)]
            else:
                A = [0] * (d + 2)
                for k, a in enumerate(L):
                    A[k] += a
                    A[k + 2] -= q * a
                b0 = [1] + [0] * n
                for e in degs:
                    for k in range(2 * e, n + 1):
                        b0[k] += b0[k - 2 * e]
            symbols: dict[Poly, int] = {}
            for h, support in hs[: upto[g - d]]:
                b = list(b0)
                for P in support:
                    e = len(P) - 1
                    if P not in symbols:
                        symbols[P] = s**e * _symbol(f, P, q) if d else 1
                    psi = symbols[P]
                    if psi:
                        for k in range(e, n + 1):
                            b[k] -= psi * b[k - e]
                T[d + 2 * (len(h) - 1)] += sum(a * b[n - k] for k, a in enumerate(A))
    return T


def dual_mean(q: int, g: int) -> float:
    T = dual_coefficients(q, g)
    n = 2 * g + 1
    total = sum((1 if j == g else 2) * t * q ** (-j / 2) for j, t in enumerate(T))
    return total / (q**n - q ** (n - 1))
