from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from zeta_explorer.gf import Poly, field, pmul


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
