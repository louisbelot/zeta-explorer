from __future__ import annotations

import functools
import math

from zeta_explorer.cubic import OMEGA, Eis, logderiv, newton
from zeta_explorer.funcfield import (
    Curve,
)
from zeta_explorer.gf import (
    Poly,
    deg,
    minimal_polys,
    pmul,
    primes_dividing,
)
from zeta_explorer.hyperelliptic import irreducible_counts


def _squarefree(p: int, N: int, primes: list[Poly]) -> dict[int, list[tuple[Poly, frozenset]]]:
    out: dict[int, list[tuple[Poly, frozenset]]] = {}

    def grow(f: Poly, used: frozenset, start: int):
        out.setdefault(deg(f), []).append((f, used))
        for i in range(start, len(primes)):
            if deg(f) + deg(primes[i]) > N:
                break
            grow(pmul(f, primes[i], p), used | {i}, i + 1)

    grow((1,), frozenset(), 0)
    return out


def cubefree(p: int, N: int, avoid: Poly) -> list[Poly]:
    bad = set(primes_dividing(avoid, p))
    primes = [P for d in range(1, N + 1) for P in sorted(minimal_polys(p, d)[1])
              if P not in bad]
    sf = _squarefree(p, N, primes)
    out = []
    for a in range(N + 1):
        if (2 * N - a) % 3:
            continue
        for D1, s1 in sf.get(a, []):
            for D2, s2 in sf.get(N - a, []):
                if not s1 & s2:
                    out.append(pmul(D1, pmul(D2, D2, p), p))
    return out


def _prime_sums(q: int, N: int, d: int, L: list[Eis] | None,
                bad: list[tuple[int, int]]) -> list[list[Eis]]:
    if L is None:
        counts = irreducible_counts(q, N)
        pi = [[Eis(c - (r == d)) for r, c in enumerate(counts)] for _ in range(3)]
    else:
        S = logderiv(L, N)
        lam = [[Eis(q**m - d * (m % d == 0)) for m in range(1, N + 1)], S, [s.conj() for s in S]]
        pi = [[Eis(0)] * (N + 1) for _ in range(3)]
        for m in range(1, N + 1):
            for s in range(3):
                acc = lam[s][m - 1]
                for r in range(1, m):
                    if m % r == 0:
                        acc = acc - pi[s * (m // r) % 3][r] * r
                pi[s][m] = acc.div(m)
    for r, e in bad:
        if r <= N:
            for s in range(3):
                pi[s][r] = pi[s][r] - OMEGA[s * e % 3]
    return pi


def weight(q: int, N: int, d: int, L: list[Eis] | None, bad: list[tuple[int, int]]) -> int:
    pi = _prime_sums(q, N, d, L, bad)
    acc = Eis(0)
    for t in range(3):
        lam = []
        for m in range(1, N + 1):
            v = Eis(0)
            for r in range(1, m + 1):
                if m % r:
                    continue
                k = m // r
                inner = sum((pi[(2 * i - k) % 3][r].rot(t * r * i) * math.comb(k, i)
                             for i in range(k + 1)), Eis(0))
                v = v + inner * (r * (-1) ** (k + 1))
            lam.append(v)
        acc = acc + newton(lam, N)[N].rot(-2 * N * t)
    W = acc.div(3)
    if W.b:
        raise ArithmeticError("weight is not rational")
    return W.a


@functools.cache
def _bad(E: Curve) -> tuple[Poly, ...]:
    return tuple(primes_dividing(E.disc, E.p))


def family_size(E: Curve, N: int) -> int:
    return weight(E.p, N, 0, None, [(deg(R), 0) for R in _bad(E)])
