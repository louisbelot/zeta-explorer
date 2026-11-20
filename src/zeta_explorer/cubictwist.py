from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.cubic import OMEGA, Eis, char_sums, index, logderiv, lpolys, newton
from zeta_explorer.funcfield import (
    Curve,
    Twist,
    from_power_sums,
    j0_traces,
    trace_inf,
    traces,
)
from zeta_explorer.gf import (
    Poly,
    deg,
    evaluate,
    field,
    minimal_polys,
    pmul,
    primes_dividing,
)
from zeta_explorer.hyperelliptic import irreducible_counts
from zeta_explorer.onelevel import ADD, inf_trace, local_trace, places


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


@functools.cache
def _principal(E: Curve, N: int, d: int, P_bad: bool) -> int:
    if P_bad:
        return family_size(E, N)
    return weight(E.p, N, d, None, [(deg(R), 0) for R in _bad(E)])


@functools.cache
def _character_weights(E: Curve, N: int, d: int) -> list[int]:
    q, pl = E.p, places(E, d)
    F = field(q, d)
    alpha = minimal_polys(q, d)[0]
    idx = index(F)
    sym = {R: idx[F.polyval(R, alpha)] for R in _bad(E)}
    out = []
    for i, L in enumerate(lpolys(q, pl.P)):
        bad = [(deg(R), int(sym[R][i])) for R in _bad(E) if R != pl.P[i]]
        out.append(weight(q, N, d, L, bad))
    return out


def _divisors(m: int) -> list[int]:
    return [d for d in range(1, m + 1) if m % d == 0]


def exact_traces(E: Curve, N: int, n: int, offdiagonal: bool = True) -> list[Fraction]:
    q, total = E.p, family_size(E, N)
    bad = set(_bad(E))
    out = []
    for m in range(1, n + 1):
        s = inf_trace(E, m) * total
        for d in _divisors(m):
            pl, k = places(E, d), m // d
            if k % 3 and not offdiagonal:
                continue
            w = (_character_weights(E, N, d) if k % 3 else
                 [_principal(E, N, d, P in bad) for P in pl.P])
            s += d * sum(local_trace(int(a), int(t), q**d, k) * c
                         for a, t, c in zip(pl.a, pl.kind, w, strict=True) if t != ADD)
        out.append(Fraction(-s, q**m * total))
    return out


def cubic_twist(E: Curve, D: Poly) -> Curve:
    if E.a2 or E.a4:
        raise ValueError("cubic twists need y^2 = x^3 + B")
    return Curve(E.p, (), (), pmul(E.a6, pmul(D, D, E.p), E.p))


def cubic_sums(E: Curve, Ds: list[Poly], m: int) -> np.ndarray:
    F = field(E.p, m)
    reps, inv = F.orbits
    size = np.bincount(inv)
    a = j0_traces(E.p, m)
    B = F.polyval(E.a6, reps)
    out = np.empty(len(Ds), dtype=np.int64)
    chunk = max(1, 2**22 // len(reps))
    for i in range(0, len(Ds), chunk):
        Dx = evaluate(F, Ds[i:i + chunk], reps)
        out[i:i + chunk] = a[F.mul(B, F.mul(Dx, Dx))] @ size
    return out + np.array([trace_inf(cubic_twist(E, D), m) for D in Ds])


def orthogonal_family(E: Curve, N: int) -> list[Twist]:
    Ds = cubefree(E.p, N, E.disc)
    M = cubic_twist(E, Ds[0]).degree
    S: list[list[int]] = [[] for _ in Ds]
    out: list = [None] * len(Ds)
    todo, m = list(range(len(Ds))), 0
    while todo:
        m += 1
        for i, s in zip(todo, cubic_sums(E, [Ds[i] for i in todo], m), strict=True):
            S[i].append(int(s))
        if m > M // 2:
            for i in todo:
                out[i] = from_power_sums(S[i], M, E.p)
            todo = [i for i in todo if out[i] is None]
    return [Twist(D, L) for D, L in zip(Ds, out, strict=True)]


@dataclass(frozen=True)
class CharTwist:
    D: Poly
    q: int
    c: tuple[Eis, ...]

    def roots(self) -> np.ndarray:
        u = np.array([complex(ck) / self.q**k for k, ck in enumerate(self.c)])
        return 1 / np.roots(u[::-1])

    def traces(self, n: int) -> list[Eis]:
        return logderiv(list(self.c), n)


def unitary_sums(E: Curve, Ds: list[Poly], m: int) -> list[Eis]:
    inf = trace_inf(E, m)
    return [s + inf for s in char_sums(E.p, m, Ds, traces(E, m))]


def unitary_family(E: Curve, N: int) -> list[CharTwist]:
    Ds = cubefree(E.p, N, E.disc)
    M = E.degree + 2 * N
    S = [unitary_sums(E, Ds, m) for m in range(1, M + 1)]
    return [CharTwist(D, E.p, tuple(newton([S[m][i] for m in range(M)], M)))
            for i, D in enumerate(Ds)]


def mean_traces(sums: list[list[int | Eis]], q: int) -> list[Fraction]:
    out = []
    for m, col in enumerate(sums, 1):
        s = sum(col, Eis(0))
        if s.b:
            raise ArithmeticError("family mean is not real")
        out.append(Fraction(-s.a, q**m * len(col)))
    return out

