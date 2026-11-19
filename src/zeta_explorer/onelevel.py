from __future__ import annotations

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.funcfield import Curve, LPoly, Twist, _newton, trace_inf, traces
from zeta_explorer.gf import Poly, deg, field, minimal_polys, pmul, primes_dividing
from zeta_explorer.hypermoments import _symbol, lpoly

GOOD, MULT, ADD = 0, 1, 2


def power_sums(c: list[int] | tuple[int, ...], n: int) -> list[int]:
    c = list(c) + [0] * n
    p: list[int] = []
    for m in range(1, n + 1):
        p.append(-m * c[m] - sum(c[k] * p[m - k - 1] for k in range(1, m)))
    return p


def frob_traces(L: LPoly, n: int) -> list[Fraction]:
    return [Fraction(x, L.q**m) for m, x in enumerate(power_sums(L.c, n), 1)]


def family_traces(fam: list[Twist], n: int, eps: int | None = None) -> list[Fraction]:
    rows = [frob_traces(t.L, n) for t in fam if eps is None or t.L.eps == eps]
    return [sum(col, Fraction(0)) / len(rows) for col in zip(*rows, strict=True)]


def rmt_trace(dim: int, det: int, n: int) -> int:
    if n < 1:
        raise ValueError("n must be >= 1")
    even, low = n % 2 == 0, n <= dim - 2
    if dim % 2 == 0:
        if det == 1:
            return int(even and low)
        return 2 - low if even else 0
    v = 1 - (not even and low)
    return v if det == 1 else (-1) ** n * v


@dataclass(frozen=True)
class Places:
    P: list[Poly]
    a: np.ndarray
    kind: np.ndarray
    chi: dict[Poly, np.ndarray]


@functools.cache
def places(E: Curve, d: int) -> Places:
    q = E.p
    F = field(q, d)
    alpha, P = minimal_polys(q, d)
    bad = F.polyval(E.disc, alpha) == 0
    mult = bad & (F.polyval(E.c4, alpha) != 0)
    kind = np.where(bad, np.where(mult, MULT, ADD), GOOD)
    chi = {Q: F.chi[F.polyval(Q, alpha)].astype(np.int64) for Q in primes_dividing(E.disc, q)}
    return Places(P, traces(E, d)[alpha].astype(np.int64), kind, chi)


def local_trace(a: int, kind: int, qd: int, k: int) -> int:
    if kind == ADD:
        return 0
    if kind == MULT:
        return a**k
    s0, s1 = 2, a
    for _ in range(k - 1):
        s0, s1 = s1, a * s1 - qd * s0
    return s1


def inf_trace(E: Curve, k: int) -> int:
    return local_trace(trace_inf(E, 1), (GOOD, MULT, ADD)[E.f_inf], E.p, k)


def twisted_infinity(E: Curve, N: int) -> Curve:
    return E.twist((0,) * N + (1,))


def multiplicative(E: Curve) -> list[Poly]:
    out = []
    for Q in primes_dividing(E.disc, E.p):
        pl = places(E, deg(Q))
        if pl.kind[pl.P.index(Q)] == MULT:
            out.append(Q)
    return out


def _mul(a: list[int], b: list[int], N: int) -> list[int]:
    out = [0] * (N + 1)
    for i, x in enumerate(a[: N + 1]):
        if x:
            for j, y in enumerate(b[: N + 1 - i]):
                out[i + j] += x * y
    return out


def _divide(a: list[int], e: int, s: int) -> list[int]:
    out = list(a)
    for j in range(e, len(out)):
        out[j] -= s * out[j - e]
    return out


def char_sum(q: int, N: int, L, fdegs: list[int], avoid: list[tuple[int, int]]) -> int:
    if L is None:
        a = [1, q] + [q**j - q ** (j - 1) for j in range(2, N + 1)]
        a = a[: N + 1]
    else:
        a = _mul([int(x) for x in L], [1, 0, -q], N)
        for e in fdegs:
            a = _divide(a, 2 * e, -1)
    for e, s in avoid:
        if s:
            a = _divide(a, e, s)
    return a[N]


def _jacobi(Q: Poly, f: list[Poly], q: int) -> int:
    return math.prod(_symbol(Q, R, q) for R in f)


def _prod(f: list[Poly], q: int) -> Poly:
    return functools.reduce(lambda x, y: pmul(x, y, q), f, (1,))


@dataclass(frozen=True)
class Weights:
    total: int
    even: list[int]
    odd: list[int]


def weights(E: Curve, N: int, d: int, sign: int | None = None) -> Weights:
    q, pl = E.p, places(E, d)
    S = primes_dividing(E.disc, q)
    base = [(deg(Q), 1) for Q in S]
    total = char_sum(q, N, None, [], base)
    coprime = char_sum(q, N, None, [], base + [(d, 1)])
    L = lpoly(q, pl.P)
    even = [coprime if k == GOOD else total for k in pl.kind]
    odd = [0 if k == ADD else char_sum(q, N, L[i], [d], [(deg(Q), int(pl.chi[Q][i])) for Q in S])
           for i, k in enumerate(pl.kind)]
    if sign is None:
        return Weights(total, even, odd)
    M = multiplicative(E)
    if not M:
        raise ValueError("no multiplicative places to split the family by")
    mdeg = [deg(Q) for Q in M]
    LM = lpoly(q, [_prod(M, q)])[0]
    chiM = {Q: _jacobi(Q, M, q) for Q in S}
    avoidM = [(deg(Q), chiM[Q]) for Q in S]
    totM = char_sum(q, N, LM, mdeg, avoidM)
    good = [i for i, k in enumerate(pl.kind) if k == GOOD]
    LPM = lpoly(q, [pmul(pl.P[i], _prod(M, q), q) for i in good]) if good else []
    evenM, oddM = [totM] * len(pl.P), [0] * len(pl.P)
    for j, i in enumerate(good):
        evenM[i] = char_sum(q, N, LM, mdeg, avoidM + [(d, _jacobi(pl.P[i], M, q))])
        oddM[i] = char_sum(q, N, LPM[j], [d, *mdeg],
                           [(deg(Q), int(pl.chi[Q][i]) * chiM[Q]) for Q in S])
    for i, P in enumerate(pl.P):
        if P in M:
            rest = [Q for Q in M if Q != P]
            if rest:
                oddM[i] = char_sum(q, N, lpoly(q, [_prod(rest, q)])[0], [deg(Q) for Q in rest],
                                   [(deg(Q), _jacobi(Q, rest, q)) for Q in S])
            else:
                oddM[i] = total

    def half(x: list[int], y: list[int]) -> list[int]:
        return [(a + sign * b) // 2 for a, b in zip(x, y, strict=True)]

    return Weights((total + sign * totM) // 2, half(even, evenM), half(odd, oddM))


def _divisors(m: int) -> list[int]:
    return [d for d in range(1, m + 1) if m % d == 0]


def exact_traces(E: Curve, N: int, n: int, sign: int | None = None) -> list[Fraction]:
    q, Einf = E.p, twisted_infinity(E, N)
    W = {d: weights(E, N, d, sign) for d in range(1, n + 1)}
    total = W[1].total
    out = []
    for m in range(1, n + 1):
        s = inf_trace(Einf, m) * total
        for d in _divisors(m):
            pl, w, k = places(E, d), W[d], m // d
            col = w.even if k % 2 == 0 else w.odd
            s += d * sum(local_trace(int(a), int(t), q**d, k) * c
                         for a, t, c in zip(pl.a, pl.kind, col, strict=True))
        out.append(Fraction(-s, q**m * total))
    return out


def limit_traces(E: Curve, N: int, n: int) -> list[Fraction]:
    q, Einf = E.p, twisted_infinity(E, N)
    out = []
    for m in range(1, n + 1):
        s = Fraction(inf_trace(Einf, m))
        if m % 2 == 0:
            for d in _divisors(m // 2):
                pl, qd = places(E, d), q**d
                for a, t in zip(pl.a, pl.kind, strict=True):
                    w = Fraction(qd, qd + 1) if t == GOOD else 1
                    s += d * local_trace(int(a), int(t), qd, m // d) * w
        out.append(-s / q**m)
    return out


def _sym2_local(a: int, kind: int, qd: int, k: int) -> int:
    if kind == ADD:
        raise ValueError("potentially good additive fiber")
    if kind == MULT:
        return 1
    return local_trace(a, GOOD, qd, k) ** 2 - qd**k


def _sym2_inf(E: Curve, m: int) -> int:
    f = E.f_inf
    if f == 0:
        return inf_trace(E, m) ** 2 - E.p**m
    if f == 2 and 3 * deg(E.c4) <= deg(E.disc):
        raise ValueError("potentially good fiber at infinity")
    return 1


def sym2_sums(E: Curve, n: int) -> list[int]:
    q, out = E.p, []
    for m in range(1, n + 1):
        s = _sym2_inf(E, m)
        for d in _divisors(m):
            pl = places(E, d)
            s += d * sum(_sym2_local(int(a), int(t), q**d, m // d)
                         for a, t in zip(pl.a, pl.kind, strict=True))
        out.append(s)
    return out


def sym2_degree(E: Curve) -> int:
    S = primes_dividing(E.disc, E.p)
    if len(multiplicative(E)) < len(S):
        raise ValueError("additive fiber at a finite place")
    _sym2_inf(E, 1)
    bad = sum(deg(Q) for Q in S) + (E.f_inf > 0)
    return 2 * bad - 6


def sym2_lpoly(E: Curve) -> tuple[int, ...]:
    d = sym2_degree(E)
    return tuple(_newton(sym2_sums(E, d))) if d > 0 else (1,)


def ms_traces(E: Curve, N: int, n: int) -> list[Fraction]:
    q, Einf = E.p, twisted_infinity(E, N)
    p2 = power_sums(sym2_lpoly(E), n // 2)
    out = []
    for m in range(1, n + 1):
        if m % 2:
            out.append(Fraction(-inf_trace(Einf, m), q**m))
            continue
        h = m // 2
        D = Fraction(0)
        for d in _divisors(h):
            pl, qd = places(E, d), q**d
            for a, t in zip(pl.a, pl.kind, strict=True):
                a1 = Fraction(local_trace(int(a), int(t), qd, m // d), q**h)
                if t == GOOD:
                    D += d * a1 / (qd + 1)
                else:
                    D -= d * (a1 - Fraction(_sym2_local(int(a), int(t), qd, h // d), q**h) + 1)
        tail = Fraction(_sym2_inf(E, h) - inf_trace(Einf, m), q**m)
        out.append(1 + Fraction(p2[h - 1], q**m) + D / q**h + tail)
    return out


def density_gap(tr: list[Fraction], dim: int, fhat: Callable[[float], float]) -> float:
    return 2 / dim * sum(fhat(m / dim) * (float(t) - (m % 2 == 0))
                         for m, t in enumerate(tr, 1))
