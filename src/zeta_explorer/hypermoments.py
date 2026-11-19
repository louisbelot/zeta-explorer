from __future__ import annotations

import itertools
import math
from collections import defaultdict
from collections.abc import Iterator

import mpmath as mp
import numpy as np

from zeta_explorer.funcfield import _newton
from zeta_explorer.gf import (
    Poly,
    character_sums,
    field,
    minimal_polys,
    pdivmod,
    pmul,
    ppowmod,
)
from zeta_explorer.hyperelliptic import irreducible_counts


def irreducibles(q: int, d: int) -> list[Poly]:
    return sorted(minimal_polys(q, d)[1], key=lambda P: P[::-1])


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


def _zeta(q: int, s: float) -> float:
    return 1 / (1 - q ** (1 - s))


def _terms(q: int) -> int:
    return int(60 / math.log10(q))


def _log_D(q: int, z, w, pi: list[int]):
    out = mp.mpf(0)
    for d in range(1, len(pi)):
        P, zd, wd = mp.mpf(q) ** d, z**d, w**d
        num = (-wd**2 - wd**3 / P + wd / (P * P * zd) + P * (w * w * z) ** d
               + (w * w * z) ** d - P * P * (w * z * z) ** d + (w**3 * z) ** d
               - P * P * (w * w * z * z) ** d)
        out += pi[d] * mp.log(1 + num / ((P * P * zd - 1) * (1 + wd)))
    return out


def florea_D(q: int, z, w) -> complex:
    pi = irreducible_counts(q, _terms(q))
    with mp.workdps(30 + len(pi)):
        return complex(_log_D(q, mp.mpmathify(z), mp.mpmathify(w), pi))


def printed_dlog(q: int) -> float:
    pi = irreducible_counts(q, _terms(q))
    out = 0.0
    for d in range(1, len(pi)):
        P = float(q) ** d
        c = P ** (1 / 3)
        out -= pi[d] * d * (P - 1) * (c + 1) / ((c - 1) * (c**4 + P) ** 2)
    return out


def euler_dlog(q: int) -> float:
    pi = irreducible_counts(q, _terms(q))
    with mp.workdps(30 + len(pi)):
        z0 = mp.mpf(q) ** (-mp.mpf(4) / 3)
        return float(z0 * mp.diff(lambda z: _log_D(q, z, q * z, pi), z0))


def florea_R(q: int, x: float, dlog: float | None = None) -> float:
    if q % 4 != 1:
        raise ValueError("Florea's secondary term assumes q = 1 mod 4")
    pi = irreducible_counts(q, _terms(q))
    log_D = 0.0
    for d in range(1, len(pi)):
        P = float(q) ** d
        c = P ** (1 / 3)
        log_D += pi[d] * math.log1p(-(c**4 + c**2 + c + 1) / (c**4 + P) ** 2)
    if dlog is None:
        dlog = printed_dlog(q)
    pre = _zeta(q, 5 / 3) * _zeta(q, 7 / 3) / (9 * q ** (2 / 3) * _zeta(q, 4 / 3) ** 2)
    const = (_zeta(q, 7 / 3) / q ** (4 / 3) * (-2.5 - 2 * q ** (1 / 3) - 2 * q + q ** (4 / 3) / 2)
             - 2 / q ** (4 / 3) * dlog)
    return pre * math.exp(log_D) * (x / 2 + const)


def florea_secondary(q: int, g: int, dlog: float | None = None) -> float:
    n = 2 * g + 1
    return q ** (n / 3) * florea_R(q, n, dlog) / (q**n - q ** (n - 1))


def _depth(q: int) -> int:
    return int(25 / math.log10(q))


def _arithmetic(q: int, beta: list, pi: list[int]):
    out = mp.mpf(0)
    k = len(beta)
    for d in range(1, len(pi)):
        P = mp.mpf(q) ** d
        x = [P ** (-mp.mpf(1) / 2 - b) for b in beta]
        plus = minus = mp.mpf(1)
        for xi in x:
            plus /= 1 - xi
            minus /= 1 + xi
        loc = 1 + P / (P + 1) * ((plus + minus) / 2 - 1)
        for i in range(k):
            for j in range(i, k):
                loc *= 1 - P ** (-1 - beta[i] - beta[j])
        out += pi[d] * mp.log(loc)
    return out


def _A(q: int, beta: list, pi: list[int]):
    L = mp.log(q)
    zeta = mp.mpf(1)
    for i in range(len(beta)):
        for j in range(i, len(beta)):
            zeta /= 1 - mp.exp(-(beta[i] + beta[j]) * L)
    return zeta * mp.exp(_arithmetic(q, beta, pi))


def arithmetic_factor(q: int, k: int) -> float:
    with mp.workdps(30):
        return float(mp.exp(_arithmetic(q, [mp.mpf(0)] * k, irreducible_counts(q, _depth(q)))))


def recipe(q: int, gs: list[int], k: int, points: int = 32) -> list[float]:
    with mp.workdps(60):
        L = mp.log(q)
        pi = irreducible_counts(q, _depth(q))
        v = [mp.mpf(i + 1) for i in range(k)]
        r = mp.mpf(1) / (40 * L * k * max(max(gs), 1))
        acc = [mp.mpc(0)] * len(gs)
        for m in range(points):
            alpha = [r * mp.expjpi(2 * mp.mpf(m) / points) * vi for vi in v]
            for eps in itertools.product((1, -1), repeat=k):
                A = _A(q, [e * a for e, a in zip(eps, alpha, strict=True)], pi)
                flip = sum(a for e, a in zip(eps, alpha, strict=True) if e < 0)
                for i, g in enumerate(gs):
                    acc[i] += A * mp.exp(-2 * g * L * flip)
        return [float(mp.re(a) / points) for a in acc]
