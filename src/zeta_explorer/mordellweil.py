from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from zeta_explorer.bsd import selmer_rank, torsion
from zeta_explorer.funcfield import Twist
from zeta_explorer.gf import Poly, deg, padd, pmul, pscale, psub
from zeta_explorer.heights import (
    Klass,
    Model,
    Point,
    _reduced,
    delta,
    lstar,
    model,
    points,
    psqrt,
    sqclass,
    torsion_classes,
)

HERMITE = {1: Fraction(1), 2: Fraction(4, 3), 3: Fraction(2)}


def lift(M: Model, P: Point) -> tuple[Poly, Poly, Poly]:
    w = psqrt(P.B, M.p)
    C = psqrt(M.cubic(P.A, P.B), M.p)
    if w is None or C is None:
        raise ArithmeticError("not a point")
    return P.A, w, C


def add(M: Model, P: Point, Q: Point, sign: int = 1) -> Point:
    p = M.p
    A1, w1, C1 = lift(M, P)
    A2, w2, C2 = lift(M, Q)
    w1s, w2s = pmul(w1, w1, p), pmul(w2, w2, p)
    dx = psub(pmul(A2, w1s, p), pmul(A1, w2s, p), p)
    if not dx:
        raise ArithmeticError("equal x")
    dy = psub(pscale(pmul(C2, pmul(w1s, w1, p), p), sign, p), pmul(C1, pmul(w2s, w2, p), p), p)
    dx2 = pmul(dx, dx, p)
    den = pmul(pmul(w1s, w2s, p), dx2, p)
    a2 = pscale(padd(padd(M.e[0], M.e[1], p), M.e[2], p), -1, p)
    num = pmul(dy, dy, p)
    for t in (pmul(a2, den, p), pmul(pmul(A1, w2s, p), dx2, p), pmul(pmul(A2, w1s, p), dx2, p)):
        num = psub(num, t, p)
    return _reduced(num, den, p)


def pairing(M: Model, P: Point, Q: Point) -> Fraction:
    return (M.height(add(M, P, Q)) - M.height(P) - M.height(Q)) / 2


def det(G: list[list[Fraction]]) -> Fraction:
    G = [row[:] for row in G]
    n, out = len(G), Fraction(1)
    for c in range(n):
        r = next((r for r in range(c, n) if G[r][c]), None)
        if r is None:
            return Fraction(0)
        if r != c:
            G[c], G[r] = G[r], G[c]
            out = -out
        out *= G[c][c]
        for r in range(c + 1, n):
            f = G[r][c] / G[c][c]
            G[r] = [a - f * b for a, b in zip(G[r], G[c], strict=True)]
    return out


def gram(M: Model, Ps: list[Point]) -> list[list[Fraction]]:
    n = len(Ps)
    G = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        G[i][i] = M.height(Ps[i])
        for j in range(i):
            G[i][j] = G[j][i] = pairing(M, Ps[i], Ps[j])
    return G


def kmul(a: Klass, b: Klass, p: int) -> Klass:
    return sqclass(pmul(a[0], b[0], p), p), sqclass(pmul(a[1], b[1], p), p)


def cosets(M: Model) -> list[Klass]:
    tors = torsion_classes(M)
    seen, out = set(tors), []
    for c in sorted(M.selmer, key=lambda c: (deg(c[0]) + deg(c[1]), c)):
        if c in seen:
            continue
        out.append(c)
        seen |= {kmul(c, t, M.p) for t in tors}
    return out


def outside(M: Model, bound: int) -> Fraction:
    return 2 * M.chi + 2 * (bound + 1) - M.max_contr


def span(M: Model, cls: list[Klass]) -> set[Klass]:
    S = set(torsion_classes(M))
    for c in cls:
        S |= {kmul(c, s, M.p) for s in S}
    return S


@dataclass(frozen=True)
class Basis:
    points: tuple[Point, ...]
    gram: tuple[tuple[Fraction, ...], ...]
    bound: int
    proof: str

    @property
    def regulator(self) -> Fraction:
        return det([list(r) for r in self.gram])


def minima(M: Model, cands: list[Point], rank: int) -> list[Point]:
    chosen: list[Point] = []
    for P in sorted(cands, key=M.height):
        if det(gram(M, [*chosen, P])) > 0:
            chosen.append(P)
            if len(chosen) == rank:
                break
    return chosen


def saturated(M: Model, chosen: list[Point], cands: list[Point], bound: int) -> str | None:
    r = len(chosen)
    out = outside(M, bound)
    if M.height(chosen[-1]) < out:
        return "minima"
    lam = min(min(map(M.height, cands)), out)
    index2 = HERMITE[r] * det(gram(M, chosen)) / lam**r
    odd = len(span(M, [delta(M, P) for P in chosen])) == 4 * 2**r
    return "index" if index2 < 4 or (odd and index2 < 9) else None


def reduce(M: Model, Ps: list[Point]) -> list[Point]:
    Ps, changed = list(Ps), True
    while changed:
        changed = False
        Ps.sort(key=M.height)
        for i in range(len(Ps)):
            for j in range(i):
                best = min((add(M, Ps[i], Ps[j], s) for s in (1, -1)), key=M.height)
                if M.height(best) < M.height(Ps[i]):
                    Ps[i], changed = best, True
    return Ps


def basis(M: Model, rank: int, max_bound: int = 3) -> Basis | None:
    if not 1 <= rank <= 3:
        raise ValueError("minima need not span past rank 3")
    reps = cosets(M)
    for bound in range(max_bound + 1):
        cands = list({P for c in reps for P in points(M, bound, c)})
        chosen = minima(M, cands, rank)
        if len(chosen) < rank:
            continue
        proof = saturated(M, chosen, cands, bound)
        if proof:
            Ps = reduce(M, chosen)
            return Basis(tuple(Ps), tuple(map(tuple, gram(M, Ps))), bound, proof)
    return None


@dataclass(frozen=True)
class MordellWeil:
    D: Poly
    q: int
    chi: int
    rank: int
    tamagawa: int
    torsion: int
    selmer: int
    lstar: Fraction
    basis: Basis
    sha: Fraction

    @property
    def sha2(self) -> int:
        return self.selmer - self.rank - 2


def solve(t: Twist, max_bound: int = 3) -> MordellWeil | None:
    q, D = t.L.q, t.D
    M = model(D, q)
    B = basis(M, t.rank, max_bound)
    if B is None:
        return None
    tam = math.prod(f.c for f in M.fibers.values())
    tors = torsion(D, q)
    ls = lstar(t.L)
    sha = ls * Fraction(q) ** (M.chi - 1) * tors**2 / (tam * B.regulator)
    return MordellWeil(D, q, M.chi, t.rank, tam, tors, selmer_rank(M.e, q), ls, B, sha)

