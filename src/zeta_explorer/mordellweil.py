from __future__ import annotations

from fractions import Fraction

from zeta_explorer.gf import Poly, padd, pmul, pscale, psub
from zeta_explorer.heights import Model, Point, _reduced, psqrt


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
