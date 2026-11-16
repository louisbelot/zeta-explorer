from __future__ import annotations

from dataclasses import dataclass

from zeta_explorer.gf import (
    Poly,
    deg,
    norm,
    pdivmod,
    pmul,
    ppowmod,
    primes_dividing,
    psub,
)

Place = Poly | None
Roots = tuple[Poly, Poly, Poly]


def local(g: Poly, P: Place, p: int) -> tuple[int, int]:
    if not g:
        raise ZeroDivisionError("class of zero")
    if P is None:
        return -deg(g), 1 if pow(g[-1], (p - 1) // 2, p) == 1 else -1
    v = 0
    while True:
        quo, rem = pdivmod(g, P, p)
        if rem:
            break
        g, v = quo, v + 1
    return v, 1 if ppowmod(rem, (p ** deg(P) - 1) // 2, P, p) == (1,) else -1


def klass(g: Poly, P: Place, p: int) -> tuple[int, int]:
    v, chi = local(g, P, p)
    return v % 2, (1 - chi) // 2


@dataclass(frozen=True)
class Fiber:
    kind: str
    c: int
    disc: int
    cond: int


def kodaira(e: Roots, P: Place, p: int) -> Fiber:
    pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
    vals = [local(psub(e[i], e[j], p), P, p)[0] for i, j, _ in pairs]
    m, M = min(vals), max(vals)
    m0, M0 = m % 2, M - (m - m % 2)
    if m0 == 1:
        n = 2 * (M0 - 1)
        return Fiber(f"I{n}*", 4, n + 6, 2)
    if M0 == 0:
        return Fiber("good", 1, 0, 0)
    i, _, k = pairs[vals.index(M)]
    split = local(psub(e[i], e[k], p), P, p)[1] == 1
    return Fiber(f"I{2 * M0}", 2 * M0 if split else 2, 2 * M0, 1)


def legendre_roots(D: Poly, p: int) -> Roots:
    return (), norm(D, p), pmul((0, 1), D, p)


def bad_places(e: Roots, p: int) -> list[Place]:
    out = {P for i, j in ((0, 1), (0, 2), (1, 2)) for P in primes_dividing(psub(e[i], e[j], p), p)}
    return [*sorted(out, key=lambda P: (deg(P), P[::-1])), None]


def place_degree(P: Place) -> int:
    return 1 if P is None else deg(P)


def fibers(e: Roots, p: int) -> dict[Place, Fiber]:
    return {P: kodaira(e, P, p) for P in bad_places(e, p)}
