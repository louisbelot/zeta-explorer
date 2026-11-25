from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from zeta_explorer.funcfield import Curve, Twist, legendre, trace_inf, twisted_legendre
from zeta_explorer.gf import Poly, deg, field, pdivmod, primes_dividing
from zeta_explorer.onelevel import GOOD, MULT, multiplicative, places, sym2_lpoly


def _val(a: Poly, P: Poly, p: int) -> int:
    k = 0
    while a:
        quo, rem = pdivmod(a, P, p)
        if rem:
            break
        a, k = quo, k + 1
    return k


def check_curve(E: Curve) -> None:
    for P in primes_dividing(E.disc, E.p):
        if not E.c4 or 3 * _val(E.c4, P, E.p) >= _val(E.disc, P, E.p):
            raise ValueError("potentially good fiber")
    if E.f_inf == 2 and 3 * deg(E.c4) <= deg(E.disc):
        raise ValueError("potentially good fiber at infinity")


def sym2(E: Curve) -> tuple[int, ...]:
    if E == twisted_legendre(E.p):
        return sym2_lpoly(legendre(E.p))
    return sym2_lpoly(E)


def chi_M(D: Poly, M: list[Poly], p: int) -> int:
    s = 1
    for P in M:
        F = field(p, deg(P))
        root = next(int(a) for a in range(F.q) if F.polyval(P, a) == 0)
        s *= int(F.chi[F.polyval(D, root)])
    return s


def root_number(fam: list[Twist], M: list[Poly]) -> int:
    p = fam[0].L.q
    eps = {t.L.eps * chi_M(t.D, M, p) for t in fam}
    if len(eps) != 1:
        raise ArithmeticError("root number is not eps0 chi_D(M)")
    return eps.pop()


def _g(lam: np.ndarray, kind: np.ndarray, X) -> np.ndarray:
    good = 1 / (1 - lam * X + X * X)
    mult = 1 / (1 - lam * X)
    return np.where(kind == GOOD, good, np.where(kind == MULT, mult, 1))


def _sym2_inv(lam: np.ndarray, kind: np.ndarray, z, Q: int) -> np.ndarray:
    good = (1 - z) * (1 - (lam * lam - 2) * z + z * z)
    return np.where(kind == GOOD, good, 1 - z / Q)


@dataclass(frozen=True)
class Local:
    d: int
    lam: np.ndarray
    kind: np.ndarray
    odd: np.ndarray


def local(E: Curve, d: int) -> Local:
    pl = places(E, d)
    M = set(multiplicative(E))
    odd = np.array([P in M for P in pl.P], dtype=bool)
    lam = pl.a / np.sqrt(float(E.p) ** d)
    return Local(d, lam, pl.kind, odd)


def _factor(L: Local, Q: int, x, y, twisted: bool) -> np.ndarray:
    X, Y = x**L.d / np.sqrt(Q), y**L.d / np.sqrt(Q)
    G = _g(L.lam, L.kind, X) * _g(L.lam, L.kind, Y)
    H = _g(L.lam, L.kind, -X) * _g(L.lam, L.kind, -Y)
    even, odd = (G + H) / 2, (G - H) / 2
    good = 1 + (even - 1) / (1 + 1 / Q)
    return np.where(L.kind == GOOD, good, np.where(twisted & L.odd, odd, even))


def raw2(E: Curve, x, y, d_max: int, twisted: bool = False) -> complex:
    out = 1
    for d in range(1, d_max + 1):
        out *= np.prod(_factor(local(E, d), E.p**d, x, y, twisted))
    return out


def _inf(E: Curve, w) -> complex:
    q = E.p
    if E.f_inf:
        return 1 - w / q**2
    lam = trace_inf(E, 1) / np.sqrt(q)
    z = w / q
    return (1 - z) * (1 - (lam * lam - 2) * z + z * z)


def _poly(c: tuple[int, ...], T) -> complex:
    return sum(ck * T**k for k, ck in enumerate(c))


def h2(E: Curve, x, y, d_max: int, twisted: bool = False) -> complex:
    q, c = E.p, sym2(E)
    out = 1
    for w in (x * x, y * y, x * y):
        out *= _poly(c, w / q**2) * _inf(E, w)
    s = 0
    for d in range(1, d_max + 1):
        L, Q = local(E, d), q**d
        r = _factor(L, Q, x, y, twisted) * (1 - (x * y) ** d / Q)
        for w in (x * x, y * y, x * y):
            r = r * _sym2_inv(L.lam, L.kind, w**d / Q, Q)
        s += np.log(r).sum()
    return out * np.exp(s)
