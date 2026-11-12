from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.gf import (
    Poly,
    correlate,
    cubic_sums,
    deg,
    distinct_roots,
    field,
    is_prime,
    monics,
    norm,
    padd,
    pgcd,
    pmul,
    pscale,
    psub,
    squarefree,
)
from zeta_explorer.ranktwo import ks2
from zeta_explorer.rmt import char_poly_at_one, o_minus_angles, so_angles


def _cubic_disc(a2: Poly, a4: Poly, a6: Poly, p: int) -> Poly:
    terms = [(1, (a2, a2, a4, a4)), (-4, (a4, a4, a4)), (-4, (a2, a2, a2, a6)),
             (18, (a2, a4, a6)), (-27, (a6, a6))]
    out: Poly = ()
    for c, f in terms:
        out = padd(out, pscale(functools.reduce(lambda x, y: pmul(x, y, p), f), c, p), p)
    return out


@dataclass(frozen=True)
class Curve:
    p: int
    a2: Poly
    a4: Poly
    a6: Poly

    def __post_init__(self):
        for name in ("a2", "a4", "a6"):
            object.__setattr__(self, name, norm(getattr(self, name), self.p))
        if self.p == 2 or not is_prime(self.p):
            raise ValueError("p must be an odd prime")
        if not self.disc:
            raise ValueError("singular generic fiber")

    @property
    def disc(self) -> Poly:
        return _cubic_disc(self.a2, self.a4, self.a6, self.p)

    @property
    def c4(self) -> Poly:
        return psub(pmul(self.a2, self.a2, self.p), pscale(self.a4, 3, self.p), self.p)

    def twist(self, D: Poly) -> Curve:
        p = self.p
        D2 = pmul(D, D, p)
        return Curve(p, pmul(D, self.a2, p), pmul(D2, self.a4, p), pmul(pmul(D2, D, p), self.a6, p))

    @property
    def infinity(self) -> tuple[int, int, int]:
        coeffs = (self.a2, self.a4, self.a6)
        k = max(math.ceil(deg(a) / i) for a, i in zip(coeffs, (2, 4, 6), strict=True) if a)
        return tuple(a[i * k] if deg(a) == i * k else 0
                     for a, i in zip(coeffs, (2, 4, 6), strict=True))

    @property
    def f_inf(self) -> int:
        a2, a4, a6 = ((c,) for c in self.infinity)
        if _cubic_disc(a2, a4, a6, self.p):
            return 0
        return 1 if psub(pmul(a2, a2, self.p), pscale(a4, 3, self.p), self.p) else 2

    @functools.cached_property
    def f_finite(self) -> int:
        return distinct_roots(self.disc, self.p) + distinct_roots(
            pgcd(self.disc, self.c4, self.p), self.p)

    @property
    def degree(self) -> int:
        n = self.f_finite + self.f_inf - 4
        if n < 0:
            raise ValueError("constant curve: L is not a polynomial of degree deg N - 4")
        return n


def legendre(p: int) -> Curve:
    return Curve(p, (-1, -1), (0, 1), ())


def ulmer(p: int, d: int) -> Curve:
    return Curve(p, (pow(4, -1, p),), (), (0,) * d + (-1,))


@functools.cache
def traces(E: Curve, m: int) -> np.ndarray:
    F = field(E.p, m)
    if E == legendre(E.p):
        x = np.arange(F.q)
        return -correlate(F, F.chi[F.mul(x, F.add(x, F.neg(np.ones_like(x))))], F.chi)
    reps, inv = F.orbits
    A, B, C = (F.polyval(a, reps) for a in (E.a2, E.a4, E.a6))
    return -cubic_sums(F, A, B, C)[inv]


@functools.cache
def _cubic_trace(p: int, m: int, cubic: tuple[int, int, int]) -> int:
    return -int(cubic_sums(field(p, m), *cubic)[0])


def trace_inf(E: Curve, m: int) -> int:
    return _cubic_trace(E.p, m, E.infinity)


def power_sums(E: Curve, M: int) -> list[int]:
    return [int(traces(E, m).sum()) + trace_inf(E, m) for m in range(1, M + 1)]


def _newton(S: list[int]) -> list[int]:
    c = [1]
    for k in range(1, len(S) + 1):
        num = sum(S[m - 1] * c[k - m] for m in range(1, k + 1))
        if num % k:
            raise ArithmeticError("power sums do not come from an integer polynomial")
        c.append(num // k)
    return c


def _sign(c: list[int], N: int, q: int) -> int | None:
    eps = None
    for k in range(N // 2 + 1, min(N, len(c) - 1) + 1):
        a, b = c[k], q ** (2 * k - N) * c[N - k]
        if b == 0:
            if a:
                raise ArithmeticError("functional equation fails")
            continue
        if abs(a) != abs(b):
            raise ArithmeticError("functional equation fails")
        e = 1 if a == b else -1
        if eps is not None and e != eps:
            raise ArithmeticError("inconsistent root number")
        eps = e
    if N % 2 == 0 and len(c) > N // 2 and c[N // 2]:
        if eps == -1:
            raise ArithmeticError("odd functional equation with a central coefficient")
        eps = 1
    return eps


def _divide(c: list[int], a: int) -> tuple[list[int], int]:
    out, acc = [], 0
    for ck in c:
        acc = ck + a * acc
        out.append(acc)
    return out[:-1], out[-1]


@dataclass(frozen=True)
class LPoly:
    q: int
    c: tuple[int, ...]
    eps: int

    @property
    def N(self) -> int:
        return len(self.c) - 1

    def value(self, T: Fraction) -> Fraction:
        return sum((Fraction(ck) * T**k for k, ck in enumerate(self.c)), Fraction(0))

    @property
    def central(self) -> Fraction:
        return self.value(Fraction(1, self.q))

    def order(self, sign: int = 1) -> int:
        c, r = list(self.c), 0
        while len(c) > 1:
            quo, rem = _divide(c, sign * self.q)
            if rem:
                break
            c, r = quo, r + 1
        return r

    def deflated(self) -> list[int]:
        c = list(self.c)
        for sign in (1, -1):
            for _ in range(self.order(sign)):
                c = _divide(c, sign * self.q)[0]
        return c

    def roots(self) -> np.ndarray:
        c = self.deflated()
        if len(c) == 1:
            return np.zeros(0, dtype=complex)
        u = np.array([ck / self.q**k for k, ck in enumerate(c)])
        return 1 / np.roots(u[::-1])

    def pairs(self) -> np.ndarray:
        return np.sort(np.abs(np.angle(self.roots())))[::2]

    def so_angles(self) -> np.ndarray:
        if self.N % 2:
            raise ValueError("odd degree")
        r, s = self.order(1), self.order(-1)
        if self.eps == -1:
            r, s = r - 1, s - 1
        return np.sort(np.concatenate([np.zeros(r // 2), np.full(s // 2, np.pi), self.pairs()]))


def from_power_sums(S: list[int], N: int, q: int) -> LPoly | None:
    c = _newton(S[:N])
    eps = 1 if N == 0 else _sign(c, N, q)
    if eps is None:
        return None
    c += [0] * (N + 1 - len(c))
    for k in range(N // 2 + 1, N + 1):
        c[k] = eps * q ** (2 * k - N) * c[N - k]
    return LPoly(q, tuple(c), eps)


def lpoly(E: Curve) -> LPoly:
    N = E.degree
    M = min(N, N // 2 + 1)
    while True:
        L = from_power_sums(power_sums(E, M), N, E.p)
        if L is not None:
            return L
        M += 1


@dataclass(frozen=True)
class Twist:
    D: Poly
    L: LPoly

    @property
    def rank(self) -> int:
        return self.L.order(1)

    @property
    def w(self) -> int:
        if self.L.N % 2:
            raise ValueError("odd degree")
        v = self.L.central * self.L.q ** (self.L.N // 2)
        if v.denominator != 1:
            raise ArithmeticError("central value is not in q^(-N/2) Z")
        return int(v)


def twist_set(p: int, n: int, avoid: Poly) -> list[Poly]:
    return [D for D in monics(p, n) if squarefree(D, p) and deg(pgcd(D, avoid, p)) == 0]


def _twisted_sums(E: Curve, Ds: list[Poly], m: int) -> np.ndarray:
    F = field(E.p, m)
    reps, inv = F.orbits
    a = (traces(E, m)[reps] * np.bincount(inv)).astype(np.float64)
    P = np.stack([F.power(reps, j) for j in range(len(Ds[0]))])
    digits = [F.dig[P, j].astype(np.float64) for j in range(F.m)]
    C = np.array(Ds, dtype=np.float64)
    chunk = max(1, 2**22 // len(reps))
    out = np.empty(len(Ds), dtype=np.int64)
    for i in range(0, len(Ds), chunk):
        rows = C[i:i + chunk]
        code = np.zeros((len(rows), len(reps)), dtype=np.int64)
        for d, w in zip(digits, F.pw, strict=True):
            code += np.fmod(rows @ d, E.p).astype(np.int64) * w
        out[i:i + chunk] = np.rint(F.chi[code] @ a).astype(np.int64)
    return out + [trace_inf(E.twist(D), m) for D in Ds]


def family(E: Curve, n: int) -> list[Twist]:
    p = E.p
    Ds = twist_set(p, n, E.disc)
    N = E.f_finite + 2 * n + E.twist(Ds[0]).f_inf - 4
    S: list[list[int]] = [[] for _ in Ds]
    out: list[LPoly | None] = [None] * len(Ds)
    todo, m = list(range(len(Ds))), 0
    while todo:
        m += 1
        for i, s in zip(todo, _twisted_sums(E, [Ds[i] for i in todo], m), strict=True):
            S[i].append(int(s))
        if m > N // 2:
            for i in todo:
                out[i] = from_power_sums(S[i], N, p)
            todo = [i for i in todo if out[i] is None]
    return [Twist(D, L) for D, L in zip(Ds, out, strict=True)]


def lowest(fam: list[Twist], eps: int = 1) -> np.ndarray:
    rank = 0 if eps == 1 else 1
    return np.array([t.L.so_angles()[0] for t in fam if t.L.eps == eps and t.rank == rank])


def central_k(fam: list[Twist]) -> np.ndarray:
    q, e = fam[0].L.q, fam[0].L.N // 4
    k2 = [Fraction(t.w, 4 * q**e) for t in fam if t.L.eps == 1]
    k = [math.isqrt(int(v)) for v in k2]
    if any(v.denominator != 1 or r * r != v for v, r in zip(k2, k, strict=True)):
        raise ArithmeticError("a central value is not 4 q^(N/4) k^2")
    return np.array(k)


@dataclass(frozen=True)
class Fit:
    q: int
    N: int
    even: int
    odd: int
    rank2: float
    predicted: float
    ks_so: float
    ks_excised: float
    ks_odd: float


def compare(fam: list[Twist], cut: float, size: int = 200_000, seed: int = 1) -> Fit:
    q, N = fam[0].L.q, fam[0].L.N
    if N % 2:
        raise ValueError("odd degree")
    even, odd = lowest(fam, 1), lowest(fam, -1)
    H = so_angles(N // 2, size, seed)
    keep = char_poly_at_one(H) >= cut
    G = o_minus_angles(N // 2, size, seed + 1)
    return Fit(q, N, len(even), len(odd),
               float(np.mean([t.rank > 0 for t in fam if t.L.eps == 1])), float(1 - keep.mean()),
               ks2(even, H[:, 0]), ks2(even, H[keep, 0]), ks2(odd, G[:, 0]))
