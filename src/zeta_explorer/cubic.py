from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np

from zeta_explorer.gf import GF, Poly, correlate, evaluate, field, minimal_polys, pdivmod, ppowmod


@dataclass(frozen=True)
class Eis:
    a: int
    b: int = 0

    def __add__(self, o: Eis | int) -> Eis:
        o = _eis(o)
        return Eis(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self) -> Eis:
        return Eis(-self.a, -self.b)

    def __sub__(self, o: Eis | int) -> Eis:
        return self + -_eis(o)

    def __mul__(self, o: Eis | int) -> Eis:
        if isinstance(o, int):
            return Eis(self.a * o, self.b * o)
        a, b, c, d = self.a, self.b, o.a, o.b
        return Eis(a * c - b * d, a * d + b * c - b * d)

    __rmul__ = __mul__

    def __pow__(self, k: int) -> Eis:
        out = Eis(1)
        for _ in range(k):
            out = out * self
        return out

    def div(self, k: int) -> Eis:
        if self.a % k or self.b % k:
            raise ArithmeticError(f"{self} is not divisible by {k}")
        return Eis(self.a // k, self.b // k)

    def conj(self) -> Eis:
        return Eis(self.a - self.b, -self.b)

    def rot(self, e: int) -> Eis:
        return self * OMEGA[e % 3]

    @property
    def norm(self) -> int:
        return self.a * self.a - self.a * self.b + self.b * self.b

    @property
    def tr(self) -> int:
        return 2 * self.a - self.b

    def __complex__(self) -> complex:
        return self.a + self.b * complex(-0.5, 3**0.5 / 2)


def _eis(x: Eis | int) -> Eis:
    return x if isinstance(x, Eis) else Eis(x)


OMEGA = (Eis(1), Eis(0, 1), Eis(-1, -1))


def classes(n: np.ndarray | list[int]) -> Eis:
    return Eis(int(n[0]) - int(n[2]), int(n[1]) - int(n[2]))


def zeta3(p: int) -> int:
    if p % 3 != 1:
        raise ValueError("cubic characters need p = 1 mod 3")
    return min(z for z in range(2, p) if pow(z, 3, p) == 1)


@functools.cache
def _index(p: int, m: int) -> np.ndarray:
    F = field(p, m)
    c = int(F.exp[(F.q - 1) // 3])
    s = 1 if c == zeta3(p) else 2
    out = np.full(F.q, -1, dtype=np.int64)
    out[F.exp] = (s * np.arange(F.q - 1)) % 3
    return out


def index(F: GF) -> np.ndarray:
    return _index(F.p, F.m)


def symbol(A: Poly, P: Poly, p: int) -> int | None:
    r = pdivmod(A, P, p)[1]
    if not r:
        return None
    v = ppowmod(r, (p ** (len(P) - 1) - 1) // 3, P, p)
    z = zeta3(p)
    return {(1,): 0, (z,): 1, (z * z % p,): 2}[v]


def _orbit_sums(F: GF, Ds: list[Poly], w: np.ndarray) -> list[Eis]:
    reps = F.orbits[0]
    idx = index(F)
    w = np.asarray(w, dtype=np.int64)
    chunk = max(1, 2**22 // len(reps))
    out = []
    for i in range(0, len(Ds), chunk):
        e = idx[evaluate(F, Ds[i:i + chunk], reps)]
        n = np.stack([(e == j).astype(np.int64) @ w for j in range(3)], axis=1)
        out += [classes(row) for row in n]
    return out


def char_sums(p: int, m: int, Ds: list[Poly], w: np.ndarray | None = None) -> list[Eis]:
    F = field(p, m)
    size = np.bincount(F.orbits[1])
    return _orbit_sums(F, Ds, size if w is None else size * np.asarray(w)[F.orbits[0]])


def newton(S: list[Eis], n: int) -> list[Eis]:
    c = [Eis(1)]
    for k in range(1, n + 1):
        c.append(sum((S[m - 1] * c[k - m] for m in range(1, k + 1)), Eis(0)).div(k))
    return c


def lpolys(p: int, Ps: list[Poly]) -> list[list[Eis]]:
    d = len(Ps[0]) - 1
    S = [char_sums(p, m, Ps) for m in range(1, d)]
    return [newton([S[m][i] for m in range(d - 1)], d - 1) for i in range(len(Ps))]


def logderiv(c: list[Eis], n: int) -> list[Eis]:
    c = c + [Eis(0)] * n
    out: list[Eis] = []
    for m in range(1, n + 1):
        out.append(c[m] * m - sum((c[k] * out[m - k - 1] for k in range(1, m)), Eis(0)))
    return out


@dataclass(frozen=True)
class Jacobi:
    P: list[Poly]
    lam: list[Eis]


def jacobi(p: int, B: Poly, d: int) -> Jacobi:
    F = field(p, d)
    alpha, Ps = minimal_polys(p, d)
    x = np.arange(F.q)
    sq = np.bincount(F.power(x, 2), minlength=F.q)
    idx = index(F)
    corr = [correlate(F, (idx == j).astype(np.int64), sq) for j in range(3)]
    c = F.neg(F.polyval(B, alpha))
    return Jacobi(Ps, [classes([corr[j][ci] for j in range(3)]) for ci in c])
