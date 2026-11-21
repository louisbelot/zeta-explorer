from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np

from zeta_explorer.cubic import Eis, index, newton
from zeta_explorer.funcfield import Curve
from zeta_explorer.gf import Poly, deg, evaluate, field


@dataclass(frozen=True)
class Gauss:
    a: int
    b: int = 0

    def __add__(self, o: Gauss | int) -> Gauss:
        o = _gauss(o)
        return Gauss(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self) -> Gauss:
        return Gauss(-self.a, -self.b)

    def __sub__(self, o: Gauss | int) -> Gauss:
        return self + -_gauss(o)

    def __mul__(self, o: Gauss | int) -> Gauss:
        if isinstance(o, int):
            return Gauss(self.a * o, self.b * o)
        return Gauss(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __pow__(self, k: int) -> Gauss:
        out = Gauss(1)
        for _ in range(k):
            out = out * self
        return out

    def div(self, k: int) -> Gauss:
        if self.a % k or self.b % k:
            raise ArithmeticError(f"{self} is not divisible by {k}")
        return Gauss(self.a // k, self.b // k)

    def conj(self) -> Gauss:
        return Gauss(self.a, -self.b)

    @property
    def norm(self) -> int:
        return self.a * self.a + self.b * self.b

    @property
    def tr(self) -> int:
        return 2 * self.a

    def __complex__(self) -> complex:
        return complex(self.a, self.b)


def _gauss(x: Gauss | int) -> Gauss:
    return x if isinstance(x, Gauss) else Gauss(x)


ZETA = {4: Gauss(0, 1), 6: Eis(1, 1)}


@functools.cache
def units(n: int) -> tuple:
    return tuple(ZETA[n] ** e for e in range(n))


def classes(n: int, counts) -> Eis | Gauss:
    z = units(n)
    return sum((z[e] * int(c) for e, c in enumerate(counts)), z[0] * 0)


def curve(n: int, p: int, c: Poly) -> Curve:
    if n == 6:
        return Curve(p, (), (), c)
    if n == 4:
        return Curve(p, (), c, ())
    raise ValueError("n must be 4 (y^2 = x^3 + cx) or 6 (y^2 = x^3 + c)")


def jacobi(p: int, n: int, m: int = 1) -> Eis | Gauss:
    F = field(p, m)
    s = np.arange(1, F.q)
    e = index(F, n)[s] * (2 if n == 6 else 1) % n
    w = F.chi[F.add(np.ones_like(s), F.neg(s))].astype(np.int64)
    return classes(n, np.bincount(e, weights=w, minlength=n).astype(np.int64))


def frobenius(p: int, n: int, m: int = 1) -> Eis | Gauss:
    J = jacobi(p, n, m)
    if n == 4 and (p**m - 1) // 4 % 2:
        return J
    return -J


def _class_counts(F, n: int, cs: list[Poly]) -> np.ndarray:
    reps, inv = F.orbits
    size = np.bincount(inv)
    idx = index(F, n)
    chunk = max(1, 2**22 // len(reps))
    out = np.zeros((len(cs), n), dtype=np.int64)
    for i in range(0, len(cs), chunk):
        e = idx[evaluate(F, cs[i:i + chunk], reps)]
        for j in range(n):
            out[i:i + chunk, j] = (e == j).astype(np.int64) @ size
    return out


def hecke_sums(p: int, n: int, cs: list[Poly], m: int) -> list[Eis | Gauss]:
    counts = _class_counts(field(p, m), n, cs)
    lead = index(field(p, 1), n)
    z = units(n)
    out = []
    for c, row in zip(cs, counts, strict=True):
        s = classes(n, np.roll(row[::-1], 1))
        if deg(c) % n == 0:
            s = s + z[-m * int(lead[c[-1]]) % n]
        out.append(s)
    return out


def hecke(p: int, n: int, cs: list[Poly]) -> list[list[Eis | Gauss]]:
    h = max(curve(n, p, c).degree for c in cs) // 2
    S = [hecke_sums(p, n, cs, m) for m in range(1, h + 1)]
    out = []
    for i in range(len(cs)):
        H = newton([S[m][i] for m in range(h)], h)
        while len(H) > 1 and H[-1] == H[0] * 0:
            H.pop()
        out.append(H)
    return out


def split(H: list[Eis | Gauss], pi: Eis | Gauss) -> list[int]:
    a = [Hk * pi**k for k, Hk in enumerate(H)]
    h = len(H) - 1
    out = []
    for k in range(2 * h + 1):
        s = sum((a[i] * a[k - i].conj() for i in range(max(0, k - h), min(k, h) + 1)), a[0] * 0)
        if s.b:
            raise ArithmeticError("H(πu)H̄(π̄u) is not real")
        out.append(s.a)
    return out


def angles(H: list[Eis | Gauss], pi: Eis | Gauss, q: int) -> np.ndarray:
    if len(H) == 1:
        return np.zeros(0)
    c = np.array([complex(Hk * pi**k) / q**k for k, Hk in enumerate(H)])
    return np.angle(1 / np.roots(c[::-1]))
