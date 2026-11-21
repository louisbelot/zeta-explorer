from __future__ import annotations

import functools
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.cubic import Eis, index, newton
from zeta_explorer.funcfield import Curve, j0_traces, j1728_traces, trace_inf
from zeta_explorer.gf import Poly, deg, evaluate, field, minimal_polys, pmul, primes_dividing


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


def powerfree(p: int, N: int, n: int, avoid: Poly) -> list[Poly]:
    bad = set(primes_dividing(avoid, p))
    primes = [P for d in range(1, N + 1) for P in sorted(minimal_polys(p, d)[1]) if P not in bad]
    out = []

    def grow(f: Poly, start: int):
        if deg(f) == N:
            out.append(f)
            return
        for i in range(start, len(primes)):
            P = primes[i]
            if deg(f) + deg(P) > N:
                break
            g = f
            for _ in range(n - 1):
                g = pmul(g, P, p)
                if deg(g) > N:
                    break
                grow(g, i + 1)

    grow((1,), 0)
    return out


def family_sums(p: int, n: int, B: Poly, Ds: list[Poly], m: int) -> np.ndarray:
    F = field(p, m)
    reps, inv = F.orbits
    size = np.bincount(inv)
    table = (j0_traces if n == 6 else j1728_traces)(p, m)
    cs = [pmul(B, D, p) for D in Ds]
    out = np.empty(len(cs), dtype=np.int64)
    chunk = max(1, 2**22 // len(reps))
    for i in range(0, len(cs), chunk):
        out[i:i + chunk] = table[evaluate(F, cs[i:i + chunk], reps)] @ size
    return out + np.array([trace_inf(curve(n, p, c), m) for c in cs])


def family_traces(p: int, n: int, B: Poly, N: int, top: int) -> list[Fraction]:
    Ds = powerfree(p, N, n, B)
    return [Fraction(-int(family_sums(p, n, B, Ds, m).sum()), p**m * len(Ds))
            for m in range(1, top + 1)]


def _series(N: int, num: list, dens: list[tuple[int, int]]) -> list:
    s = list(num[:N + 1]) + [num[0] * 0] * (N + 1 - len(num))
    for a, e in dens:
        for i in range(e, N + 1):
            s[i] = s[i] + s[i - e] * a
    return s


def _times(f: list, g: list, N: int) -> list:
    out = [f[0] * 0] * (N + 1)
    for i, fi in enumerate(f[:N + 1]):
        for j, gj in enumerate(g[:N + 1 - i]):
            out[i + j] = out[i + j] + fi * gj
    return out


def _one_minus(a, e: int) -> list:
    return [a * 0 + 1] + [a * 0] * (e - 1) + [-a]


@dataclass(frozen=True)
class Places:
    P: list[Poly]
    B: np.ndarray
    R: dict[Poly, np.ndarray]
    counts: list[np.ndarray]


@functools.cache
def places(p: int, n: int, B: Poly, d: int) -> Places:
    alpha, Ps = minimal_polys(p, d)
    F = field(p, d)
    idx = index(F, n)
    sym = {R: idx[F.polyval(R, alpha)] for R in primes_dividing(B, p)}
    counts = [_class_counts(field(p, r), n, Ps) for r in range(1, d)]
    return Places(Ps, idx[F.polyval(B, alpha)], sym, counts)


def char_lpoly(p: int, n: int, pl: Places, i: int, j: int) -> list:
    d = len(pl.P[i]) - 1
    z = units(n)
    S = [classes(n, [sum(int(c[i][e]) for e in range(n) if e * j % n == t) for t in range(n)])
         for c in pl.counts]
    L = newton(S, d - 1) if S else [z[0]]
    eps = (-1) ** (j * ((p - 1) // n) * d % 2)
    return [c * eps**k for k, c in enumerate(L)]


def _size_series(p: int, n: int, B: Poly, N: int) -> list:
    R = primes_dividing(B, p)
    num = [1]
    for r in R:
        num = _times(num, _one_minus(1, deg(r)), N)
    num = _times(num, _one_minus(p, n), N)
    return _series(N, num, [(p, 1)] + [(1, n * deg(r)) for r in R])


def family_size(p: int, n: int, B: Poly, N: int) -> int:
    return _size_series(p, n, B, N)[N]


def weight(p: int, n: int, B: Poly, N: int, pl: Places, i: int, k: int) -> Eis | Gauss:
    d = len(pl.P[i]) - 1
    z = units(n)
    if k % n == 0:
        num = _times(_size_series(p, n, B, N), _one_minus(1, d), N)
        return z[0] * _series(N, num, [(1, n * d)])[N]
    j = -k % n
    num = char_lpoly(p, n, pl, i, j)
    for r, e in pl.R.items():
        num = _times(num, _one_minus(z[e[i] * j % n], deg(r)), N)
    num = _times(num, _one_minus(z[0] * p, n), N)
    dens = [(1, n * deg(r)) for r in pl.R] + [(1, n * d)]
    return _series(N, num, dens)[N]


def _divisors(m: int) -> list[int]:
    return [d for d in range(1, m + 1) if m % d == 0]


def exact_traces(p: int, n: int, B: Poly, N: int, top: int,
                 offdiagonal: bool = True) -> list[Fraction]:
    total = family_size(p, n, B, N)
    z, pi = units(n), frobenius(p, n)
    lead = int(index(field(p, 1), n)[B[-1]])
    out = []
    for m in range(1, top + 1):
        s = z[0] * 0
        if (deg(B) + N) % n == 0:
            s = s + z[-m * lead % n] * total
        for d in _divisors(m):
            k = m // d
            if k % n and not offdiagonal:
                continue
            pl = places(p, n, B, d)
            for i in range(len(pl.P)):
                if pl.B[i] < 0:
                    continue
                s = s + z[-k * int(pl.B[i]) % n] * weight(p, n, B, N, pl, i, k) * d
        out.append(Fraction(-(pi**m * s).tr, p**m * total))
    return out


def limit_traces(p: int, n: int, B: Poly, N: int, top: int) -> list[Fraction]:
    z, pi = units(n), frobenius(p, n)
    lead = int(index(field(p, 1), n)[B[-1]])
    bad = primes_dividing(B, p)
    out = []
    for m in range(1, top + 1):
        s = Fraction(0)
        if (deg(B) + N) % n == 0:
            s += (z[-m * lead % n] * pi**m).tr
        for d in _divisors(m):
            if (m // d) % n == 0:
                good = len(minimal_polys(p, d)[1]) - sum(deg(R) == d for R in bad)
                s += d * good * (1 - Fraction(1, p**d)) / (1 - Fraction(1, p**(n * d))) * (pi**m).tr
        out.append(-s / p**m)
    return out
