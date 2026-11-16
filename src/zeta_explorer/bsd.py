from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.funcfield import Twist, legendre, traces
from zeta_explorer.gf import (
    Poly,
    deg,
    field,
    norm,
    pdivmod,
    pmul,
    ppowmod,
    primes_dividing,
    psub,
)
from zeta_explorer.ranktwo import ks2
from zeta_explorer.rmt import char_poly_at_one, so_angles

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


def _rank2(rows: list[list[int]], ncols: int) -> int:
    rows = [int("".join(map(str, r)), 2) for r in rows if any(r)]
    rank = 0
    for bit in range(ncols - 1, -1, -1):
        pivot = next((r for r in rows if r >> bit & 1), None)
        if pivot is None:
            continue
        rows = [r ^ pivot if r >> bit & 1 else r for r in rows if r != pivot]
        rank += 1
    return rank


def _span(vs: list[tuple[int, ...]]) -> int:
    return _rank2([list(v) for v in vs], len(vs[0])) if vs else 0


def _delta(x: tuple[Poly, Poly], e: Roots, P: Place, p: int) -> tuple[int, ...] | None:
    a, b = x
    parts = [psub(a, pmul(b, ei, p), p) for ei in e]
    if not all(parts):
        return None
    cb = klass(b, P, p)
    cls = [klass(t, P, p) for t in parts]
    if any((sum(c[j] for c in cls) + cb[j]) % 2 for j in range(2)):
        return None
    return tuple((cls[i][j] + cb[j]) % 2 for i in range(2) for j in range(2))


def _torsion_images(e: Roots, P: Place, p: int) -> list[tuple[int, ...]]:
    d = {(i, j): psub(e[i], e[j], p) for i in range(3) for j in range(3) if i != j}
    first = (pmul(d[0, 1], d[0, 2], p), d[0, 1])
    second = (d[1, 0], pmul(d[1, 0], d[1, 2], p))
    return [tuple(c for g in img for c in klass(g, P, p)) for img in (first, second)]


def local_image(e: Roots, P: Place, p: int, tries: int = 20_000,
                seed: int = 0) -> list[tuple[int, ...]]:
    found = [v for v in _torsion_images(e, P, p) if any(v)]
    if _span(found) == 2:
        return _basis(found)
    rng = np.random.default_rng(seed)
    width = 3 * (1 if P is None else deg(P)) + 2
    for _ in range(tries):
        a = norm([int(c) for c in rng.integers(0, p, rng.integers(1, width + 4))], p)
        if not a:
            continue
        j = int(rng.integers(0, 4))
        b = (1,) if P is None else functools.reduce(lambda x, y: pmul(x, y, p), [P] * j, (1,))
        if P is None and j:
            b = norm([int(c) for c in rng.integers(0, p, 2 * j)] + [1], p)
        v = _delta((a, b), e, P, p)
        if v is not None and any(v) and _span(found + [v]) > _span(found):
            found.append(v)
            if _span(found) == 2:
                return _basis(found)
    raise RuntimeError("local image not found")


def _basis(vs: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for v in vs:
        if _span(out + [v]) > len(out):
            out.append(v)
    return out


def _annihilator(W: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    cand = [tuple((k >> i) & 1 for i in range(4)) for k in range(1, 16)]
    out: list[tuple[int, ...]] = []
    for a in cand:
        if all(sum(x * y for x, y in zip(a, w, strict=True)) % 2 == 0 for w in W):
            if _span(out + [a]) > len(out):
                out.append(a)
    return out


def nonsquare(p: int) -> int:
    return next(g for g in range(2, p) if pow(g, (p - 1) // 2, p) == p - 1)


def selmer_rank(e: Roots, p: int) -> int:
    S = bad_places(e, p)
    gens = [(nonsquare(p),), *[P for P in S if P is not None]]
    n = len(gens)
    rows: list[list[int]] = []
    for P in S:
        A = _annihilator(local_image(e, P, p))
        res = [klass(g, P, p) for g in gens]
        for a in A:
            row = [0] * (2 * n)
            for comp in range(2):
                for k, r in enumerate(res):
                    row[comp * n + k] = (a[2 * comp] * r[0] + a[2 * comp + 1] * r[1]) % 2
            rows.append(row)
    return 2 * n - _rank2(rows, 2 * n)


def point_counts(D: Poly, p: int, m: int) -> np.ndarray:
    F = field(p, m)
    t = F.orbits[0]
    s, size = t, np.zeros(len(t), dtype=np.int64)
    for k in range(1, m + 1):
        s = F.frobenius(s)
        size[(size == 0) & (s == t)] = k
    x = F.polyval(D, t)
    good = (size == m) & (x != 0) & (t != 0) & (F.polyval((p - 1, 1), t) != 0)
    a = traces(legendre(p), m)[t[good]]
    return F.q + 1 - F.chi[x[good]].astype(np.int64) * a


def torsion(D: Poly, p: int, max_degree: int = 6) -> int:
    g = 0
    for m in range(1, max_degree + 1):
        g = math.gcd(g, *map(int, point_counts(D, p, m)))
        odd = g >> ((g & -g).bit_length() - 1) if g else 0
        if odd == 1:
            return 4
    raise ArithmeticError(f"odd torsion not ruled out, gcd {g}")


@dataclass(frozen=True)
class Arith:
    D: Poly
    q: int
    N: int
    rank: int
    fibers: dict[Place, Fiber]
    tamagawa: int
    height: int
    selmer: int
    torsion: int
    sha: Fraction | None

    @property
    def omega(self) -> int:
        return sum(1 for P in self.fibers if P is not None and self.fibers[P].kind == "I0*")

    @property
    def sha2(self) -> int:
        return self.selmer - self.rank - 2

    @property
    def tamagawa2(self) -> int:
        if self.N % 4:
            raise ValueError("k is defined for N divisible by 4")
        v = _v2(Fraction(self.tamagawa, 4 * self.torsion**2))
        if v % 2:
            raise ArithmeticError("odd 2-adic valuation in the Tamagawa factor")
        return v // 2

    @property
    def floor2(self) -> int:
        if self.sha2 % 2:
            raise ArithmeticError("Sha[2] is not a square")
        return self.tamagawa2 + self.sha2 // 2


def _v2(x: Fraction) -> int:
    n, d, v = x.numerator, x.denominator, 0
    while n % 2 == 0:
        n, v = n // 2, v + 1
    while d % 2 == 0:
        d, v = d // 2, v - 1
    return v


def arith(t: Twist) -> Arith:
    q, D = t.L.q, t.D
    e = legendre_roots(D, q)
    fib = fibers(e, q)
    tam = math.prod(f.c for f in fib.values())
    disc = sum(f.disc * place_degree(P) for P, f in fib.items())
    if disc % 12:
        raise ArithmeticError("degree of the minimal discriminant is not a multiple of 12")
    tors = torsion(D, q)
    sha = None
    if t.rank == 0:
        sha = t.L.central * tors**2 * Fraction(q) ** (disc // 12 - 1) / tam
    return Arith(D, q, t.L.N, t.rank, fib, tam, disc // 12, selmer_rank(e, q), tors, sha)


def v2k(t: Twist) -> int:
    k = math.isqrt(t.w // (4 * t.L.q ** (t.L.N // 4)))
    return (k & -k).bit_length() - 1


def ks_mix(x: np.ndarray, first: np.ndarray, lam: np.ndarray, cuts: np.ndarray) -> float:
    grid = np.sort(np.concatenate([x, first]))
    model = np.zeros(len(grid))
    vals, counts = np.unique(cuts, return_counts=True)
    for c, w in zip(vals, counts, strict=True):
        s = np.sort(first[lam >= c])
        model += w * np.searchsorted(s, grid, side="right") / len(s)
    fx = np.searchsorted(np.sort(x), grid, side="right") / len(x)
    return float(np.max(np.abs(fx - model / len(cuts))))


@dataclass(frozen=True)
class Group:
    a: int
    twists: int
    vanish: float
    predicted: float
    lowest: np.ndarray
    model: np.ndarray


@dataclass(frozen=True)
class Excision:
    q: int
    N: int
    rank0: int
    exact: float
    vanish: float
    predicted_one: float
    predicted_local: float
    ks_so: float
    ks_one: float
    ks_tamagawa: float
    ks_selmer: float
    groups: tuple[Group, ...]


def excise(fam: list[Twist], scale: float = 0.25, size: int = 200_000,
           seed: int = 1) -> Excision:
    q, N = fam[0].L.q, fam[0].L.N
    if N % 4:
        raise ValueError("N must be divisible by 4")
    even = [t for t in fam if t.L.eps == 1]
    info = [arith(t) for t in even]
    a = np.array([A.tamagawa2 for A in info])
    zero = np.array([t.rank == 0 for t in even])
    x = np.array([t.L.so_angles()[0] for t in even if t.rank == 0])
    one = 4 * float(q) ** (-N / 4)
    base = scale * one
    H = so_angles(N // 2, size, seed)
    lam, first = char_poly_at_one(H), H[:, 0]
    tam = base * 4.0**a
    sel = np.array([base * 4.0**A.floor2 for A in info if A.rank == 0])
    exact = np.mean([v2k(t) == A.floor2 for t, A in zip(even, info, strict=True) if t.rank == 0])
    groups = tuple(
        Group(int(g), int(np.sum(a == g)), float(1 - zero[a == g].mean()),
              float(np.mean(lam < base * 4.0**g)),
              x[a[zero] == g], first[lam >= base * 4.0**g])
        for g in np.unique(a))
    return Excision(
        q, N, len(x), float(exact), float(1 - zero.mean()),
        float(np.mean(lam < one)), float(np.mean([np.mean(lam < c) for c in tam])),
        ks2(x, first), ks_mix(x, first, lam, np.full(len(x), one)),
        ks_mix(x, first, lam, tam[zero]), ks_mix(x, first, lam, sel), groups)
