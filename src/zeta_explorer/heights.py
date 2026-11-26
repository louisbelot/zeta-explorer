from __future__ import annotations

import functools
import itertools
import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.bsd import (
    Fiber,
    Place,
    Roots,
    _annihilator,
    fibers,
    klass,
    legendre_roots,
    local,
    local_image,
    nonsquare,
    place_degree,
    selmer_rank,
    torsion,
)
from zeta_explorer.funcfield import LPoly, Twist, _divide
from zeta_explorer.gf import (
    GF,
    Poly,
    deg,
    field,
    monics,
    norm,
    padd,
    pdivmod,
    pgcd,
    pmul,
    primes_dividing,
    pscale,
    psub,
)

MAX_CONTR = {"good": Fraction(0), "I2": Fraction(1, 2), "I0*": Fraction(1),
             "I2*": Fraction(3, 2)}


def psqrt(f: Poly, p: int) -> Poly | None:
    if not f:
        return ()
    if deg(f) % 2:
        return None
    n = deg(f) // 2
    r = next((s for s in range(1, p) if s * s % p == f[-1]), None)
    if r is None:
        return None
    s = [0] * (n + 1)
    s[n] = r
    inv = pow(2 * r, -1, p)
    for k in range(1, n + 1):
        acc = sum(s[i] * s[2 * n - k - i] for i in range(n - k + 1, n))
        s[n - k] = (f[2 * n - k] - acc) * inv % p
    s = norm(s, p)
    return s if pmul(s, s, p) == f else None


def _prod(fs, p: int) -> Poly:
    return functools.reduce(lambda a, b: pmul(a, b, p), fs, (1,))


def _val(g: Poly, P: Place, p: int) -> float:
    return math.inf if not g else local(g, P, p)[0]


@dataclass(frozen=True)
class Point:
    A: Poly
    B: Poly


@dataclass(frozen=True)
class Model:
    p: int
    e: Roots
    chi: int
    fibers: dict[Place, Fiber]

    def cubic(self, A: Poly, B: Poly) -> Poly:
        return _prod([psub(A, pmul(ei, B, self.p), self.p) for ei in self.e], self.p)

    def dcubic(self, A: Poly, B: Poly) -> Poly:
        p, out = self.p, ()
        for i in range(3):
            out = padd(out, _prod([psub(A, pmul(self.e[j], B, p), p)
                                   for j in range(3) if j != i], p), p)
        return out

    @functools.cached_property
    def disc(self) -> Poly:
        p, e = self.p, self.e
        d = _prod([psub(e[i], e[j], p) for i, j in ((0, 1), (0, 2), (1, 2))], p)
        return pscale(pmul(d, d, p), 16, p)

    @functools.cached_property
    def selmer(self) -> list[tuple[Poly, Poly]]:
        return selmer_classes(self)

    @property
    def max_contr(self) -> Fraction:
        return sum((place_degree(P) * MAX_CONTR[f.kind] for P, f in self.fibers.items()),
                   Fraction(0))

    def double(self, P: Point) -> Point:
        p, (A, B), e = self.p, (P.A, P.B), self.e
        a2 = pscale(padd(padd(e[0], e[1], p), e[2], p), -1, p)
        a4 = padd(padd(pmul(e[0], e[1], p), pmul(e[0], e[2], p), p), pmul(e[1], e[2], p), p)
        a6 = pscale(_prod(e, p), -1, p)
        A2, B2 = pmul(A, A, p), pmul(B, B, p)
        terms = [pmul(A2, A2, p), pscale(_prod([a4, A2, B2], p), -2, p),
                 pscale(_prod([a6, A, B2, B], p), -8, p),
                 pmul(psub(pmul(a4, a4, p), pscale(pmul(a2, a6, p), 4, p), p), pmul(B2, B2, p), p)]
        num = functools.reduce(lambda x, y: padd(x, y, p), terms, ())
        den = pscale(pmul(B, self.cubic(A, B), p), 4, p)
        return _reduced(num, den, p)

    def translate(self, P: Point, i: int) -> Point:
        p, e = self.p, self.e
        j, k = (m for m in range(3) if m != i)
        c = pmul(psub(e[i], e[j], p), psub(e[i], e[k], p), p)
        den = psub(P.A, pmul(e[i], P.B, p), p)
        return _reduced(padd(pmul(e[i], den, p), pmul(c, P.B, p), p), den, p)

    def is_torsion(self, P: Point) -> bool:
        return not self.cubic(P.A, P.B)

    def _shift(self, P: Place) -> int:
        return self.chi if P is None else 0

    def contr(self, P: Point, v: Place) -> Fraction:
        p, s = self.p, self._shift(v)
        vB = _val(P.B, v, p)
        vx = _val(P.A, v, p) - vB + 2 * s
        vf = _val(self.cubic(P.A, P.B), v, p) - 3 * vB + 6 * s
        Q = self.double(P)
        vx2 = _val(Q.A, v, p) - _val(Q.B, v, p) + 2 * s
        vD = _val(self.disc, v, p) + 12 * s
        lam = (Fraction(max(0, -vx2), 2) + Fraction(vD, 12) - Fraction(vf, 2)
               + Fraction(vD, 4)) / 4
        return max(0, -vx) + Fraction(vD, 6) - 2 * lam

    def singular(self, P: Point, v: Place) -> bool:
        p, s = self.p, self._shift(v)
        vB = _val(P.B, v, p)
        vx = _val(P.A, v, p) - vB + 2 * s
        vf = _val(self.cubic(P.A, P.B), v, p) - 3 * vB + 6 * s
        vd = _val(self.dcubic(P.A, P.B), v, p) - 2 * vB + 4 * s
        return vx >= 0 and vf > 0 and vd > 0

    def dot_o(self, P: Point) -> Fraction:
        return Fraction(deg(P.B), 2) + Fraction(max(0, deg(P.A) - deg(P.B) - 2 * self.chi), 2)

    def height(self, P: Point) -> Fraction:
        if self.is_torsion(P):
            return Fraction(0)
        bad = sum((place_degree(v) * self.contr(P, v) for v, f in self.fibers.items()
                   if f.kind != "good"), Fraction(0))
        return 2 * self.chi + 2 * self.dot_o(P) - bad


def _reduced(num: Poly, den: Poly, p: int) -> Point:
    g = pgcd(num, den, p)
    num, den = pdivmod(num, g, p)[0], pdivmod(den, g, p)[0]
    c = pow(den[-1], -1, p)
    return Point(pscale(num, c, p), pscale(den, c, p))


def model(D: Poly, p: int) -> Model:
    e = legendre_roots(D, p)
    fib = fibers(e, p)
    for f in fib.values():
        if f.kind not in MAX_CONTR:
            raise ValueError(f"no height table for {f.kind}")
    disc = sum(f.disc * place_degree(P) for P, f in fib.items())
    return Model(p, e, disc // 12, fib)


def _kernel(rows: list[list[int]], ncols: int) -> list[list[int]]:
    pivots: dict[int, list[int]] = {}
    for r in rows:
        r = r[:]
        for c, pr in pivots.items():
            if r[c]:
                r = [(a + b) % 2 for a, b in zip(r, pr, strict=True)]
        lead = next((c for c in range(ncols) if r[c]), None)
        if lead is None:
            continue
        for c, pr in pivots.items():
            if pr[lead]:
                pivots[c] = [(a + b) % 2 for a, b in zip(pr, r, strict=True)]
        pivots[lead] = r
    free = [c for c in range(ncols) if c not in pivots]
    out = []
    for f in free:
        v = [0] * ncols
        v[f] = 1
        for c, pr in pivots.items():
            v[c] = pr[f]
        out.append(v)
    return out


def selmer_classes(M: Model) -> list[tuple[Poly, Poly]]:
    p, e = M.p, M.e
    S = list(M.fibers)
    gens: list[Poly] = [(nonsquare(p),), *[P for P in S if P is not None]]
    n = len(gens)
    rows: list[list[int]] = []
    for P in S:
        res = [klass(g, P, p) for g in gens]
        for a in _annihilator(local_image(e, P, p)):
            row = [0] * (2 * n)
            for comp in range(2):
                for k, r in enumerate(res):
                    row[comp * n + k] = (a[2 * comp] * r[0] + a[2 * comp + 1] * r[1]) % 2
            rows.append(row)
    basis = _kernel(rows, 2 * n)
    out = []
    for bits in itertools.product((0, 1), repeat=len(basis)):
        v = [sum(b * vec[c] for b, vec in zip(bits, basis, strict=True)) % 2
             for c in range(2 * n)]
        out.append((_prod([g for g, b in zip(gens, v[:n], strict=True) if b], p),
                    _prod([g for g, b in zip(gens, v[n:], strict=True) if b], p)))
    return out


@functools.cache
def _alphas(p: int, size: int = 24) -> tuple[GF, np.ndarray]:
    m = next(m for m in range(1, 8) if p**m >= 4 * size)
    F = field(p, m)
    reps = F.orbits[0]
    return F, reps[reps != 0][:size]


@functools.cache
def _u_table(p: int, top: int) -> tuple[list[Poly], np.ndarray]:
    F, a = _alphas(p)
    dig = [F.dig[F.power(a, j)] for j in range(top + 1)]
    polys, vals = [], []
    for d in range(top + 1):
        C = np.array(list(itertools.product(range(p), repeat=d)), dtype=np.int64).reshape(p**d, d)
        acc = np.broadcast_to(dig[d], (len(C), *dig[d].shape)).copy()
        for j in range(d):
            acc += C[:, j, None, None] * dig[j][None]
        vals.append(F.code(acc))
        polys += [(*map(int, c), 1) for c in C]
    return polys, np.concatenate(vals)


def _sub(F: GF, a, b) -> np.ndarray:
    return F.add(a, F.neg(b))


Klass = tuple[Poly, Poly]


def sqclass(g: Poly, p: int) -> Poly:
    odd = [P for P in primes_dividing(g, p) if local(g, P, p)[0] % 2]
    c = 1 if pow(g[-1], (p - 1) // 2, p) == 1 else nonsquare(p)
    return _prod([(c,), *odd], p)


def delta(M: Model, P: Point) -> Klass:
    p = M.p
    return (sqclass(pmul(psub(P.A, pmul(M.e[0], P.B, p), p), P.B, p), p),
            sqclass(pmul(psub(P.A, pmul(M.e[1], P.B, p), p), P.B, p), p))


def torsion_classes(M: Model) -> set[Klass]:
    p, e = M.p, M.e
    d = {(i, j): psub(e[i], e[j], p) for i in range(3) for j in range(3) if i != j}
    pairs = [(pmul(d[0, 1], d[0, 2], p), d[0, 1]), (d[1, 0], pmul(d[1, 0], d[1, 2], p)),
             (d[2, 0], d[2, 1])]
    return {((1,), (1,)), *((sqclass(a, p), sqclass(b, p)) for a, b in pairs)}


def points(M: Model, bound: int, cls: Klass) -> list[Point]:
    p, (d1, d2) = M.p, cls
    d3 = pmul(d1, d2, p)
    F, a = _alphas(p)
    top = (2 * M.chi + 2 * bound - deg(d1)) // 2
    if top < 0:
        return []
    us, U = _u_table(p, top)
    U2 = F.power(U, 2)
    ev = [F.polyval(f, a)[None] for f in (d1, d2, d3, M.e[1], M.e[2])]
    found = []
    for lam in range(1, (p + 1) // 2):
        A = F.mul(F.mul(lam * lam % p, ev[0]), U2)
        for j in range(bound + 1):
            for w in monics(p, j):
                B = F.power(F.polyval(w, a), 2)[None]
                x1 = F.mul(ev[1], _sub(F, A, F.mul(ev[3], B)))
                x2 = F.mul(ev[2], _sub(F, A, F.mul(ev[4], B)))
                keep = ((F.chi[x1] != -1) & (F.chi[x2] != -1)).all(axis=1)
                for i in np.flatnonzero(keep):
                    du = pmul(d1, us[i], p)
                    if j and deg(pgcd(du, w, p)) > 0:
                        continue
                    P = Point(pscale(pmul(du, us[i], p), lam * lam, p), pmul(w, w, p))
                    c = M.cubic(P.A, P.B)
                    if c and psqrt(c, p) is not None and delta(M, P) == cls:
                        found.append(P)
    return found


@dataclass(frozen=True)
class Generator:
    point: Point
    height: Fraction
    bound: int
    cls: Klass


def generator(M: Model, max_bound: int = 4) -> Generator | None:
    tors = torsion_classes(M)
    classes = sorted(set(M.selmer) - tors, key=lambda c: (deg(c[0]), c))
    for bound in range(max_bound + 1):
        for cls in classes:
            hs = [(M.height(P), P) for P in points(M, bound, cls)]
            if not hs:
                continue
            h, P = min(hs, key=lambda x: x[0])
            if h / 9 < 2 * M.chi + 2 * (bound + 1) - M.max_contr:
                return Generator(P, h, bound, cls)
    return None


def lstar(L: LPoly) -> Fraction:
    c = list(L.c)
    for _ in range(L.order(1)):
        c = _divide(c, L.q)[0]
    return sum((Fraction(ck, L.q**k) for k, ck in enumerate(c)), Fraction(0))


@dataclass(frozen=True)
class RankOne:
    D: Poly
    q: int
    chi: int
    tamagawa: int
    torsion: int
    selmer: int
    lstar: Fraction
    generator: Generator
    sha: Fraction

    @property
    def sha2(self) -> int:
        return self.selmer - 3


def rank_one(t: Twist, max_bound: int = 4) -> RankOne | None:
    if t.rank != 1:
        raise ValueError("analytic rank is not 1")
    q, D = t.L.q, t.D
    M = model(D, q)
    G = generator(M, max_bound)
    if G is None:
        return None
    tam = math.prod(f.c for f in M.fibers.values())
    tors = torsion(D, q)
    ls = lstar(t.L)
    sha = ls * Fraction(q) ** (M.chi - 1) * tors**2 / (tam * G.height)
    return RankOne(D, q, M.chi, tam, tors, selmer_rank(M.e, q), ls, G, sha)
