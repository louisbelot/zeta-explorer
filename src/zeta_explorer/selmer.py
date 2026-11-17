from __future__ import annotations

import functools
import math
from dataclasses import dataclass

import numpy as np

from zeta_explorer.bsd import Block, klass, legendre_roots, local_image, nonsquare, selmer_dim
from zeta_explorer.funcfield import Twist
from zeta_explorer.gf import Poly, deg, primes_dividing


def poonen_rains(n: int, p: int = 2, terms: int = 80) -> float:
    c = math.prod(1 / (1 + p ** (-j)) for j in range(terms))
    return c * math.prod(p / (p**j - 1) for j in range(1, n + 1))


def equilibrium(n: int, parity: int, p: int = 2) -> float:
    return 2 * poonen_rains(n, p) if n % 2 == parity % 2 else 0.0


def markov(size: int, p: int = 2) -> np.ndarray:
    M = np.zeros((size, size))
    for n in range(size):
        up = p ** (-n)
        if n + 1 < size:
            M[n, n + 1] = up
        if n:
            M[n, n - 1] = 1 - up
    return M


def chain(start: int, steps: int, size: int = 40, p: int = 2) -> np.ndarray:
    x = np.zeros(size)
    x[start] = 1.0
    M = markov(size, p)
    for _ in range(steps):
        x = x @ M
    return x


Bit = tuple[int, int]


def _mul(*xs: Bit) -> Bit:
    return sum(x[0] for x in xs) % 2, sum(x[1] for x in xs) % 2


@dataclass(frozen=True)
class Symbols:
    q: int
    degs: tuple[int, ...]
    at0: tuple[int, ...]
    at1: tuple[int, ...]
    mutual: tuple[tuple[int, ...], ...]

    @property
    def omega(self) -> int:
        return len(self.degs)

    @property
    def n(self) -> int:
        return sum(self.degs)

    def sign(self, d: int) -> int:
        return (self.q - 1) // 2 * d % 2

    def at_t(self, i: int) -> int:
        return (self.at0[i] + self.sign(self.degs[i])) % 2

    def at_t1(self, i: int) -> int:
        return (self.at1[i] + self.sign(self.degs[i])) % 2


def symbols(D: Poly, q: int) -> Symbols:
    P = primes_dividing(D, q)
    res = [[klass(Pj, Pi, q)[1] if i != j else 0 for j, Pj in enumerate(P)]
           for i, Pi in enumerate(P)]
    return Symbols(q, tuple(deg(x) for x in P),
                   tuple(klass(x, (0, 1), q)[1] for x in P),
                   tuple(klass(x, (q - 1, 1), q)[1] for x in P),
                   tuple(map(tuple, res)))


def random_symbols(q: int, degs: tuple[int, ...], rng: np.random.Generator) -> Symbols:
    w = len(degs)
    up = np.triu(rng.integers(0, 2, (w, w)), 1)
    sign = np.outer(degs, degs) * ((q - 1) // 2) % 2
    M = (up + up.T + np.tril(sign, -1)) % 2
    return Symbols(q, tuple(degs), tuple(map(int, rng.integers(0, 2, w))),
                   tuple(map(int, rng.integers(0, 2, w))), tuple(map(tuple, M.tolist())))


@dataclass(frozen=True)
class Tables:
    q: int
    t: tuple[list[tuple[int, ...]], list[tuple[int, ...]]]
    t1: tuple[list[tuple[int, ...]], list[tuple[int, ...]]]
    inf: tuple[list[tuple[int, ...]], list[tuple[int, ...]]]


@functools.cache
def tables(q: int) -> Tables:
    c = (1, nonsquare(q))
    t = tuple(local_image(legendre_roots((x,), q), (0, 1), q) for x in c)
    t1 = tuple(local_image(legendre_roots((x,), q), (q - 1, 1), q) for x in c)
    inf = tuple(local_image(legendre_roots(D, q), None, q) for D in ((1,), (1, 1)))
    return Tables(q, t, t1, inf)


def symbolic_rank(S: Symbols) -> int:
    T, w, minus = tables(S.q), S.omega, S.sign(1)
    gens_t = [(0, 1), (1, 0), (0, minus), *[(0, b) for b in S.at0]]
    gens_t1 = [(0, 1), (0, 0), (1, 0), *[(0, b) for b in S.at1]]
    gens_inf = [(0, 1), (1, 0), (1, 0), *[(d % 2, 0) for d in S.degs]]
    blocks: list[Block] = [
        (T.t[sum(S.at0) % 2], gens_t),
        (T.t1[sum(S.at1) % 2], gens_t1),
        (T.inf[S.n % 2], gens_inf),
    ]
    for i, d in enumerate(S.degs):
        rest = sum(S.mutual[i]) % 2
        neg, tb, t1b = S.sign(d), S.at_t(i), S.at_t1(i)
        d01, d02 = (1, (neg + rest) % 2), (1, (neg + tb + rest) % 2)
        d10, d12 = (1, rest), (1, (rest + neg + t1b) % 2)
        W = [v for v in (_mul(d01, d02) + d01, d10 + _mul(d10, d12)) if any(v)]
        gens = [(0, d % 2), (0, tb), (0, t1b),
                *[(1, 0) if j == i else (0, S.mutual[i][j]) for j in range(w)]]
        blocks.append((W, gens))
    return selmer_dim(blocks)


def histogram(s: list[int] | np.ndarray, size: int = 8) -> np.ndarray:
    h = np.bincount(np.minimum(np.asarray(s, dtype=np.int64), size - 1), minlength=size)
    return h / max(1, h.sum())


def model(degs: list[tuple[int, ...]], q: int, draws: int = 4, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.array([symbolic_rank(random_symbols(q, d, rng)) - 2 for d in degs
                     for _ in range(draws)])


@dataclass(frozen=True)
class Law:
    q: int
    n: int
    twists: int
    observed: np.ndarray
    predicted: np.ndarray
    mean: float
    mean_model: float

    def given(self, parity: int, which: str = "observed") -> np.ndarray:
        h = getattr(self, which).copy()
        h[np.arange(len(h)) % 2 != parity] = 0
        return h / h.sum()

    @property
    def tv(self) -> float:
        return 0.5 * float(np.abs(self.observed - self.predicted).sum())


def law(fam: list[Twist], draws: int = 4, seed: int = 0, size: int = 8) -> Law:
    q = fam[0].L.q
    sym = [symbols(t.D, q) for t in fam]
    s = np.array([symbolic_rank(S) - 2 for S in sym])
    if np.any(s % 2 != np.array([t.rank % 2 for t in fam])):
        raise ArithmeticError("Selmer parity differs from the rank parity")
    m = model([S.degs for S in sym], q, draws, seed)
    return Law(q, deg(fam[0].D), len(fam), histogram(s, size), histogram(m, size),
               float(np.mean(2.0**s)), float(np.mean(2.0**m)))
