from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.funcfield import Twist


def vp(x: np.ndarray, p: int, K: int) -> np.ndarray:
    x = np.asarray(x) % p**K
    v = np.zeros(x.shape, dtype=np.int64)
    alive = x != 0
    v[~alive] = K
    while alive.any():
        alive &= x % p == 0
        v[alive] += 1
        x = np.where(alive, x // p, x)
    return v


def _inverse(u: np.ndarray, m: int, phi: int) -> np.ndarray:
    out, base, e = np.ones_like(u), u % m, phi - 1
    while e:
        if e & 1:
            out = out * base % m
        base, e = base * base % m, e >> 1
    return out


def det_valuation(M: np.ndarray, p: int, K: int) -> np.ndarray:
    m = p**K
    if m >= 2**31:
        raise ValueError("p^K too large for int64 products")
    A = np.array(M, dtype=np.int64) % m
    S, n, _ = A.shape
    idx = np.arange(S)
    v = np.zeros(S, dtype=np.int64)
    for c in range(n):
        val = vp(A[:, c:, c], p, K)
        r = np.argmin(val, axis=1) + c
        pv = val[idx, r - c]
        top = A[idx, c].copy()
        A[idx, c] = A[idx, r]
        A[idx, r] = top
        v += pv
        live = pv < K
        scale = np.where(live, p ** np.minimum(pv, K - 1), 1)
        unit = np.where(live, A[idx, c, c] // scale, 1)
        f = (A[:, c + 1:, c] // scale[:, None]) * _inverse(unit, m, m - m // p)[:, None] % m
        f[~live] = 0
        A[:, c + 1:, :] = (A[:, c + 1:, :] - f[:, :, None] * A[:, c, None, :]) % m
    return np.minimum(v, K)


def alternating(n: int, p: int, K: int, size: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    U = np.triu(rng.integers(0, p**K, (size, n, n)), 1)
    return (U - U.transpose(0, 2, 1)) % p**K


def alternating_count(n: int, r: int, q: int) -> int:
    num = math.prod(q ** (n - i) - 1 for i in range(2 * r))
    den = math.prod(q ** (2 * i + 2) - 1 for i in range(r))
    return num * q ** (r * (r - 1)) // den


def nonsingular(n: int, p: int) -> float:
    return math.prod(1 - p ** (1 - 2 * j) for j in range(1, n // 2 + 1))


def divides(p: int, rank: int = 0, terms: int = 60) -> float:
    return 1 - math.prod(1 - p ** (1 - 2 * j) for j in range(rank + 1, rank + terms))


def kernel_law(k: int, p: int, terms: int = 60) -> float:
    tail = math.prod(1 - p ** (1 - 2 * j) for j in range(k + 1, k + terms))
    return p ** (-k * (2 * k - 1)) * tail / math.prod(1 - p ** (-2 * j) for j in range(1, k + 1))


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        raise ZeroDivisionError("not a unit")
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def nonresidue(p: int) -> int:
    return next(g for g in range(2, p) if legendre(g, p) == -1)


def square_class(x: Fraction, p: int) -> tuple[int, int]:
    if x == 0:
        raise ZeroDivisionError("class of zero")
    n, d, v = x.numerator, x.denominator, 0
    while n % p == 0:
        n, v = n // p, v + 1
    while d % p == 0:
        d, v = d // p, v - 1
    return v % 2, legendre(n * d, p)


@dataclass(frozen=True)
class Space:
    n: int
    p: int
    K: int
    disc: int

    @property
    def gram(self) -> np.ndarray:
        g = np.ones(self.n, dtype=np.int64)
        if self.disc == -1:
            g[-1] = nonresidue(self.p)
        return g

    def norm(self, v: np.ndarray) -> np.ndarray:
        return (v * v * self.gram).sum(axis=-1) % self.p**self.K

    def anisotropic(self, chi: int) -> np.ndarray:
        for a in range(self.p):
            for b in range(self.p):
                v = np.zeros(self.n, dtype=np.int64)
                v[0], v[1] = a, b
                Q = int(self.norm(v))
                if Q % self.p and legendre(Q, self.p) == chi:
                    return v
        raise ValueError("no vector of that norm class")


def _reflect(g: np.ndarray, v: np.ndarray, V: Space) -> np.ndarray:
    m = V.p**V.K
    Q = V.norm(v)
    c = 2 * _inverse(Q, m, m - m // V.p) % m
    w = np.einsum("si,sij->sj", v * V.gram, g) % m
    return (g - (c[:, None] * v)[:, :, None] * w[:, None, :] % m) % m


def special_orthogonal(V: Space, size: int, spinor: int | None = None, seed: int = 0,
                       steps: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    m = V.p**V.K
    steps = steps or 2 * (2 * V.n + 6)
    g = np.broadcast_to(np.eye(V.n, dtype=np.int64), (size, V.n, V.n)).copy()
    theta = np.ones(size, dtype=np.int64)
    chi = np.array([0] + [legendre(a, V.p) for a in range(1, V.p)])
    for _ in range(steps):
        v = rng.integers(0, m, (size, V.n))
        bad = V.norm(v) % V.p == 0
        while bad.any():
            v[bad] = rng.integers(0, m, (int(bad.sum()), V.n))
            bad = V.norm(v) % V.p == 0
        theta *= chi[V.norm(v) % V.p]
        g = _reflect(g, v, V)
    if spinor is not None:
        flip = theta != spinor
        if flip.any():
            k = int(flip.sum())
            for chi in (1, -1):
                g[flip] = _reflect(g[flip], np.tile(V.anisotropic(chi), (k, 1)), V)
            theta[flip] = spinor
    return g, theta


def orthogonal_valuations(V: Space, size: int, spinor: int | None = None,
                          seed: int = 0) -> np.ndarray:
    g, _ = special_orthogonal(V, size, spinor, seed)
    eye = np.eye(V.n, dtype=np.int64)
    return det_valuation((eye - g) % V.p**V.K, V.p, V.K)


def first_step(p: int, terms: int = 60) -> float:
    return (1 - divides(p, terms=terms)) * p / (p * p - 1)


def central_integer(t: Twist) -> int:
    v = t.L.value(Fraction(1, t.L.q)) * t.L.q**t.L.N
    if v.denominator != 1:
        raise ArithmeticError("q^N L(1/q) is not an integer")
    return int(v)


def valuation(x: int, p: int) -> int:
    if x == 0:
        raise ZeroDivisionError("valuation of zero")
    v = 0
    while x % p == 0:
        x, v = x // p, v + 1
    return v


@dataclass(frozen=True)
class OddSha:
    q: int
    N: int
    p: int
    twists: int
    observed: np.ndarray
    orthogonal: np.ndarray
    bklpr: np.ndarray


def _buckets(v: np.ndarray) -> np.ndarray:
    return np.array([np.mean(v == 0), np.mean(v == 2), np.mean(v >= 4)])


def odd_sha(fam: list[Twist], p: int, size: int = 20_000, K: int = 6, seed: int = 0) -> OddSha:
    q, N = fam[0].L.q, fam[0].L.N
    if q % p == 0 or p == 2:
        raise ValueError("p must be odd and prime to q")
    v = np.array([valuation(central_integer(t), p) for t in fam if t.L.eps == 1 and t.rank == 0])
    model = orthogonal_valuations(Space(2 * (N // 2), p, K, 1), size, seed=seed)
    c = 1 - divides(p)
    b = first_step(p)
    return OddSha(q, N, p, len(v), _buckets(v), _buckets(model), np.array([c, b, 1 - c - b]))
