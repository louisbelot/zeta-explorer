from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath
import numpy as np

from zeta_explorer.characters import factor, kronecker
from zeta_explorer.explicit import primes

CONVEXITY = {2: "2.97655", 4: "5.13781", 8: "9.49562"}
J0 = "1.912"
TOL = 1e-7


@dataclass(frozen=True)
class Constants:
    r: float
    phi: float
    E: float


def constants(r: float, dps: int = 30) -> Constants:
    with mpmath.workdps(dps):
        mpf, pi, log = mpmath.mpf, mpmath.pi, mpmath.log
        r = mpf(r)
        R = r + mpf(7) / 8
        C = {k: mpf(v) for k, v in CONVEXITY.items()}
        th = [pi / 2] + [mpmath.acos((r + mpf(o) / 8) / R) for o in (1, 2, 4, 6)] + [mpf(0)]
        s = [mpmath.sin(x) for x in th]
        c = [mpmath.cos(x) for x in th]

        def A(n, m):
            return s[m] - s[n]

        def B(n, m):
            return (th[n] - th[m]) / 2 + (s[n] * c[n] - s[m] * c[m]) / 2

        def weight(shift, n, m):
            return (r + shift) / (pi * R) * A(n, m) + B(n, m) / pi

        phi = r * A(1, 5) / (pi * R) + B(0, 1) / (pi * (8 * r + 1)) + B(1, 5) / pi
        rho = -2 * (A(0, 1) / (pi * R) + B(0, 1) / (pi * (r + mpf(1) / 8)))
        K = (2 * log(C[8]) + log(2 * r + 2) / 8) * B(0, 1) / (pi * (r + mpf(1) / 8))
        arc = abs(2 + r - R * mpmath.expjpi(th[3] / pi))
        K += 2 * (log(arc) - log(2 * pi)) * weight(mpf(1) / 2, 3, 5)
        for k in (2, 3):
            a = mpf(2) ** (k - 4)
            u = 2 ** (6 - k) * log(C[2 ** (4 - k)]) + 2 * log(mpf(15) / 8 + a + 2 * r)
            v = 2 ** (6 - k) * log(C[2 ** (5 - k)]) + log(mpf(15) / 8 + a + 2 * r)
            K += u * weight(a / 2, k - 1, k) - v * weight(a, k - 1, k)
        for k in (4, 5):
            a = mpf(2) ** (3 - k)
            u = 2 ** (k - 1) * log(C[2 ** (k - 3)]) + 2 * log(mpf(7) / 8 + a)
            v = 2 ** (k - 1) * log(C[2 ** (k - 2)]) + log(mpf(7) / 8 + a)
            K += v * weight(1 - a, k - 1, k) - u * weight(1 - a / 2, k - 1, k)
        K += mpf(J0) / (pi * R) - log(2 * pi) + mpmath.euler / 2 + 1
        K += mpmath.digamma((3 + r) / 2) / 2
        E = K + rho * log(1 + 1 / r)
        return Constants(float(r), float(phi), float(E))


def bound(q: np.ndarray | int, k: Constants, c: float = 0.2) -> np.ndarray:
    L = np.log(np.asarray(q, dtype=float))
    r = k.r
    R = r + 7 / 8
    return c / (r * (r * L + c)) + (r + c / L) / R**2 + k.phi * L + k.E


def fundamental(X: int) -> np.ndarray:
    sq = np.ones(X + 1, dtype=bool)
    sq[0] = False
    for p in primes(math.isqrt(X)):
        sq[p * p :: p * p] = False
    m = np.flatnonzero(sq)
    m4 = m[m <= X // 4]
    out = np.concatenate(
        [
            m[m % 4 == 1],
            4 * m4[np.isin(m4 % 4, (2, 3))],
            -m[m % 4 == 3],
            -4 * m4[np.isin(m4 % 4, (1, 2))],
        ]
    )
    out = out[out != 1]
    return out[np.argsort(np.abs(out) * 2 + (out < 0), kind="stable")]


def nth_primes(n: int) -> np.ndarray:
    if n < 6:
        return primes(13)[:n]
    x = int(n * (math.log(n) + math.log(math.log(n)))) + 10
    return primes(x)[:n]


def legendre_table(p: int) -> np.ndarray:
    t = np.full(p, -1, dtype=np.int8)
    t[0] = 0
    t[np.arange(1, (p + 1) // 2, dtype=np.int64) ** 2 % p] = 1
    return t


def kronecker_column(d: np.ndarray, p: int) -> np.ndarray:
    if p == 2:
        return np.select([d % 2 == 0, np.isin(d % 8, (1, 7))], [0, 1], -1).astype(np.int8)
    if 32 * d.size > p:
        return legendre_table(p)[d % p]
    base = d % p
    out = np.ones_like(base)
    e = (p - 1) // 2
    while e:
        if e & 1:
            out = out * base % p
        base = base * base % p
        e >>= 1
    return np.where(out == p - 1, -1, out).astype(np.int8)


def weights(ps: np.ndarray, sigma: float) -> np.ndarray:
    lp = np.log(ps)
    x = ps.astype(float) ** sigma
    return np.stack([2 * lp / (x * x - 1), lp / (x - 1), 2 * lp / (x - 1)], axis=1)


def full_sum(d: int, sigma: float, dps: int = 30) -> mpmath.mpf:
    vals = [kronecker(d, a) for a in range(abs(d))]
    with mpmath.workdps(dps):
        s = mpmath.mpf(sigma)
        lz = mpmath.zeta(s, derivative=1) / mpmath.zeta(s)
        return -lz - mpmath.dirichlet(s, vals, 1) / mpmath.dirichlet(s, vals)


@dataclass
class Scan:
    d: np.ndarray
    k: Constants
    c: float
    needed: np.ndarray
    margin: np.ndarray
    n_primes: int

    @property
    def failed(self) -> np.ndarray:
        return self.d[self.needed == 0]


def scan(
    d: np.ndarray,
    k: Constants,
    n_primes: int,
    c: float = 0.2,
    ps: np.ndarray | None = None,
) -> Scan:
    d = np.asarray(d, dtype=np.int64)
    if ps is None:
        ps = nth_primes(n_primes)
    ps = ps[:n_primes]
    w = weights(ps, 1 + k.r)
    target = bound(np.abs(d), k, c) + TOL
    S = np.zeros(d.size)
    needed = np.zeros(d.size, dtype=np.int64)
    margin = np.zeros(d.size)
    live = np.arange(d.size)
    for i, p in enumerate(ps):
        if not live.size:
            break
        S[live] += w[i][kronecker_column(d[live], int(p)) + 1]
        done = S[live] > target[live]
        hit = live[done]
        needed[hit] = i + 1
        margin[hit] = S[hit] - target[hit] + TOL
        live = live[~done]
    margin[live] = S[live] - target[live] + TOL
    return Scan(d, k, c, needed, margin, len(ps))


def character_table(d: int) -> np.ndarray:
    q = abs(d)
    n = np.arange(q, dtype=np.int64)
    t = np.ones(q, dtype=np.int8)
    for p, e in factor(q):
        if p == 2:
            if e == 3:
                t *= np.where(np.isin(n % 8, (1, 7)), 1, -1).astype(np.int8)
            if (d >> e) % 4 == 3:
                t *= np.where(n % 4 == 1, 1, -1).astype(np.int8)
            t *= (n % 2).astype(np.int8)
        else:
            t *= legendre_table(p)[n % p]
    return t


def scan_one(
    d: int, k: Constants, ps: np.ndarray, w: np.ndarray, c: float = 0.2, chunk: int = 1 << 15
) -> tuple[int, float]:
    tab = character_table(d)
    target = float(bound(abs(d), k, c)) + TOL
    S, lo = 0.0, 0
    while lo < len(ps):
        hi = min(lo + chunk, len(ps))
        chi = tab[ps[lo:hi] % abs(d)]
        part = S + np.cumsum(w[np.arange(lo, hi), chi + 1])
        i = int(np.searchsorted(part, target, side="right"))
        if i < hi - lo:
            return lo + i + 1, float(part[i]) - target + TOL
        S, lo, chunk = float(part[-1]), hi, 2 * chunk
    return 0, S - target + TOL
