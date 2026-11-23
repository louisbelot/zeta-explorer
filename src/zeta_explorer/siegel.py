from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath
import numpy as np
from flint import acb, arb, ctx, dirichlet_char

from zeta_explorer.characters import factor, kronecker, totient
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


def resum(d: int, r: float, n: int, bits: int = 128) -> arb:
    old = ctx.prec
    ctx.prec = bits
    try:
        sigma = 1 + arb(r)
        S = arb(0)
        for p in nth_primes(n):
            chi = int(kronecker_column(np.array([d]), int(p))[0])
            lp = arb(int(p)).log()
            x = arb(int(p)) ** sigma
            S += lp / (x - 1) + (chi * lp / (x - chi) if chi else 0)
        return S
    finally:
        ctx.prec = old


def flint_character(d: int) -> dirichlet_char:
    q = abs(d)
    n = np.arange(1, q + 1, dtype=np.int64)
    tab = character_table(d)
    units = np.flatnonzero(tab)
    for m in n[n * n % q == 1]:
        chi = dirichlet_char(q, int(m))
        if not chi.is_primitive():
            continue
        if all((-1 if chi.chi_exponent(int(x)) else 1) == tab[x] for x in units):
            return chi
    raise ValueError(f"no primitive real character for {d}")


def positive_on(d: int, lo: float, hi: float, bits: int = 64, depth: int = 30) -> bool:
    if lo <= 0:
        raise ValueError("lo must be positive")
    chi = flint_character(d)
    old = ctx.prec
    ctx.prec = bits
    try:
        M = arb(totient(abs(d))) / (arb.const_e() * arb(lo))
        stack = [(lo, hi, 0)]
        while stack:
            a, b, k = stack.pop()
            m = (a + b) / 2
            v = chi.l(acb(m)).real
            if v > M * (arb(m) - arb(a)) and v > M * (arb(b) - arb(m)):
                continue
            if v < 0 or k == depth:
                return False
            stack += [(a, m, k + 1), (m, b, k + 1)]
        return True
    finally:
        ctx.prec = old


PASSES = ((1.6, 60_000), (1.3, 5_000_000), (3.5, 5_000_000))


@dataclass
class Verification:
    X: int
    c: float
    d: np.ndarray
    stage: np.ndarray
    needed: np.ndarray
    margin: np.ndarray

    @property
    def ok(self) -> bool:
        return bool((self.stage > 0).all())

    def counts(self) -> list[int]:
        return [int((self.stage == i).sum()) for i in range(1, self.stage.max() + 1)]


def verify(
    X: int, c: float = 0.2, passes: tuple[tuple[float, int], ...] = PASSES
) -> Verification:
    d = fundamental(X)
    d = d[np.abs(d) >= 3]
    stage = np.zeros(d.size, dtype=np.int64)
    needed = np.zeros(d.size, dtype=np.int64)
    margin = np.zeros(d.size)
    ps = nth_primes(max(n for _, n in passes))
    live = np.arange(d.size)
    for j, (lam, n) in enumerate(passes, start=1):
        k = constants(lam / (10 * math.log(10)))
        if j == 1:
            res = scan(d[live], k, n, c, ps)
            got, gap = res.needed, res.margin
        else:
            w = weights(ps[:n], 1 + k.r)
            pairs = [scan_one(int(x), k, ps[:n], w, c) for x in d[live]]
            got = np.array([a for a, _ in pairs], dtype=np.int64)
            gap = np.array([b for _, b in pairs])
        hit = got > 0
        stage[live[hit]] = j
        needed[live[hit]] = got[hit]
        margin[live[hit]] = gap[hit]
        live = live[~hit]
    for i in live:
        q = abs(int(d[i]))
        if positive_on(int(d[i]), 1 - c / math.log(q), 1.0):
            stage[i] = len(passes) + 1
    return Verification(X, c, d, stage, needed, margin)
