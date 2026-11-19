from __future__ import annotations

import functools
import itertools

import numpy as np

Poly = tuple[int, ...]


def norm(a, p: int) -> Poly:
    out = [x % p for x in a]
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def deg(a: Poly) -> int:
    return len(a) - 1


def padd(a: Poly, b: Poly, p: int) -> Poly:
    n = max(len(a), len(b))
    return norm([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)], p)


def pscale(a: Poly, c: int, p: int) -> Poly:
    return norm([c * x for x in a], p)


def psub(a: Poly, b: Poly, p: int) -> Poly:
    return padd(a, pscale(b, -1, p), p)


def pmul(a: Poly, b: Poly, p: int) -> Poly:
    if not a or not b:
        return ()
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return norm(out, p)


def pdivmod(a: Poly, b: Poly, p: int) -> tuple[Poly, Poly]:
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    r = list(a)
    inv = pow(b[-1], -1, p)
    q = [0] * max(len(a) - len(b) + 1, 0)
    for k in range(len(q) - 1, -1, -1):
        c = r[k + len(b) - 1] * inv % p
        q[k] = c
        for j, y in enumerate(b):
            r[k + j] -= c * y
    return norm(q, p), norm(r, p)


def monic(a: Poly, p: int) -> Poly:
    return pscale(a, pow(a[-1], -1, p), p) if a else ()


def pgcd(a: Poly, b: Poly, p: int) -> Poly:
    while b:
        a, b = b, pdivmod(a, b, p)[1]
    return monic(a, p)


def ppowmod(a: Poly, e: int, f: Poly, p: int) -> Poly:
    out, a = (1,), pdivmod(a, f, p)[1]
    while e:
        if e & 1:
            out = pdivmod(pmul(out, a, p), f, p)[1]
        a = pdivmod(pmul(a, a, p), f, p)[1]
        e >>= 1
    return out


def deriv(a: Poly, p: int) -> Poly:
    return norm([i * a[i] for i in range(1, len(a))], p)


def peval(a: Poly, x: int, p: int) -> int:
    out = 0
    for c in reversed(a):
        out = (out * x + c) % p
    return out


def squarefree(a: Poly, p: int) -> bool:
    return deg(pgcd(a, deriv(a, p), p)) == 0


def distinct_roots(a: Poly, p: int) -> int:
    f = monic(a, p)
    if deg(f) < 1:
        return 0
    h, count, total = (0, 1), {}, 0
    for d in range(1, deg(f) + 1):
        h = ppowmod(h, p, f, p)
        g = deg(pgcd(f, psub(h, (0, 1), p), p))
        count[d] = (g - sum(e * count[e] for e in count if d % e == 0 and e < d)) // d
        total += d * count[d]
    return total


def _split(g: Poly, d: int, p: int, rng: np.random.Generator) -> list[Poly]:
    if deg(g) == d:
        return [g]
    e = (p**d - 1) // 2
    while True:
        r = norm([int(c) for c in rng.integers(0, p, deg(g))] + [1], p)
        h = pgcd(g, psub(ppowmod(r, e, g, p), (1,), p), p)
        if 0 < deg(h) < deg(g):
            return _split(h, d, p, rng) + _split(pdivmod(g, h, p)[0], d, p, rng)


def factor(a: Poly, p: int) -> list[Poly]:
    f = monic(a, p)
    if not squarefree(f, p):
        raise ValueError("factor expects a squarefree polynomial")
    rng = np.random.default_rng(p)
    out, h, d = [], (0, 1), 0
    while deg(f) >= 2 * (d + 1):
        d += 1
        h = ppowmod(h, p, f, p)
        g = pgcd(f, psub(h, (0, 1), p), p)
        if deg(g) > 0:
            out += _split(g, d, p, rng)
            f = pdivmod(f, g, p)[0]
            h = pdivmod(h, f, p)[1]
    if deg(f) > 0:
        out.append(f)
    return sorted(out, key=lambda P: (deg(P), P[::-1]))


def primes_dividing(a: Poly, p: int) -> list[Poly]:
    f, out = monic(a, p), set()
    while deg(f) > 0:
        df = deriv(f, p)
        if not df:
            f = f[::p]
            continue
        r = pdivmod(f, pgcd(f, df, p), p)[0]
        for P in factor(r, p) if deg(r) > 0 else []:
            out.add(P)
            while not pdivmod(f, P, p)[1]:
                f = pdivmod(f, P, p)[0]
    return sorted(out, key=lambda P: (deg(P), P[::-1]))


def monics(p: int, n: int):
    for low in itertools.product(range(p), repeat=n):
        yield (*low, 1)


def is_prime(p: int) -> bool:
    return p > 1 and all(p % d for d in range(2, int(p**0.5) + 1))


def _prime_factors(n: int) -> list[int]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    return out + [n] if n > 1 else out


def _primitive(p: int, m: int) -> Poly:
    q = p**m
    roots = [g for g in range(1, p) if all(pow(g, (p - 1) // r, p) != 1
                                           for r in _prime_factors(p - 1))]
    for f in monics(p, m):
        if (-1) ** m * f[0] % p not in roots:
            continue
        if ppowmod((0, 1), q - 1, f, p) == (1,) and all(
                ppowmod((0, 1), (q - 1) // r, f, p) != (1,) for r in _prime_factors(q - 1)):
            return f
    raise RuntimeError(f"no primitive polynomial of degree {m} over F_{p}")


def _exp_table(p: int, m: int) -> np.ndarray:
    low = _primitive(p, m)[:-1]
    A = np.zeros((m, m), dtype=np.int64)
    A[1:, :-1] = np.eye(m - 1, dtype=np.int64)
    A[:, -1] = [-c % p for c in low]
    q = p**m
    d = np.eye(1, m, dtype=np.int64)
    while len(d) < q - 1:
        d = np.vstack([d, (d @ A.T) % p])
        A = (A @ A) % p
    return d[:q - 1] @ (p ** np.arange(m, dtype=np.int64))


class GF:
    def __init__(self, p: int, m: int):
        if p == 2 or not is_prime(p):
            raise ValueError("p must be an odd prime")
        if m < 1:
            raise ValueError("m must be >= 1")
        self.p, self.m, self.q = p, m, p**m
        self.pw = p ** np.arange(m, dtype=np.int64)
        self.dig = (np.arange(self.q, dtype=np.int64)[:, None] // self.pw) % p
        self.exp = _exp_table(p, m)
        self.log = np.zeros(self.q, dtype=np.int64)
        self.log[self.exp] = np.arange(self.q - 1)
        self.chi = np.zeros(self.q, dtype=np.int8)
        self.chi[self.exp] = 1 - 2 * (np.arange(self.q - 1) % 2)

    def __repr__(self) -> str:
        return f"GF({self.p}^{self.m})"

    def code(self, digits) -> np.ndarray:
        return (np.asarray(digits) % self.p) @ self.pw

    def add(self, a, b) -> np.ndarray:
        return self.code(self.dig[a] + self.dig[b])

    def neg(self, a) -> np.ndarray:
        return self.code(-self.dig[a])

    def mul(self, a, b) -> np.ndarray:
        a, b = np.asarray(a), np.asarray(b)
        out = self.exp[(self.log[a] + self.log[b]) % (self.q - 1)]
        return np.where((a == 0) | (b == 0), 0, out)

    def power(self, a, k: int) -> np.ndarray:
        a = np.asarray(a)
        if k == 0:
            return np.ones_like(a)
        return np.where(a == 0, 0, self.exp[(self.log[a] * k) % (self.q - 1)])

    def frobenius(self, a) -> np.ndarray:
        return self.power(a, self.p)

    def degrees(self, a) -> np.ndarray:
        a = np.asarray(a)
        s, size = a, np.zeros(a.shape, dtype=np.int64)
        for k in range(1, self.m + 1):
            s = self.frobenius(s)
            size[(size == 0) & (s == a)] = k
        return size

    def polyval(self, c: Poly, t) -> np.ndarray:
        t = np.asarray(t)
        acc = np.zeros((*t.shape, self.m), dtype=np.int64)
        for j, cj in enumerate(c):
            if cj:
                acc += cj * self.dig[self.power(t, j)]
        return self.code(acc)

    @functools.cached_property
    def orbits(self) -> tuple[np.ndarray, np.ndarray]:
        cur = rep = np.arange(self.q)
        for _ in range(self.m - 1):
            cur = self.frobenius(cur)
            rep = np.minimum(rep, cur)
        return np.unique(rep, return_inverse=True)


@functools.cache
def field(p: int, m: int) -> GF:
    return GF(p, m)


def minimal_polys(q: int, d: int) -> tuple[np.ndarray, list[Poly]]:
    F = field(q, d)
    reps, inv = F.orbits
    alpha = reps[np.bincount(inv) == d]
    c = np.zeros((len(alpha), d + 1), dtype=np.int64)
    c[:, 0] = 1
    beta = alpha
    for _ in range(d):
        shifted = np.zeros_like(c)
        shifted[:, 1:] = c[:, :-1]
        c = F.add(shifted, F.neg(F.mul(c, beta[:, None])))
        beta = F.frobenius(beta)
    if np.any(c >= q):
        raise ArithmeticError("minimal polynomial outside the prime field")
    return alpha, [tuple(int(x) for x in row) for row in c]


def cubic_sums(F: GF, A, B, C) -> np.ndarray:
    A, B, C = (np.atleast_1d(np.asarray(v, dtype=np.int64)) for v in (A, B, C))
    x = np.arange(F.q)
    x2, x3 = F.power(x, 2), F.dig[F.power(x, 3)]
    rows = max(1, 2**22 // (F.q * F.m))
    out = np.empty(len(A), dtype=np.int64)
    for i in range(0, len(A), rows):
        s = slice(i, i + rows)
        d = (x3 + F.dig[F.mul(A[s, None], x2)] + F.dig[F.mul(B[s, None], x)]
             + F.dig[C[s]][:, None, :])
        out[s] = F.chi[F.code(d)].sum(axis=1, dtype=np.int64)
    return out


def correlate(F: GF, g: np.ndarray, h: np.ndarray) -> np.ndarray:
    shape = (F.p,) * F.m
    G = np.fft.fftn(g.reshape(shape))
    H = np.fft.fftn(h.reshape(shape))
    return np.rint(np.fft.ifftn(G * np.conj(H)).real).astype(np.int64).reshape(F.q)


def character_sums(F: GF, Ds: list[Poly], w: np.ndarray) -> np.ndarray:
    reps = F.orbits[0]
    P = np.stack([F.power(reps, j) for j in range(len(Ds[0]))])
    digits = [F.dig[P, j].astype(np.float64) for j in range(F.m)]
    C = np.array(Ds, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    chunk = max(1, 2**22 // len(reps))
    out = np.empty(len(Ds), dtype=np.int64)
    for i in range(0, len(Ds), chunk):
        rows = C[i:i + chunk]
        code = np.zeros((len(rows), len(reps)), dtype=np.int64)
        for d, pw in zip(digits, F.pw, strict=True):
            code += np.fmod(rows @ d, F.p).astype(np.int64) * pw
        out[i:i + chunk] = np.rint(F.chi[code] @ w).astype(np.int64)
    return out
