from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np


def factor(n: int) -> list[tuple[int, int]]:
    if n < 1:
        raise ValueError("n must be >= 1")
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
        p += 1 if p == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def totient(n: int) -> int:
    return math.prod((p - 1) * p ** (e - 1) for p, e in factor(n))


def conrey_generator(p: int) -> int:
    if p == 2 or factor(p) != [(p, 1)]:
        raise ValueError("p must be an odd prime")
    primes = [r for r, _ in factor(p - 1)]
    for g in range(2, p):
        if all(pow(g, (p - 1) // r, p) != 1 for r in primes) and pow(g, p - 1, p * p) != 1:
            return g
    raise AssertionError("unreachable")


@dataclass(frozen=True, eq=False)
class Character:
    q: int
    num: np.ndarray
    den: int

    def __post_init__(self):
        g = math.gcd(self.den, *self.num[self.num >= 0].tolist())
        if g > 1:
            object.__setattr__(self, "num", np.where(self.num >= 0, self.num // g, -1))
            object.__setattr__(self, "den", self.den // g)

    def __call__(self, m: int) -> complex:
        return complex(self.values[m % self.q])

    @cached_property
    def values(self) -> np.ndarray:
        out = np.zeros(self.q, dtype=complex)
        ok = self.num >= 0
        out[ok] = root_of_unity(self.num[ok], self.den)
        return out

    @property
    def order(self) -> int:
        return self.den

    @property
    def is_real(self) -> bool:
        return self.den <= 2

    @property
    def is_principal(self) -> bool:
        return self.den == 1

    @property
    def parity(self) -> int:
        return 0 if self.q < 3 or self.num[self.q - 1] == 0 else 1

    @cached_property
    def conductor(self) -> int:
        m = np.arange(self.q)
        for f in sorted(d for d in range(1, self.q + 1) if self.q % d == 0):
            if np.all(self.num[(m % f == 1 % f) & (self.num >= 0)] == 0):
                return f
        raise AssertionError("unreachable")

    @property
    def is_primitive(self) -> bool:
        return self.conductor == self.q

    def primitive(self) -> Character:
        f = self.conductor
        if f == self.q:
            return self
        num = np.full(f, -1, dtype=np.int64)
        for r in range(f):
            if math.gcd(r, f) == 1:
                m = next(m for m in range(r, self.q + f, f) if math.gcd(m, self.q) == 1)
                num[r] = self.num[m % self.q]
        return Character(f, num, self.den)

    def conj(self) -> Character:
        return Character(self.q, np.where(self.num >= 0, -self.num % self.den, -1), self.den)

    def __mul__(self, other: Character) -> Character:
        if other.q != self.q:
            raise ValueError("characters need the same modulus")
        den = math.lcm(self.den, other.den)
        num = self.num * (den // self.den) + other.num * (den // other.den)
        return Character(self.q, np.where(self.num >= 0, num % den, -1), den)

    def gauss_sum(self) -> complex:
        a = np.arange(self.q)
        return complex(np.sum(self.values * root_of_unity(a, self.q)))

    def root_number(self) -> complex:
        if not self.is_primitive:
            raise ValueError("the root number is defined for primitive characters")
        return self.gauss_sum() / (1j**self.parity * math.sqrt(self.q))


def root_of_unity(k: np.ndarray, n: int) -> np.ndarray:
    k = np.asarray(k) % n
    out = np.exp(2j * np.pi * k / n)
    quarter = (4 * k) % n == 0
    out[quarter] = np.array([1, 1j, -1, -1j])[(4 * k[quarter]) // n]
    return out


def conrey(q: int, n: int) -> Character:
    if q < 1 or not 1 <= n <= max(q, 1) or math.gcd(n, q) != 1:
        raise ValueError("need q >= 1 and 1 <= n <= q with gcd(n, q) = 1")
    parts = [_conrey_prime_power(p, e, n % p**e) for p, e in factor(q)]
    den = math.lcm(1, *(d for _, d in parts))
    num = np.zeros(q, dtype=np.int64)
    m = np.arange(q)
    for (p, e), (part, d) in zip(factor(q), parts, strict=True):
        num += part[m % p**e] * (den // d)
    coprime = np.gcd(m, q) == 1
    return Character(q, np.where(coprime, num % den, -1), den)


def _conrey_prime_power(p: int, e: int, n: int) -> tuple[np.ndarray, int]:
    pe = p**e
    if p == 2:
        if e == 1:
            return np.zeros(2, dtype=np.int64), 1
        den = 2 ** (e - 1)
        eps, idx = _two_adic_log(e)
        both = (eps == -1) & (eps[n] == -1)
        return (both * 2 ** (e - 2) + 2 * idx * idx[n]) % den, den
    g = conrey_generator(p)
    phi = pe - pe // p
    idx = np.zeros(pe, dtype=np.int64)
    x = 1
    for k in range(phi):
        idx[x] = k
        x = x * g % pe
    return idx * idx[n] % phi, phi


def _two_adic_log(e: int) -> tuple[np.ndarray, np.ndarray]:
    pe = 2**e
    eps = np.zeros(pe, dtype=np.int64)
    idx = np.zeros(pe, dtype=np.int64)
    x = 1
    for k in range(max(1, pe // 4)):
        eps[x], idx[x] = 1, k
        eps[pe - x], idx[pe - x] = -1, k
        x = x * 5 % pe
    return eps, idx


def characters(q: int) -> list[Character]:
    return [conrey(q, n) for n in range(1, max(q, 2)) if math.gcd(n, q) == 1]


def primitive_characters(q: int) -> list[tuple[int, Character]]:
    out = []
    for n in range(1, max(q, 2)):
        if math.gcd(n, q) == 1 and (chi := conrey(q, n)).is_primitive:
            out.append((n, chi))
    return out


def kronecker(a: int, n: int) -> int:
    if n < 0:
        raise ValueError("n must be >= 0")
    if n == 0:
        return 1 if abs(a) == 1 else 0
    out = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            out = -out
    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                out = -out
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            out = -out
        a %= n
    return out if n == 1 else 0


def is_fundamental(d: int) -> bool:
    if d in (0, 1):
        return False
    if d % 4 == 1:
        return _squarefree(abs(d))
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (2, 3) and _squarefree(abs(m))
    return False


def _squarefree(n: int) -> bool:
    return all(e == 1 for _, e in factor(n))


def quadratic(d: int) -> Character:
    if not is_fundamental(d):
        raise ValueError(f"{d} is not a fundamental discriminant")
    q = abs(d)
    num = np.array([{1: 0, -1: 1, 0: -1}[kronecker(d, m)] for m in range(q)], dtype=np.int64)
    return Character(q, num, 2)


def class_number(d: int) -> int:
    if d >= 0 or d % 4 not in (0, 1):
        raise ValueError("class_number counts forms of negative discriminant")
    h = 0
    a = 1
    while 3 * a * a <= -d:
        for b in range(-a + 1, a + 1):
            if (b * b - d) % (4 * a):
                continue
            c = (b * b - d) // (4 * a)
            if c < a or (c == a and b < 0) or math.gcd(math.gcd(a, b), c) != 1:
                continue
            h += 1
        a += 1
    return h


def unit_count(d: int) -> int:
    return {-3: 6, -4: 4}.get(d, 2)


def fundamental_unit(d: int) -> tuple[int, int]:
    if d <= 1 or not is_fundamental(d):
        raise ValueError("need a positive fundamental discriminant")
    if d % 4 == 1:
        P, Q, D = 1, 2, d
    else:
        P, Q, D = 0, 1, d // 4
    root = math.isqrt(D)
    h0, h1, k0, k1 = 0, 1, 1, 0
    while True:
        a = (P + root) // Q
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        if abs(_norm(h1, -k1, d)) == 1:
            return (2 * h1 - k1, k1) if d % 4 == 1 else (2 * h1, k1)
        P = a * Q - P
        Q = (D - P * P) // Q


def _norm(x: int, y: int, d: int) -> int:
    if d % 4 == 1:
        return x * x + x * y - y * y * (d - 1) // 4
    return x * x - y * y * (d // 4)
