from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.afe import _legendre
from zeta_explorer.characters import factor, is_fundamental, quadratic
from zeta_explorer.explicit import primes

_CUT = 36.0


@dataclass(frozen=True)
class Curve:
    a1: int
    a2: int
    a3: int
    a4: int
    a6: int

    @property
    def b(self) -> tuple[int, int, int, int]:
        a1, a2, a3, a4, a6 = self.a1, self.a2, self.a3, self.a4, self.a6
        b2 = a1 * a1 + 4 * a2
        b4 = 2 * a4 + a1 * a3
        b6 = a3 * a3 + 4 * a6
        b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
        return b2, b4, b6, b8

    @property
    def c4(self) -> int:
        b2, b4, _, _ = self.b
        return b2 * b2 - 24 * b4

    @property
    def discriminant(self) -> int:
        b2, b4, b6, b8 = self.b
        return -b2 * b2 * b8 - 8 * b4**3 - 27 * b6 * b6 + 9 * b2 * b4 * b6

    @cached_property
    def bad_primes(self) -> list[int]:
        if self.discriminant == 0:
            raise ValueError("singular curve")
        bad = [p for p, _ in factor(abs(self.discriminant))]
        additive = [p for p in bad if self.c4 % p == 0]
        if additive:
            raise ValueError(f"additive reduction at {additive}, only semistable curves")
        return bad

    @property
    def conductor(self) -> int:
        return math.prod(self.bad_primes)

    def ap(self, p: int) -> int:
        if p == 2:
            pts = sum(
                (y * y + self.a1 * x * y + self.a3 * y - x**3 - self.a2 * x * x - self.a4 * x
                 - self.a6) % 2 == 0
                for x in range(2) for y in range(2)
            )
            return 2 - pts
        b2, b4, b6, _ = self.b
        x = np.arange(p, dtype=np.int64)
        f = ((((4 * x + b2) % p) * x + 2 * b4) % p * x + b6) % p
        chi = np.full(p, -1, dtype=np.int64)
        chi[x * x % p] = 1
        chi[0] = 0
        return -int(chi[f].sum())

    @cached_property
    def root_number(self) -> int:
        return -math.prod(-self.ap(p) for p in self.bad_primes)


@functools.cache
def _coefficients(E: Curve, M: int) -> np.ndarray:
    a = np.ones(M + 1, dtype=np.int64)
    a[0] = 0
    bad = set(E.bad_primes)
    for p in primes(M).tolist():
        ap = E.ap(p)
        powers = [1, ap]
        while p ** len(powers) <= M:
            prev = powers[-1] * ap if p in bad else ap * powers[-1] - p * powers[-2]
            powers.append(prev)
        m = np.arange(p, M + 1, p)
        v = np.ones(len(m), dtype=np.int64)
        r = m // p
        while np.any(hit := r % p == 0):
            v[hit] += 1
            r[hit] //= p
        a[m] *= np.array(powers, dtype=np.int64)[v]
    return a


@dataclass(frozen=True)
class EllipticL:
    E: Curve
    d: int = 1

    def __post_init__(self):
        if self.d != 1 and not is_fundamental(self.d):
            raise ValueError(f"{self.d} is not a fundamental discriminant")
        if math.gcd(self.d, self.E.conductor) != 1:
            raise ValueError("the twist needs gcd(d, N) = 1")

    @property
    def conductor(self) -> int:
        return self.E.conductor * self.d * self.d

    @cached_property
    def root_number(self) -> int:
        w = self.E.root_number
        if self.d == 1:
            return w
        return w * int(self._chi[-self.E.conductor % abs(self.d)])

    @property
    def A(self) -> float:
        return math.sqrt(self.conductor) / (2 * math.pi)

    def an(self, M: int) -> np.ndarray:
        a = _coefficients(self.E, 1 << max(10, (M - 1).bit_length()))[: M + 1]
        if self.d == 1:
            return a
        return a * self._chi[np.arange(M + 1) % abs(self.d)]

    @cached_property
    def _chi(self) -> np.ndarray:
        return quadratic(self.d).values.real.round().astype(np.int64)

    def __repr__(self) -> str:
        e = (self.E.a1, self.E.a2, self.E.a3, self.E.a4, self.E.a6)
        return f"EllipticL({list(e)}, d={self.d})"


def exp1(x: ArrayLike) -> float | np.ndarray:
    xs = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xs <= 0):
        raise ValueError("exp1 needs x > 0")
    out = np.empty_like(xs)
    small = xs < 1
    xx = xs[small]
    term = -xx
    total = term.copy()
    for k in range(2, 40):
        term = -term * xx / k
        total += term / k
    out[small] = -np.euler_gamma - np.log(xx) - total
    w = xs[~small].astype(complex)
    out[~small] = (np.exp(-w) * _legendre(0.0, w)).real
    return float(out[0]) if np.ndim(x) == 0 else out


def central(L: EllipticL) -> float:
    M = math.ceil(_CUT * L.A) + 1
    n = np.arange(1, M + 1)
    a = L.an(M)[1:] / n
    if L.root_number == 1:
        return float(2 * np.sum(a * np.exp(-n / L.A)))
    return float(2 * np.sum(a * exp1(n / L.A)))


def twists(E: Curve, lo: int, hi: int, sign: int = 0) -> list[EllipticL]:
    N = E.conductor
    ds = [d for m in range(lo, hi + 1) for d in (m, -m)
          if (sign == 0 or d * sign > 0) and is_fundamental(d) and math.gcd(d, N) == 1]
    return [EllipticL(E, d) for d in ds]
