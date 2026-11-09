from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import (
    class_number,
    fundamental_unit,
    is_fundamental,
    kronecker,
    quadratic,
    reduced_forms,
    unit_count,
)
from zeta_explorer.lfunction import dirichlet_l, l_theta, l_zeros
from zeta_explorer.riemann_siegel import theta
from zeta_explorer.spacing import _gap_interpolant, first_zeros
from zeta_explorer.zeta import zeta

Form = tuple[int, int, int]


def _check(D: int) -> None:
    if not is_fundamental(D):
        raise ValueError(f"{D} is not a fundamental discriminant")


def ideal_counts(D: int, N: int) -> np.ndarray:
    _check(D)
    chi = np.array([kronecker(D, m) for m in range(abs(D))], dtype=np.int64)
    a = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1):
        a[d::d] += chi[d % abs(D)]
    return a


def dedekind_zeta(s: ArrayLike, D: int) -> complex | np.ndarray:
    _check(D)
    return zeta(s) * dirichlet_l(s, quadratic(D))


def form_counts(D: int, N: int) -> np.ndarray:
    if D >= 0:
        raise ValueError("definite forms need D < 0")
    _check(D)
    r = np.zeros(N + 1, dtype=np.int64)
    for a, b, c in reduced_forms(D):
        ymax = math.isqrt(4 * a * N // -D) + 1
        for y in range(-ymax, ymax + 1):
            disc = y * y * D + 4 * a * N
            if disc < 0:
                continue
            w = math.isqrt(disc)
            for x in range((-b * y - w) // (2 * a) - 1, (-b * y + w) // (2 * a) + 2):
                v = a * x * x + b * x * y + c * y * y
                if v <= N:
                    r[v] += 1
    r[0] = 0
    return r


def reduced_indefinite(D: int) -> list[Form]:
    if D <= 0:
        raise ValueError("indefinite forms need D > 0")
    _check(D)
    out = []
    for b in range(D % 2 or 2, math.isqrt(D) + 1, 2):
        m = (D - b * b) // 4
        for a in range(1, m + 1):
            if m % a or not _is_reduced(a, b, D):
                continue
            out += [(a, b, -m // a), (-a, b, m // a)]
    return out


def _is_reduced(a: int, b: int, D: int) -> bool:
    lo = D < (2 * a + b) ** 2
    hi = 2 * a - b < 0 or (2 * a - b) ** 2 < D
    return 0 < b and b * b < D and lo and hi


def rho(f: Form, D: int) -> Form:
    a, b, c = f
    m = 2 * abs(c)
    r = math.isqrt(D)
    b2 = r - (r + b) % m
    return c, b2, (b2 * b2 - D) // (4 * c)


def cycles(D: int) -> list[list[Form]]:
    left = set(reduced_indefinite(D))
    out = []
    while left:
        start = min(left)
        cyc, f = [], start
        while True:
            cyc.append(f)
            left.discard(f)
            f = rho(f, D)
            if f == start:
                break
            if f not in left:
                raise RuntimeError(f"rho left the reduced forms at {f}")
        out.append(cyc)
    return out


def unit_norm(D: int) -> int:
    x, y = fundamental_unit(D)
    return (x * x - D * y * y) // 4


def regulator(D: int) -> float:
    x, y = fundamental_unit(D)
    return math.log(x) + math.log1p(math.sqrt(D) * y / x) - math.log(2)


def narrow_class_number(D: int) -> int:
    return len(cycles(D))


def real_class_number(D: int) -> int:
    hp = narrow_class_number(D)
    return hp if unit_norm(D) == -1 else hp // 2


def field_class_number(D: int) -> int:
    return class_number(D) if D < 0 else real_class_number(D)


def analytic_class_number(D: int) -> float:
    L1 = dirichlet_l(1.0, quadratic(D)).real
    if D < 0:
        return unit_count(D) * math.sqrt(-D) * L1 / (2 * math.pi)
    return math.sqrt(D) * L1 / (2 * regulator(D))


def zeta_minus_one_siegel(D: int) -> Fraction:
    if D <= 0:
        raise ValueError("Siegel's formula is for real fields")
    _check(D)
    total = 0
    for b in range(-math.isqrt(D), math.isqrt(D) + 1):
        if (D - b * b) % 4 == 0:
            n = (D - b * b) // 4
            total += sum(d + n // d for d in range(1, math.isqrt(n) + 1) if n % d == 0)
            r = math.isqrt(n)
            if r * r == n:
                total -= r
    return Fraction(total, 60)


def zeta_minus_one_bernoulli(D: int) -> Fraction:
    _check(D)
    q = abs(D)
    b2 = sum(kronecker(D, a) * (Fraction(a, q) ** 2 - Fraction(a, q) + Fraction(1, 6))
             for a in range(1, q + 1))
    return q * b2 / 24


def zeros(D: int, T: float, per_gap: int = 32) -> tuple[np.ndarray, np.ndarray]:
    _check(D)
    n = math.ceil(theta(T) / math.pi + 1)
    gz = first_zeros(n + 8)
    gl = l_zeros(quadratic(D), 0.0, T, per_gap)
    return gz[gz < T], gl[gl > 0]


def unfold(D: int, gamma: ArrayLike) -> np.ndarray:
    g = np.asarray(gamma, dtype=float)
    return (np.asarray(theta(g)) + np.asarray(l_theta(g, quadratic(D)))) / math.pi


def superposition_gap(s: ArrayLike) -> np.ndarray:
    E = _gap_interpolant()
    return E(np.asarray(s, dtype=float) / 2) ** 2


def superposition_spacing(s: ArrayLike) -> np.ndarray:
    E = _gap_interpolant()
    h = np.asarray(s, dtype=float) / 2
    return (E.deriv(1)(h) ** 2 + E(h) * E.deriv(2)(h)) / 2


def superposition_cdf(s: ArrayLike) -> np.ndarray:
    E = _gap_interpolant()
    h = np.asarray(s, dtype=float) / 2
    return 1 + E(h) * E.deriv(1)(h)
