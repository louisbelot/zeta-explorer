from __future__ import annotations

import functools
import math
from fractions import Fraction

import mpmath
import numpy as np

from zeta_explorer.delta import lam, root_number, tau
from zeta_explorer.explicit import primes
from zeta_explorer.selfdual import SelfDual, multiplicative

WEIGHT = 12


def gamma_shifts(m: int, k: int = WEIGHT) -> tuple[float, ...]:
    mu = []
    for j in range((m + 1) // 2):
        kappa = (m - 2 * j) * (k - 1) / 2
        mu += [kappa, kappa + 1]
    if m % 2 == 0:
        mu.append(float((m // 2) * (k - 1) % 2))
    return tuple(sorted(mu))


def _angles(P: int) -> tuple[np.ndarray, np.ndarray]:
    p = primes(P)
    return p, np.arccos(np.clip(lam(P)[p] / 2, -1, 1))


def _local(m: int, angle: float, K: int) -> np.ndarray:
    out = np.zeros(K + 1, dtype=complex)
    out[0] = 1
    for j in range(m + 1):
        root = np.exp(1j * (m - 2 * j) * angle)
        for i in range(1, K + 1):
            out[i] += root * out[i - 1]
    return out.real


@functools.cache
def coefficients(m: int, M: int) -> np.ndarray:
    angle = dict(zip(*(x.tolist() for x in _angles(M)), strict=True))
    return multiplicative(M, lambda p, K: _local(m, angle[p], K))


def delta_power(m: int) -> SelfDual:
    if not 1 <= m <= 4:
        raise ValueError("sym^m of Delta is known to be automorphic for m <= 4")
    return SelfDual(lambda M: coefficients(m, M), 1.0, gamma_shifts(m), root_number(m))


def euler(m: int, s: complex, P: int = 1 << 15) -> complex:
    p, angle = _angles(P)
    x = p.astype(float) ** -complex(s)
    j = np.arange(m + 1)[:, None]
    return complex(np.exp(-np.sum(np.log(1 - np.exp(1j * (m - 2 * j) * angle) * x))))


@functools.cache
def petersson(dps: int = 30) -> mpmath.mpf:
    with mpmath.workdps(dps + 5):
        top = 4 + math.ceil(dps * math.log(10) / (math.pi * math.sqrt(3)))
        t = [int(v) for v in tau(top)]
        two_pi = 2 * mpmath.pi

        def strip(x):
            y = mpmath.sqrt(1 - x * x)
            total = mpmath.mpf(0)
            for m in range(1, top):
                for n in range(m, top + 1 - m):
                    c = two_pi * (m + n)
                    weight = (1 if m == n else 2) * t[m] * t[n]
                    upper = mpmath.gammainc(WEIGHT - 1, c * y) / c ** (WEIGHT - 1)
                    total += weight * mpmath.cos(two_pi * (m - n) * x) * upper
            return total

        return +(2 * mpmath.quad(strip, [0, mpmath.mpf(1) / 2]))


ZAGIER = {
    1: Fraction(2**15, 3**4 * 5**2 * 7 * 11),
    3: Fraction(2**17, 3**5 * 5**2 * 7**2 * 11 * 13),
    5: Fraction(2**16, 3**7 * 5**3 * 7**2 * 11 * 13),
    7: Fraction(2**19, 3**9 * 5**4 * 7**2 * 11 * 13 * 17),
    9: Fraction(2**23, 3**8 * 5**5 * 7**4 * 11 * 13 * 17 * 19),
    11: Fraction(2**24, 3**9 * 5**4 * 7**2 * 11 * 13 * 17 * 19 * 23 * 691),
}


def critical_ratio(s: int, P: int = 1 << 16, dps: int = 30) -> mpmath.mpf:
    if s not in ZAGIER:
        raise ValueError("critical points right of the center are s = 1, 3, ..., 11")
    with mpmath.workdps(dps + 5):
        if s < 5:
            value = mpmath.mpf(delta_power(2)(float(s)).real)
        else:
            t = tau(P)
            value = mpmath.mpf(1)
            for p in primes(P).tolist():
                x = mpmath.mpf(p) ** -s
                trace = mpmath.mpf(int(t[p])) ** 2 / mpmath.mpf(p) ** (WEIGHT - 1) - 2
                value /= (1 - x) * (1 - trace * x + x * x)
        return +(value / (mpmath.pi ** (2 * s + WEIGHT - 1) * petersson(dps)))


def rational(x: float | mpmath.mpf, rtol: float) -> Fraction:
    if isinstance(x, mpmath.mpf):
        sign, man, exp, _ = x._mpf_
        f = (-1) ** sign * Fraction(int(man)) * Fraction(2) ** exp
    else:
        f = Fraction(x)
    h0, h1, k0, k1 = 0, 1, 1, 0
    y = f
    while True:
        a = y.numerator // y.denominator
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        if abs(Fraction(h1, k1) - f) <= rtol * abs(f):
            return Fraction(h1, k1)
        y = 1 / (y - a)
